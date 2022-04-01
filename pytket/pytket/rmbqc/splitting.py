# Copyright 2019-2022 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from enum import Enum
from pytket.circuit import Circuit  # type: ignore
from utils import is_mbqc_clifford, depth_structure
from math import ceil
from typing import Tuple, List


class Splitter(Enum):
    def depth_split(c: Circuit, splits: int) -> List[Tuple[Circuit, bool]]:
        """
        Splits a pytket circuit into a given number of subcircuits of approximately
        equal depth and then returns a list of tuples of pytket circuits and booleans
        indicating whether the circuits should be converted to MBQC.
        
        :param c:        Original circuit to split.
        :param type:     Circuit
        
        :splits:         Maximum number of subcircuits to split into.
        :param type:     int
        
        :returns:        A list of tuples of circuits and booleans indicating whether to convert them or not.
        :rtype:          List[Tuple[Circuit, bool]]
        """
        depth_struct = depth_structure(c)
        depth = len(depth_struct)
        slice_size = ceil(depth / splits)
        done_depth = 0
        subc_list = []
        finish_at = 0
        while finish_at < depth:
            finish_at = min(done_depth + slice_size, depth)
            subcircuit = Circuit()
            for qubit in c.qubits:
                subcircuit.add_qubit(qubit)
            for bit in c.bits:
                subcircuit.add_bit(bit)
            for depth_list in depth_struct[done_depth:finish_at]:
                for gate in depth_list:
                    subcircuit.add_gate(Op=gate.op, args=gate.args)
            subc_list.append((subcircuit, True))
            done_depth = finish_at
        return subc_list

    def gates_split(c: Circuit, splits: int) -> List[Tuple[Circuit, bool]]:
        """
        Splits a pytket circuit into a given number of subcircuits containing an
        approximately equal number of non-Clifford gates and then returns a
        list of tuples of pytket circuits and booleans indicating whether the
        subcircuits should be converted to MBQC or not.
        
        :param c:        Original circuit to split.
        :param type:     Circuit
        
        :splits:         Maximum number of subcircuits to split into.
        :param type:     int
        
        :returns:        A list of tuples of circuits and booleans indicating whether to convert them or not.
        :rtype:          List[Tuple[Circuit, bool]]
        """
        depth_struct = depth_structure(c)  # type: ignore
        non_cliff = 0
        for d in depth_struct:
            for gate in d:
                if not is_mbqc_clifford(gate):
                    non_cliff += 1
        slice_size = ceil(non_cliff / splits)
        done_depth = 0
        subc_list = []
        for curr in range(splits):
            ncliff_total = 0
            added_depths = 0
            stop_at_next_nClifford = False  # This will be flagged true when the current 'slice' contains enough non-Clifford gates
            stopped = False
            # The following loop checks through the remaining timesteps and counts non-Clifford gates.
            # Once the number of non-Clifford gates exceeds 'slice_size', 'stop_at_next_nClifford' is flagged True.
            # However, this doesn't mean the subcircuit is over, because it is still possible to include more timesteps
            # as long as they only contain Clifford gates. So only if stop_at_next_nClifford is True AND we encounter
            # another non-Clifford gate do we flag 'stopped' and exit the loop. As long as 'stopped' is False,
            # we continue moving to the next timeslice.
            for depth in depth_struct[
                done_depth:
            ]:  # Look in all the remaining timesteps of the circuit
                for gate in depth:  # Look through each gate in each such timestep
                    if not is_mbqc_clifford(gate):
                        if (
                            stop_at_next_nClifford
                        ):  # Will only be True if current subcircuit is already saturated with non-Clifford gates.
                            stopped = True
                            break
                        else:
                            ncliff_total += 1  # Keep incrementing for every non-Clifford gate we find.
                if stopped:
                    break
                else:  # stopped is 'False' so the current timeslice can be added to the subcircuit.
                    added_depths += 1
                    if (
                        ncliff_total >= slice_size
                    ):  # If the number of non-Clifford gates exceeds the saturation point we flag this to stop next time we find one.
                        stop_at_next_nClifford = True
            subcircuit = Circuit()
            for qubit in c.qubits:
                subcircuit.add_qubit(qubit)
            for bit in c.bits:
                subcircuit.add_bit(bit)
            # The loop above was only tracking which timeslices are gonna be added to each subcircuit.
            # This is the loop where the gates are actually added to the subcircuit. The starting point
            #'done_depth' is where the previous subcircuit left off. The endpoint 'done_depth + added_depths'
            # is where the current subcircuit will stop.
            for depth_list in depth_struct[done_depth : done_depth + added_depths]:
                for gate in depth_list:
                    subcircuit.add_gate(Op=gate.op, args=gate.args)
            subc_list.append((subcircuit, True))
            if done_depth + added_depths >= len(depth_struct):
                break
            done_depth += added_depths
        return subc_list

    def clifford_split(c: Circuit, splits: int = 1) -> List[Tuple[Circuit, bool]]:
        """
        Splits a pytket circuit into a series of alternating Clifford/non-Clifford
        subcircuits. Returns a list of tuples of Circuits and booleans indicating
        whether the circuits should be converted to MBQC or not.
        
        :param c:        Original circuit to split.
        :param type:     Circuit
        
        :returns:        A list of tuples of circuits and booleans indicating whether to convert them or not.
        :rtype:          List[Tuple[Circuit, bool]]
        """
        depth_struct = depth_structure(c)  # type: ignore
        depth = len(depth_struct)
        done_depth = 0
        subc_list = []
        if depth > 0:
            clifford_circ = True
            for gate in depth_struct[
                0
            ]:  # This loop checks if the first timeslice is a Clifford circuit (i.e. doesn't contain non-Clifford gate)
                if not is_mbqc_clifford(gate):
                    clifford_circ = False
                    break  # If it contains at least one non-Clifford gate it's not Clifford, no need to keep checking.
            while done_depth < depth:  # Look through the remainder of the circuit.
                subcircuit = Circuit()  # Define a new subcircuit.
                for qubit in c.qubits:
                    subcircuit.add_qubit(qubit)
                for bit in c.bits:
                    subcircuit.add_bit(bit)
                for depth_list in depth_struct[done_depth:depth]:
                    has_non_clifford = False
                    for (
                        gate
                    ) in (
                        depth_list
                    ):  # Check if the current timeslice is Clifford or not, like we did above.
                        if not is_mbqc_clifford(gate):
                            has_non_clifford = True
                            break
                    if (
                        has_non_clifford == clifford_circ
                    ):  # If the current timeslice differs from the previous one, then that is the end of the current subcircuit.
                        if (
                            clifford_circ
                        ):  # If the previous subcircuit that has just been concluded is Clifford, it will be flagged for conversion to MBQC.
                            subc_list.append((subcircuit, True))
                        else:  # If the previous subcircuit is non-Clifford, then it will not be flagged for conversion to MBQC.
                            subc_list.append((subcircuit, False))
                        clifford_circ = (
                            not clifford_circ
                        )  # And since the new timeslice doesn't match the previous one, we are switching to the opposite type of subcircuit.
                        break  # ... and exit the current subcircuit to begin a new one of the opposite type.
                    else:  # However, if the new timeslice DOES match the previous one, then we simply include it in the subcircuit and continue iterating.
                        done_depth += 1
                        for gate in depth_list:
                            subcircuit.add_gate(Op=gate.op, args=gate.args)
            if clifford_circ:
                subc_list.append((subcircuit, True))
            else:
                subc_list.append((subcircuit, True))
        return subc_list
