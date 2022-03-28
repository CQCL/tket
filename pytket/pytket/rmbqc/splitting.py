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
from pytket.circuit import Circuit, Command  # type: ignore
from utils import is_Clifford
from math import ceil
from typing import Tuple, List


class Splitter(Enum):
    def depth_structure(c: Circuit) -> List[List[Command]]:
        """
        Converts a pytket circuit to a list containing 'x' lists, each containing
        'y' gates, where 'x' is the depth of the circuit and 'y' is the number
        of gates acting in a given timestep. Essentially we split the circuit
        into a list of 'timeslices'.
        
        :returns:       A list containing lists of gates.
        :rtype:         List[List[Command]]
        """
        qubits = c.qubits
        current_frontiers = [0] * c.n_qubits
        depth_slices: List[List[int]] = [[] for d in range(c.depth)]
        for gate in c.get_commands():
            involved_qubits = gate.qubits
            qubit_indices = []
            for qubit in involved_qubits:
                qubit_indices.append(qubits.index(qubit))
            max_frontier = max([current_frontiers[qid] for qid in qubit_indices])
            for qid in qubit_indices:
                current_frontiers[qid] = max_frontier + 1
            depth_slices[max_frontier].append(gate)
        return depth_slices
    
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
        depth_structure = Splitter.depth_structure(c)
        depth = len(depth_structure)
        slice_size = ceil(depth / splits)
        done_depth = 0
        subc_list = []
        for curr in range(splits):
            finish_at = min(done_depth + slice_size, depth)
            subcircuit = Circuit()
            for qubit in c.qubits:
                subcircuit.add_qubit(qubit)
            for bit in c.bits:
                subcircuit.add_bit(bit)
            for depth_list in depth_structure[done_depth:finish_at]:
                for gate in depth_list:
                    subcircuit.add_gate(Op=gate.op, args=gate.args)
            new_tuple = (subcircuit, True)
            subc_list.append(new_tuple)
            if finish_at >= depth:
                break
            else:
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
        depth_structure = Splitter.depth_structure()  # type: ignore
        non_cliff = 0
        for d in depth_structure:
            for gate in d:
                if not is_Clifford(gate):
                    non_cliff += 1
        slice_size = ceil(non_cliff / splits)
        done_depth = 0
        subc_list = []
        for curr in range(splits):
            ncliff_total = 0
            added_depths = 0
            stop_at_next_nClifford = False
            stopped = False
            for depth in depth_structure[done_depth:]:
                for gate in depth:
                    if not is_Clifford(gate):
                        if stop_at_next_nClifford:
                            stopped = True
                            break
                        else:
                            ncliff_total += 1
                if stopped:
                    break
                else:
                    added_depths += 1
                    if ncliff_total >= slice_size:
                        stop_at_next_nClifford = True
            subcircuit = Circuit()
            for qubit in c.qubits:
                subcircuit.add_qubit(qubit)
            for bit in c.bits:
                subcircuit.add_bit(bit)
            for depth_list in depth_structure[done_depth : done_depth + added_depths]:
                for gate in depth_list:
                    subcircuit.add_gate(Op=gate.op, args=gate.args)
            new_tuple = (subcircuit, True)
            subc_list.append(new_tuple)
            if done_depth + added_depths >= len(depth_structure):
                break
            else:
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
        depth_structure = Splitter.depth_structure()  # type: ignore
        depth = len(depth_structure)
        done_depth = 0
        subc_list = []
        if depth > 0:
            clifford_circ = True
            for gate in depth_structure[0]:
                if not is_Clifford(gate):
                    clifford_circ = False
                    break
            while done_depth < depth:
                subcircuit = Circuit()
                for qubit in c.qubits:
                    subcircuit.add_qubit(qubit)
                for bit in c.bits:
                    subcircuit.add_bit(bit)
                for depth_list in depth_structure[done_depth:depth]:
                    has_non_clifford = False
                    for gate in depth_list:
                        if not is_Clifford(gate):
                            has_non_clifford = True
                            break
                    if has_non_clifford == clifford_circ:
                        if clifford_circ:
                            new_tuple = (subcircuit, True)
                            subc_list.append(new_tuple)
                        else:
                            new_tuple = (subcircuit, False)
                            subc_list.append(new_tuple)
                            """
                            new_map = {"i":{},"o":{}}
                            for qubit in subcircuit.qubits:
                                new_map["i"][qubit] = subcircuit.qubits.index(qubit)
                                new_map["o"][qubit] = subcircuit.qubits.index(qubit)
                            new_tuple = (subcircuit.copy(),new_map)
                            """
                        clifford_circ = not clifford_circ
                        break
                    else:
                        done_depth += 1
                        for gate in depth_list:
                            subcircuit.add_gate(Op=gate.op, args=gate.args)
            if clifford_circ:
                new_tuple = (subcircuit, True)
                subc_list.append(new_tuple)
            else:
                new_tuple = (subcircuit, False)
                subc_list.append(new_tuple)
                """
                new_map = {"i":{},"o":{}}
                for qubit in subcircuit.qubits:
                    new_map["i"][qubit] = subcircuit.qubits.index(qubit)
                    new_map["o"][qubit] = subcircuit.qubits.index(qubit)
                new_tuple = (subcircuit.copy(),new_map)
                """
        return subc_list