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


class Splitter(Enum):
    def equal_splits(
        circuit: Circuit, splits: int
    ) -> List[Circuit]:
        """
        
        """ 
        return [circuit]

    def clifford_splits(
        circuit: Circuit, splits: int
    ) -> List[Circuit]:
        """
        
        """ 
        return [circuit]



    
#SPLITTING-------------------------------------------------------

    def depth_structure(self) -> list:
        """
        Converts a pytket circuit to a list containing 'x' lists, each containing
        'y' gates, where 'x' is the depth of the circuit and 'y' is the number
        of gates acting in a given timestep. Essentially we split the circuit
        into a list of 'timeslices'.
        
        :returns:       A list containing lists of gates.
        :rtype:         list
        """
        gates = self.c.get_commands()
        qubits = self.c.qubits
        depth = self.c.depth()
        qn = self.c.n_qubits
        current_frontiers = [0]*qn
        depth_slices = [[] for d in range(depth)]
        for gate in gates:
            involved_qubits = gate.qubits
            qubit_indices = []
            for qubit in involved_qubits:
                qubit_indices.append(qubits.index(qubit))
            max_frontier = max([current_frontiers[qid] for qid in qubit_indices])
            for qid in qubit_indices:
                current_frontiers[qid] = max_frontier + 1
            depth_slices[max_frontier].append(gate)
        return depth_slices

    def multi_conversion(self, n: int = 1, strategy: str = "Gates", ) -> list:
        #Currently 'strategy' takes a 'str' type parameter that is either "Depth"
        #or "Gates". Might want to consider additional strategies and switch to
        #an enum in the future.
        """
        Splits a pytket circuit into 'n' subcircuits, each subcircuit containing
        either an approximately equal depth or an approximately equal number of
        non-Clifford gates. Then converts each subcircuit to a measurement pattern,
        extracts a new circuit from the measurement pattern, and ultimately
        returns a list of tuples containing the new circuits and the dictionaries
        mapping the inputs and outputs of the new circuits to the original.
        There is also a third splitting strategy which only converts the Clifford
        parts of a circuit to MBQC.
        
        :param n:        Number of segments to attempt to split into (May return fewer).
        :param type:     int
        
        :param strategy: Splitting strategy either by "Depth", by "Gates" or only convert "Clifford" segments.
        :param type:     str
        
        :returns:        A list of tuples containing circuits and i/o maps.
        :rtype:          list
        """
        depth_structure = self.depth_structure()
        size = len(depth_structure)
        if strategy == "Gates":
            non_cliff = 0
            for d in depth_structure:
                for gate in d:
                    if not MPattern.is_Clifford(gate):
                        non_cliff += 1
            size = non_cliff
        slice_size = math.ceil(size/n)
        done_depth = 0
        output = []
        if strategy == "Depth":
            for curr in range(n):
                finish_at = min(done_depth + slice_size,size)
                subcircuit = Circuit()
                for qubit in self.c.qubits:
                    subcircuit.add_qubit(qubit)
                for bit in self.c.bits:
                    subcircuit.add_bit(bit)
                for depth_list in depth_structure[done_depth:finish_at]:
                    for gate in depth_list:
                        subcircuit.add_gate(Op=gate.op, args=gate.args)
                sub_pattern = MPattern(subcircuit)
                output.append(sub_pattern.single_conversion())
                if finish_at >= size:
                    break
                else:
                    done_depth = finish_at
        elif strategy == "Gates":
            for curr in range(n):
                ncliff_total = 0
                added_depths = 0
                stop_at_next_nClifford = False
                stopped = False
                for depth in depth_structure[done_depth:]:
                    for gate in depth:
                        if not MPattern.is_Clifford(gate):
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
                for qubit in self.c.qubits:
                    subcircuit.add_qubit(qubit)
                for bit in self.c.bits:
                    subcircuit.add_bit(bit)
                for depth_list in depth_structure[done_depth:done_depth+added_depths]:
                    for gate in depth_list:
                        subcircuit.add_gate(Op=gate.op, args=gate.args)
                sub_pattern = MPattern(subcircuit)
                output.append(sub_pattern.single_conversion())
                if done_depth+added_depths >= len(depth_structure):
                    break
                else:
                    done_depth += added_depths
        elif strategy == "Clifford":
            if size > 0:
                clifford_circ = True
                for gate in depth_structure[0]:
                    if not MPattern.is_Clifford(gate):
                        clifford_circ = False
                        break
                while done_depth < size:
                    subcircuit = Circuit()
                    for qubit in self.c.qubits:
                        subcircuit.add_qubit(qubit)
                    for bit in self.c.bits:
                        subcircuit.add_bit(bit)
                    for depth_list in depth_structure[done_depth:size]:
                        has_non_clifford = False
                        for gate in depth_list:
                            if not MPattern.is_Clifford(gate):
                                has_non_clifford = True
                                break
                        if has_non_clifford == clifford_circ:
                            if clifford_circ:
                                sub_pattern = MPattern(subcircuit)
                                output.append(sub_pattern.single_conversion())
                            else:
                                new_map = {"i":{},"o":{}}
                                for qubit in subcircuit.qubits:
                                    new_map["i"][qubit] = subcircuit.qubits.index(qubit)
                                    new_map["o"][qubit] = subcircuit.qubits.index(qubit)
                                new_tuple = (subcircuit.copy(),new_map)
                                output.append(new_tuple)
                            clifford_circ = not clifford_circ
                            break
                        else:
                            done_depth += 1
                            for gate in depth_list:
                                subcircuit.add_gate(Op=gate.op, args=gate.args)
                if clifford_circ:
                    sub_pattern = MPattern(subcircuit)
                    output.append(sub_pattern.single_conversion())
                else:
                    new_map = {"i":{},"o":{}}
                    for qubit in subcircuit.qubits:
                        new_map["i"][qubit] = subcircuit.qubits.index(qubit)
                        new_map["o"][qubit] = subcircuit.qubits.index(qubit)
                    new_tuple = (subcircuit.copy(),new_map)
                    output.append(new_tuple)
        return output



