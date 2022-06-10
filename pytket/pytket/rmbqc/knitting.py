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
from pytket.mapping import get_token_swapping_network  # type: ignore
from pytket.circuit import Circuit, Bit, Qubit, OpType  # type: ignore
from pytket.predicates import CompilationUnit  # type: ignore
from pytket.passes import RoutingPass  # type: ignore
from pytket.architecture import Architecture, FullyConnected  # type: ignore
from typing import Dict, Tuple, List, Union


class Knitter(Enum):
    @staticmethod
    def separate(
        pattern_list: List[Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]],
        arch: Architecture,
    ) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
        """
        Takes a list of tuples containing pytket circuits and their corresponding
        input/output maps as well as an architecture. It then routes each circuit on the
        architecture separately. Finally it connects the outputs of one segment to
        the inputs of the next by generating a network of SWAP gates. Then returns
        a tuple containing the new circuit object and its corresponding input/output maps.
        
        :param pattern_list:  List of tuples of circuits with their input/output maps.
                                These map the input/output qubits of the original circuit
                                to their new locations in the qubit register.
        :param type:          List[Tuple[Circuit,Dict[Qubit, int],Dict[Qubit, int]]]
        
        :param arch:          The architecture onto which to fit the circuits.
        :param type:          Architecture
        
        :returns:        A tuple containing a circuit and an input/output map dictionary.
        :rtype:          Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]
        """
        new_c = Circuit()
        for q in arch.nodes:
            new_c.add_qubit(q)
        for p in range(len(pattern_list)):
            (curr_circ, curr_inputs, curr_outputs) = pattern_list[p]
            route_map = {}
            if p == 0:
                for k in curr_inputs.keys():
                    if (
                        k in arch.nodes
                    ):  # If a qubit was already placed on an architecture node in the original circuit, we place it back to the same place.
                        temp_val = curr_inputs[k]
                        route_map[curr_circ.qubits[temp_val]] = k
                        curr_inputs[k] = k
                        if curr_outputs[k] == temp_val:
                            curr_outputs[k] = k
                        else:
                            curr_outputs[k] = curr_circ.qubits[curr_outputs[k]]
                    else:  # Otherwise we retrieve the corresponding qubit object from the register.
                        curr_inputs[k] = curr_circ.qubits[curr_inputs[k]]
                        curr_outputs[k] = curr_circ.qubits[curr_outputs[k]]
                curr_circ.rename_units(
                    route_map
                )  # We perform the partial qubit placement.
            else:  # For subsequent segments we assummed that all the qubits start unplaced and we start fresh since the segments will be linked with a SWAP network anyway.
                for k in curr_inputs.keys():
                    curr_inputs[k] = curr_circ.qubits[curr_inputs[k]]
                    curr_outputs[k] = curr_circ.qubits[curr_outputs[k]]
            cu = CompilationUnit(curr_circ)
            if arch == FullyConnected:
                raise TypeError(
                    "This method doesn't work with a FullyConnected architecture."
                )
            RoutingPass(arch).apply(cu)
            used_nodes = set()
            all_nodes = set(arch.nodes)
            unassigned_qubits = []
            for k in curr_inputs.keys():
                if (
                    cu.initial_map[curr_inputs[k]] in all_nodes
                ):  # We keep track of the architecture nodes which have been used.
                    used_nodes |= {cu.initial_map[curr_inputs[k]]}
                    used_nodes |= {cu.initial_map[curr_outputs[k]]}
                    curr_inputs[k] = cu.initial_map[curr_inputs[k]]
                    curr_outputs[k] = cu.final_map[curr_outputs[k]]
                else:  # We also keep track of all the qubits which have not been placed yet.
                    unassigned_qubits.append(k)
            for command in cu.circuit.get_commands():
                if (
                    command.op.type == OpType.CZ
                ):  # Initially 'used_nodes' only contains the nodes where the inputs were placed, but through the CZ interactions we can find other nodes occupied by logical qubits.
                    for q in command.qubits:
                        if q in used_nodes:
                            used_nodes |= set(command.qubits)
                            break
                elif (
                    command.op.type == OpType.SWAP
                ):  # If a SWAP gate is found between a qubit in used_nodes and a qubit not in used_nodes, then the latter takes the former's place.
                    if (
                        (command.qubits[0] in used_nodes)
                        and not (command.qubits[1] in used_nodes)
                    ) or (
                        not (command.qubits[0] in used_nodes)
                        and (command.qubits[1] in used_nodes)
                    ):
                        used_nodes ^= set(command.qubits)
                elif (
                    command.op.type == OpType.BRIDGE
                ):  # The BRIDGE gate is also a logical interaction between the first and third qubit, so if one is in 'used_nodes' the other must be too.
                    for q in {command.qubits[0], command.qubits[2]}:
                        if q in used_nodes:
                            used_nodes |= {command.qubits[0], command.qubits[2]}
                            break
            permutation = {
                x: x for x in all_nodes
            }  # We make note of the qubit permutations caused by SWAP gates throughout the circuit.
            for com in cu.circuit.commands_of_type(OpType.SWAP):
                permutation[com.qubits[0]] = com.qubits[1]
                permutation[com.qubits[1]] = com.qubits[0]
            unused_nodes = (
                all_nodes - used_nodes
            )  # The unused nodes are where we want to place the unassigned qubits.
            segment_circuit = cu.circuit.copy()
            for (
                uq
            ) in (
                unassigned_qubits
            ):  # Each unassigned qubit is placed on an unused node.
                temp = unused_nodes.pop()
                curr_outputs[
                    uq
                ] = temp  # However, we only know where the qubit will be at the end of the circuit.
                for (
                    k
                ) in (
                    permutation.keys()
                ):  # We must use the noted permutations to trace the previously unassigned node back to the beginning.
                    if permutation[k] == temp:
                        segment_circuit.rename_units(
                            {curr_inputs[uq]: k}
                        )  # So this is where the qubit is going to be placed at the beginning of the circuit.
                        curr_inputs[uq] = k
                        break
            new_tuple = (segment_circuit, curr_inputs, curr_outputs)
            pattern_list[p] = new_tuple
            for q in arch.nodes:
                if not q in segment_circuit.qubits:
                    segment_circuit.add_qubit(q)
            if (
                p > 0
            ):  # If this isn't the first segment then we need to use token swapping to link to the previous segment.
                prev_outputs = pattern_list[p - 1][2]
                matching_dict = {}
                for k in curr_inputs.keys():
                    matching_dict[prev_outputs[k]] = curr_inputs[k]
                swaps_as_pairs = get_token_swapping_network(arch, matching_dict)
                for pair in swaps_as_pairs:
                    new_c.SWAP(qubit_0=pair[0], qubit_1=pair[1])
            prev_bits = len(new_c.bits)
            b_map = {}
            for b in range(len(segment_circuit.bits)):
                b_map[segment_circuit.bits[b]] = Bit(prev_bits + b)
            segment_circuit.rename_units(b_map)
            new_c.add_circuit(segment_circuit, [])
        return (new_c, pattern_list[0][1], pattern_list[-1][2])

    @staticmethod
    def sequential(
        pattern_list: List[Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]],
        arch: Architecture,
    ) -> Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]:
        """
        Takes a list of tuples containing pytket circuits and their corresponding
        input/output maps as well as an architecture. It then routes the first circuit onto
        the architecture. For every subsequent circuit, the inputs are forcefully
        assigned to the location of the corresponding outputs of the previous segment
        on the architecture. The current subcircuit is then routed as best as possible
        given those initial conditions. One all subcircuits have been placed, the
        method returns a tuple containing the new circuit object and its corresponding
        input/output maps.
        
        :param pattern_list:  List of tuples of circuits with their input/output maps.
                                These map the inputs/outputs of the original circuit to
                                their new locations in the qubit register.
        :param type:          List[Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]]
        
        :param arch:          The architecture onto which to fit the circuits.
        :param type:          Architecture
        
        :returns:        A tuple containing a circuit and input/output map dictionaries.
        :rtype:          Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]
        """
        new_c = Circuit()
        for q in arch.nodes:
            new_c.add_qubit(q)
        for p in range(len(pattern_list)):
            (curr_circ, curr_inputs, curr_outputs) = pattern_list[p]
            if p > 0:
                (prev_circ, prev_inputs, prev_outputs) = pattern_list[p - 1]
                route_map = {}
                for k in curr_inputs.keys():
                    if (
                        curr_outputs[k] == curr_inputs[k]
                    ):  # If the input and output qubits in the new circuit which correspond to some qubit in the original circuit are the same, then the current output is set to the corresponding output of the previous segment.
                        curr_outputs[k] = prev_outputs[k]
                    else:  # Otherwise, we retrieve the qubit object from the qubit register.
                        curr_outputs[k] = curr_circ.qubits[curr_outputs[k]]
                    route_map[curr_circ.qubits[curr_inputs[k]]] = prev_outputs[
                        k
                    ]  # The inputs will be placed where the corresponding outputs of the previous segment were placed, so that the two circuits can be matched.
                    curr_inputs[k] = prev_outputs[k]
                curr_circ.rename_units(route_map)  # Partial placement of the qubits.
            else:
                route_map = {}
                for k in curr_inputs.keys():
                    if (
                        k in arch.nodes
                    ):  # If qubit is already placed on an architecture node then it stays there.
                        temp_val = curr_inputs[k]
                        route_map[curr_circ.qubits[temp_val]] = k
                        curr_inputs[
                            k
                        ] = k  # Map the input qubit of the original circuit to the architecture node.
                        if (
                            curr_outputs[k] == temp_val
                        ):  # And if the output corresponding to this input is the same, then it is also mapped to the same architecture node.
                            curr_outputs[k] = k
                        else:  # Otherwise, we retrieve the qubit object from the register to map the output to.
                            curr_outputs[k] = curr_circ.qubits[curr_outputs[k]]
                    else:  # If the input qubit is not an architecture node we retrieve the input and output qubit objects from the register.
                        curr_inputs[k] = curr_circ.qubits[curr_inputs[k]]
                        curr_outputs[k] = curr_circ.qubits[curr_outputs[k]]
                curr_circ.rename_units(
                    route_map
                )  # Partial placement of qubits on the architecture.
            cu = CompilationUnit(curr_circ)
            if arch == FullyConnected:
                raise TypeError(
                    "This method doesn't work with a FullyConnected architecture."
                )
            RoutingPass(arch).apply(cu)
            used_nodes = set()
            all_nodes = set(arch.nodes)
            unassigned_qubits = []
            for k in curr_inputs.keys():
                if (
                    cu.initial_map[curr_inputs[k]] in all_nodes
                ):  # We keep track of the architecture nodes which have been used.
                    used_nodes |= {cu.initial_map[curr_inputs[k]]}
                    used_nodes |= {cu.initial_map[curr_outputs[k]]}
                    curr_inputs[k] = cu.initial_map[curr_inputs[k]]
                    curr_outputs[k] = cu.final_map[curr_outputs[k]]
                else:  # We also keep track of all the qubits which have not been placed yet.
                    unassigned_qubits.append(k)
            for command in cu.circuit.get_commands():
                if (
                    command.op.type == OpType.CZ
                ):  # Initially 'used_nodes' only contains the nodes where the inputs were placed, but through the CZ interactions we can find other nodes occupied by logical qubits.
                    for q in command.qubits:
                        if q in used_nodes:
                            used_nodes |= set(command.qubits)
                            break
                elif (
                    command.op.type == OpType.SWAP
                ):  # If a SWAP gate is found between a qubit in used_nodes and a qubit not in used_nodes, then the latter takes the former's place.
                    if (
                        (command.qubits[0] in used_nodes)
                        and not (command.qubits[1] in used_nodes)
                    ) or (
                        not (command.qubits[0] in used_nodes)
                        and (command.qubits[1] in used_nodes)
                    ):
                        used_nodes ^= set(command.qubits)
                elif (
                    command.op.type == OpType.BRIDGE
                ):  # The BRIDGE gate is also a logical interaction between the first and third qubit, so if one is in 'used_nodes' the other must be too.
                    for q in {command.qubits[0], command.qubits[2]}:
                        if q in used_nodes:
                            used_nodes |= {command.qubits[0], command.qubits[2]}
                            break
            permutation = {
                x: x for x in all_nodes
            }  # We make note of the qubit permutations caused by SWAP gates throughout the circuit.
            for com in cu.circuit.commands_of_type(OpType.SWAP):
                permutation[com.qubits[0]] = com.qubits[1]
                permutation[com.qubits[1]] = com.qubits[0]
            unused_nodes = (
                all_nodes - used_nodes
            )  # The unused nodes are where we want to place the unassigned qubits.
            segment_circuit = cu.circuit.copy()
            for (
                uq
            ) in (
                unassigned_qubits
            ):  # Each unassigned qubit is placed on an unused node.
                temp = unused_nodes.pop()
                curr_outputs[
                    uq
                ] = temp  # However, we only know where the qubit will be at the end of the circuit.
                for (
                    k
                ) in (
                    permutation.keys()
                ):  # We must use the noted permutations to trace the previously unassigned node back to the beginning.
                    if permutation[k] == temp:
                        segment_circuit.rename_units(
                            {curr_inputs[uq]: k}
                        )  # So this is where the qubit is going to be placed at the beginning of the circuit.
                        curr_inputs[uq] = k
                        break
            new_tuple = (segment_circuit, curr_inputs, curr_outputs)
            pattern_list[p] = new_tuple
            for q in arch.nodes:
                if not q in segment_circuit.qubits:
                    segment_circuit.add_qubit(q)
            prev_bits = len(new_c.bits)
            b_map = {}
            for b in range(len(segment_circuit.bits)):
                b_map[segment_circuit.bits[b]] = Bit(prev_bits + b)
            segment_circuit.rename_units(b_map)
            new_c.add_circuit(segment_circuit, [])
        return (new_c, pattern_list[0][1], pattern_list[-1][2])

    @staticmethod
    def unrouted(
        pattern_list: List[Tuple[Circuit, Dict[str, Dict[Qubit, int]]]],
        arch: Architecture = FullyConnected,
    ) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
        """
        Takes a list of tuples containing pytket circuits and their corresponding
        input/output maps and joins them together into a new pytket circuit object assuming
        full connectivity. Then returns a tuple containing the new circuit object
        and its corresponding input/output map.
        
        :param pattern_list:  List of tuples of circuits with their input/output maps.
                                These map the inputs/outputs of the original circuit to
                                their new location in the qubit register.
        :param type:          List[Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]]
        
        :returns:        A tuple containing a circuit and an input/output map dictionaries.
        :rtype:          Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]
        """
        new_c = Circuit()
        prev_inputs: Dict[Qubit, Union[Qubit, int]] = {}
        prev_outputs: Dict[Qubit, Union[Qubit, int]] = {}
        for p in range(len(pattern_list)):
            (curr_circ, curr_inputs, curr_outputs) = pattern_list[p]
            if (
                len(prev_inputs) == 0
            ):  # If it's the first segment then it's just added as it is.
                new_c.add_circuit(curr_circ, [])
                prev_inputs = curr_inputs.copy()
                prev_outputs = curr_outputs.copy()
            else:
                q_map = {}
                for (
                    k
                ) in (
                    prev_outputs.keys()
                ):  # Map the current segment's inputs to the previous segment's outputs.
                    q_map[curr_circ.qubits[curr_inputs[k]]] = new_c.qubits[
                        prev_outputs[k]
                    ]
                prev_ancillas = []
                curr_ancillas = []
                for (
                    q
                ) in (
                    new_c.qubits
                ):  # Any qubit that was not an output of the previous circuit was an ancilla.
                    if not q in list(q_map.values()):
                        prev_ancillas.append(q)
                for (
                    q
                ) in (
                    curr_circ.qubits
                ):  # Any qubit that is not an input in the current circuit is also an ancilla.
                    if not q in q_map.keys():
                        curr_ancillas.append(q)
                while (
                    len(prev_ancillas) > 0 and len(curr_ancillas) > 0
                ):  # Arbitrarily assign qubits of the new segment to ancillas in the previous segment.
                    q_map[curr_ancillas.pop()] = prev_ancillas.pop()
                if (
                    len(curr_ancillas) > 0
                ):  # If we still have ancillas in the current segment but the ancillas in the previous segment have ran out, then we need to assign them to new qubit ids.
                    unused_id_pool = []
                    for q in curr_circ.qubits:
                        if q not in list(q_map.values()):
                            unused_id_pool.append(
                                q
                            )  # Keep track of the qubits we already have in the current circuit's register which haven't been mapped to.
                    for q in curr_circ.qubits:
                        if (
                            q not in q_map.keys()
                        ):  # All the qubits in current segment that haven't been assigned yet, are now assigned to one of the locations we kept track of.
                            q_map[q] = unused_id_pool.pop()
                for k in curr_outputs.keys():
                    curr_outputs[k] = q_map[
                        curr_circ.qubits[curr_outputs[k]]
                    ]  # Use the map between segments to update the output map for the current segment with the new assigned indices of the outputs.
                curr_circ.rename_units(
                    q_map
                )  # Use the map we developed to rename the qubits of the current segment to the corresponding qubits of the previous segment.
                prev_bits = len(new_c.bits)
                b_map = {}
                for b in range(len(curr_circ.bits)):
                    b_map[curr_circ.bits[b]] = Bit(prev_bits + b)
                curr_circ.rename_units(b_map)
                new_c.add_circuit(
                    curr_circ, []
                )  # Add the current segment to the rest of the new circuit.
                for k in curr_outputs.keys():
                    curr_outputs[k] = new_c.qubits.index(curr_outputs[k])
                prev_inputs = curr_inputs.copy()
                prev_outputs = curr_outputs.copy()
        for io in [
            pattern_list[0][1],
            prev_outputs,
        ]:  # For each of the input/output maps.
            for q in io.keys():  # For each qubit in that map.
                io[q] = new_c.qubits[
                    io[q]
                ]  # Switch the values of the maps from indices to the actual qubit objects.
        return (new_c, pattern_list[0][1], prev_outputs)
