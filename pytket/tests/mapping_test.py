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

from pytket.mapping import MappingManager, RoutingMethodCircuit, LexiRouteRoutingMethod, LexiLabellingMethod  # type: ignore
from pytket.architecture import Architecture  # type: ignore
from pytket import Circuit, OpType
from pytket.circuit import Node, Qubit  # type: ignore
from pytket.placement import Placement  # type: ignore
from typing import Tuple, Dict


# simple deterministic heuristic used for testing purposes
def route_subcircuit_func(
    circuit: Circuit, architecture: Architecture
) -> Tuple[Circuit, Dict[Node, Node], Dict[Node, Node]]:
    #     make a replacement circuit with identical unitds
    replacement_circuit = Circuit()
    for qb in circuit.qubits:
        replacement_circuit.add_qubit(qb)
    for bit in circuit.bits:
        replacement_circuit.add_bit(bit)

    # "place" unassigned logical qubits to physical qubits
    unused_nodes = list(architecture.nodes)
    relabelling_map = dict()

    for qb in circuit.qubits:
        if qb in unused_nodes:
            unused_nodes.remove(qb)

    for qb in circuit.qubits:
        if qb not in architecture.nodes:
            relabelling_map[qb] = unused_nodes.pop()
        else:
            #           this is so later architecture.get_distance works
            #           yes this is obviously bad, buts its a simple test heuristic so who cares?!
            relabelling_map[qb] = qb

    replacement_circuit.rename_units(relabelling_map)
    permutation_map = dict()
    for qb in replacement_circuit.qubits:
        permutation_map[qb] = qb

    #   very simple heuristic -> the first time a physically invalid CX is encountered, add a SWAP
    #   then add all remaining gates as is (using updated physical mapping)
    #   note this is possible as routing accepts partially solved problems
    max_swaps = 1
    swaps_added = 0
    for com in circuit.get_commands():
        rp_qubits = [permutation_map[relabelling_map[q]] for q in com.qubits]
        if len(com.qubits) > 2:
            raise ValueError("Command must have maximum two qubits")
        if len(com.qubits) == 1:
            replacement_circuit.add_gate(com.op.type, rp_qubits)
        if len(com.qubits) == 2:
            if swaps_added < max_swaps:
                for n in architecture.nodes:
                    if n == rp_qubits[0]:
                        n0 = n
                    if n == rp_qubits[1]:
                        n1 = n
                distance = architecture.get_distance(n0, n1)
                if distance > 1:
                    for node in architecture.get_adjacent_nodes(n0):
                        if architecture.get_distance(
                            node, n1
                        ) < architecture.get_distance(n0, n1):
                            replacement_circuit.add_gate(
                                OpType.SWAP, [rp_qubits[0], node]
                            )

                            permutation_map[rp_qubits[0]] = node
                            permutation_map[node] = rp_qubits[0]
                            rp_qubits = [
                                permutation_map[relabelling_map[q]] for q in com.qubits
                            ]
                            swaps_added += 1
                            break

            replacement_circuit.add_gate(com.op.type, rp_qubits)

    return (replacement_circuit, relabelling_map, permutation_map)


def check_subcircuit_func_true(circuit: Circuit, architecture: Architecture) -> bool:
    return True


def check_subcircuit_func_false(circuit: Circuit, architecture: Architecture) -> bool:
    return False


def test_LexiRouteRoutingMethod() -> None:
    test_c = Circuit(3).CX(0, 1).CX(0, 2).CX(1, 2)
    nodes = [Node("test", 0), Node("test", 1), Node("test", 2)]
    test_a = Architecture([[nodes[0], nodes[1]], [nodes[1], nodes[2]]])
    test_mm = MappingManager(test_a)
    test_mm.route_circuit(test_c, [LexiLabellingMethod(), LexiRouteRoutingMethod()])
    routed_commands = test_c.get_commands()

    assert routed_commands[0].op.type == OpType.CX
    assert routed_commands[0].qubits == [nodes[1], nodes[0]]
    assert routed_commands[1].op.type == OpType.CX
    assert routed_commands[1].qubits == [nodes[1], nodes[2]]
    assert routed_commands[2].op.type == OpType.SWAP
    assert routed_commands[2].qubits == [nodes[2], nodes[1]]
    assert routed_commands[3].op.type == OpType.CX
    assert routed_commands[3].qubits == [nodes[0], nodes[1]]


def test_RoutingMethodCircuit_custom() -> None:
    test_c = Circuit(3).CX(0, 1).CX(0, 2).CX(1, 2)
    nodes = [Node("test", 0), Node("test", 1), Node("test", 2)]
    test_a = Architecture([[nodes[0], nodes[1]], [nodes[1], nodes[2]]])

    test_mm = MappingManager(test_a)
    test_mm.route_circuit(
        test_c,
        [RoutingMethodCircuit(route_subcircuit_func, check_subcircuit_func_true, 5, 5)],
    )
    routed_commands = test_c.get_commands()

    assert routed_commands[0].op.type == OpType.CX
    assert routed_commands[0].qubits == [nodes[0], nodes[1]]
    assert routed_commands[1].op.type == OpType.SWAP
    assert routed_commands[1].qubits == [nodes[0], nodes[1]]
    assert routed_commands[2].op.type == OpType.CX
    assert routed_commands[2].qubits == [nodes[1], nodes[2]]
    assert routed_commands[3].op.type == OpType.SWAP
    assert routed_commands[3].qubits == [nodes[0], nodes[1]]
    assert routed_commands[4].op.type == OpType.CX
    assert routed_commands[4].qubits == [nodes[1], nodes[2]]


def test_RoutingMethodCircuit_custom_list() -> None:
    test_c = Circuit(3).CX(0, 1).CX(0, 2).CX(1, 2)
    nodes = [Node("test", 0), Node("test", 1), Node("test", 2)]
    test_a = Architecture([[nodes[0], nodes[1]], [nodes[1], nodes[2]]])

    test_mm = MappingManager(test_a)
    test_mm.route_circuit(
        test_c,
        [
            RoutingMethodCircuit(
                route_subcircuit_func, check_subcircuit_func_false, 5, 5
            ),
            LexiLabellingMethod(),
            LexiRouteRoutingMethod(),
        ],
    )
    routed_commands = test_c.get_commands()

    assert routed_commands[0].op.type == OpType.CX
    assert routed_commands[0].qubits == [nodes[1], nodes[0]]
    assert routed_commands[1].op.type == OpType.CX
    assert routed_commands[1].qubits == [nodes[1], nodes[2]]
    assert routed_commands[2].op.type == OpType.SWAP
    assert routed_commands[2].qubits == [nodes[2], nodes[1]]
    assert routed_commands[3].op.type == OpType.CX
    assert routed_commands[3].qubits == [nodes[0], nodes[1]]

    test_c = Circuit(3).CX(0, 1).CX(0, 2).CX(1, 2)
    test_mm.route_circuit(
        test_c,
        [
            RoutingMethodCircuit(
                route_subcircuit_func, check_subcircuit_func_true, 5, 5
            ),
            LexiLabellingMethod(),
            LexiRouteRoutingMethod(),
        ],
    )
    routed_commands = test_c.get_commands()
    assert routed_commands[0].op.type == OpType.CX
    assert routed_commands[0].qubits == [nodes[0], nodes[1]]
    assert routed_commands[1].op.type == OpType.SWAP
    assert routed_commands[1].qubits == [nodes[0], nodes[1]]
    assert routed_commands[2].op.type == OpType.CX
    assert routed_commands[2].qubits == [nodes[1], nodes[2]]
    assert routed_commands[3].op.type == OpType.SWAP
    assert routed_commands[3].qubits == [nodes[0], nodes[1]]
    assert routed_commands[4].op.type == OpType.CX
    assert routed_commands[4].qubits == [nodes[1], nodes[2]]


def test_basic_mapping() -> None:
    circ = Circuit(5)
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ.CX(0, 1)
    circ.CX(0, 3)
    circ.CX(2, 4)
    circ.CX(1, 4)
    circ.CX(0, 4)

    init_map = dict()
    init_map[Qubit(0)] = Node(0)
    init_map[Qubit(1)] = Node(1)
    init_map[Qubit(2)] = Node(2)
    init_map[Qubit(3)] = Node(3)
    init_map[Qubit(4)] = Node(4)
    pl = Placement(arc)
    pl.place_with_map(circ, init_map)
    MappingManager(arc).route_circuit(circ, [LexiRouteRoutingMethod(50)])
    assert circ.valid_connectivity(arc, False)
    assert len(circ.get_commands()) == 10


if __name__ == "__main__":
    test_LexiRouteRoutingMethod()
    test_RoutingMethodCircuit_custom()
    test_RoutingMethodCircuit_custom_list()
    test_basic_mapping()
