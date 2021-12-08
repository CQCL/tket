# Copyright 2019-2021 Cambridge Quantum Computing
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

from pathlib import Path
from pytket.circuit import OpType, Qubit, Node, Circuit  # type: ignore
from pytket.routing import (  # type: ignore
    NodeGraph,
    Architecture,
    LinePlacement,
    GraphPlacement,
    NoiseAwarePlacement,
    Placement,
    SquareGrid,
    FullyConnected,
    place_with_map,
    route,
)
from pytket.predicates import CompilationUnit, NoMidMeasurePredicate  # type: ignore
from pytket.passes import (  # type: ignore
    DefaultMappingPass,
    FullMappingPass,
    RoutingPass,
    PlacementPass,
    CXMappingPass,
    AASRouting,
    PauliSimp,
    CNotSynthType,
)
from pytket.qasm import circuit_from_qasm
from pytket.transform import Transform  # type: ignore
import numpy as np
import pytest  # type: ignore

import json


def test_architectures() -> None:
    basic_index_coupling = [(0, 1), (2, 1), (2, 3), (4, 3)]
    basic_index_architecture = Architecture(basic_index_coupling)
    basic_index_coupling_convert = [
        (Node(0), Node(1)),
        (Node(2), Node(1)),
        (Node(2), Node(3)),
        (Node(4), Node(3)),
    ]
    assert basic_index_architecture.coupling == basic_index_coupling_convert

    node_0 = Node("example_register", 0)
    node_1 = Node("example_register", 1)
    node_2 = Node("example_register", 2)
    node_3 = Node("example_register", 3)
    basic_uid_coupling = [(node_0, node_1), (node_1, node_2), (node_2, node_3)]
    basic_uid_architecture = Architecture(basic_uid_coupling)
    assert basic_uid_architecture.coupling == basic_uid_coupling

    square_arc = SquareGrid(2, 2, 2)
    assert square_arc.nodes[0] == Node("gridNode", [0, 0, 0])
    assert square_arc.coupling[0] == (
        Node("gridNode", [0, 0, 0]),
        Node("gridNode", [0, 1, 0]),
    )


def test_architecture_eq() -> None:
    coupling = [(1, 2), (3, 4), (0, 6), (0, 3)]
    arc = Architecture(coupling)

    assert arc != Architecture([])
    assert arc == Architecture(coupling)
    assert arc == Architecture([(Node(i), Node(j)) for (i, j) in coupling])
    assert arc != Architecture([(Node("s", i), Node("s", j)) for (i, j) in coupling])

    # only Node IDs and coupling matters
    g00, g01, g10, g11 = [
        Node("gridNode", [i, j, 0]) for i in range(2) for j in range(2)
    ]
    sq_arc = Architecture([(g00, g01), (g01, g11), (g00, g10), (g10, g11)])
    assert sq_arc == SquareGrid(2, 2)
    assert sq_arc != Architecture([(g00, g01), (g01, g11), (g00, g10)])


def test_fully_connected() -> None:
    fc = FullyConnected(3)
    assert fc.nodes == [Node("fcNode", i) for i in range(3)]
    d = fc.to_dict()
    fc1 = FullyConnected.from_dict(d)
    assert fc == fc1


def test_arch_types() -> None:
    arch = Architecture([(0, 1)])
    assert isinstance(arch, Architecture)
    assert isinstance(arch, NodeGraph)
    fc = FullyConnected(2)
    assert isinstance(fc, FullyConnected)
    assert isinstance(fc, NodeGraph)
    sg = SquareGrid(2, 2, 2)
    assert isinstance(sg, SquareGrid)
    assert isinstance(sg, NodeGraph)


def test_placements() -> None:
    test_coupling = [(0, 1), (1, 2), (1, 3), (4, 1), (4, 5)]
    test_architecture = Architecture(test_coupling)
    circ = Circuit(6)
    for pair in test_coupling:
        circ.CX(pair[0], pair[1])
    circ_qbs = circ.qubits
    base_pl = Placement(test_architecture)
    line_pl = LinePlacement(test_architecture)
    graph_pl = GraphPlacement(test_architecture)
    base_placed = circ.copy()
    line_placed = circ.copy()
    graph_placed = circ.copy()

    base_map = base_pl.get_placement_map(circ)
    line_map = line_pl.get_placement_map(circ)
    graph_map = graph_pl.get_placement_map(circ)

    assert base_map != line_map
    assert base_map != graph_map
    assert circ.qubits == circ_qbs

    base_pl.place(base_placed)
    line_pl.place(line_placed)
    graph_pl.place(graph_placed)

    assert line_placed.qubits[0] == line_map[circ_qbs[0]]
    assert line_placed.qubits[1] == line_map[circ_qbs[1]]
    assert line_placed.qubits[2] == line_map[circ_qbs[2]]

    assert base_placed.qubits[0] == base_map[circ_qbs[0]]
    assert base_placed.qubits[1] == base_map[circ_qbs[1]]
    assert base_placed.qubits[2] == base_map[circ_qbs[2]]

    assert graph_placed.qubits[0] == graph_map[circ_qbs[0]]
    assert graph_placed.qubits[1] == graph_map[circ_qbs[1]]
    assert graph_placed.qubits[2] == graph_map[circ_qbs[2]]

    assert circ_qbs != base_placed.qubits
    assert circ_qbs != line_placed.qubits
    assert circ_qbs != graph_placed.qubits

    base_placed = route(base_placed, test_architecture)
    line_placed = route(line_placed, test_architecture)
    graph_placed = route(graph_placed, test_architecture)

    assert base_placed.valid_connectivity(test_architecture, False)
    assert line_placed.valid_connectivity(test_architecture, False)
    assert graph_placed.valid_connectivity(test_architecture, False)


def test_placements_serialization() -> None:
    with open(
        Path(__file__).resolve().parent / "json_test_files" / "placements.json", "r"
    ) as f:
        dict = json.load(f)
        base_pl_serial = dict["base_placement"]
        line_pl_serial = dict["line_placement"]
        graph_pl_serial = dict["graph_placement"]
        noise_pl_serial = dict["noise_placement"]

    assert Placement.from_dict(base_pl_serial).to_dict() == base_pl_serial
    assert LinePlacement.from_dict(line_pl_serial).to_dict() == line_pl_serial
    assert GraphPlacement.from_dict(graph_pl_serial).to_dict() == graph_pl_serial
    assert NoiseAwarePlacement.from_dict(noise_pl_serial).to_dict() == noise_pl_serial


def test_placement_config() -> None:
    test_coupling = [(0, 1), (1, 2), (1, 3), (4, 1), (4, 5)]
    test_architecture = Architecture(test_coupling)
    test_pl = GraphPlacement(test_architecture)
    test_circuit = Circuit(6)
    test_circuit.CX(0, 1)
    test_circuit.CX(2, 3)
    test_circuit.CX(4, 3)
    test_circuit.CX(2, 4)
    test_circuit.CX(3, 5)
    test_circuit.CX(0, 5)
    circ1 = test_circuit.copy()
    circ2 = test_circuit.copy()
    map1 = test_pl.get_placement_map(test_circuit)
    test_pl.place(circ1)
    test_pl.modify_config(
        max_matches=1, depth_limit=0, max_interaction_edges=2, timeout=100
    )
    map2 = test_pl.get_placement_map(test_circuit)
    test_pl.place(circ2)
    assert map1 != map2
    circ1 = route(circ1, test_architecture)
    circ2 = route(circ2, test_architecture)
    assert circ1.n_gates < circ2.n_gates


def test_convert_index_mapping() -> None:
    test_circuit = Circuit(6)
    test_circuit.CX(0, 1)
    test_circuit.CX(2, 3)
    test_circuit.CX(4, 3)
    test_circuit.CX(2, 4)
    test_circuit.CX(3, 5)
    test_circuit.CX(0, 5)

    c0 = test_circuit.copy()
    c1 = test_circuit.copy()

    index_map = {0: 1, 1: 2, 2: 0, 3: 4, 4: 3}
    uid_map = {Qubit(i): Node(j) for i, j in index_map.items()}
    circ_qbs = test_circuit.qubits
    assert uid_map[circ_qbs[0]] == Node(1)
    assert uid_map[circ_qbs[1]] == Node(2)
    assert uid_map[circ_qbs[2]] == Node(0)
    assert uid_map[circ_qbs[3]] == Node(4)
    assert uid_map[circ_qbs[4]] == Node(3)

    place_with_map(test_circuit, uid_map)

    new_circ_qbs = test_circuit.qubits
    assert circ_qbs != new_circ_qbs
    assert new_circ_qbs[0] == Node(0)
    assert new_circ_qbs[1] == Node(1)
    assert new_circ_qbs[2] == Node(2)
    assert new_circ_qbs[3] == Node(3)
    assert new_circ_qbs[4] == Node(4)
    assert new_circ_qbs[5] == Qubit("unplaced", 0)

    index_map_0 = {0: 5, 1: 4, 2: 0, 3: 1, 4: 3, 5: 2}
    index_map_1 = {0: 1, 1: 2, 2: 0, 3: 4, 4: 3, 5: 5}
    uid_0 = {Qubit(i): Node(j) for i, j in index_map_0.items()}
    uid_1 = {Qubit(i): Node(j) for i, j in index_map_1.items()}
    assert uid_0 != uid_1

    place_with_map(c0, uid_0)
    place_with_map(c1, uid_1)
    assert c0 != c1


def test_basic_routing() -> None:
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
    out_circ = route(circ, arc, swap_lookahead=50)
    assert out_circ.valid_connectivity(arc, False)
    assert len(out_circ.get_commands()) == 10


def test_basic_routing_with_line_map() -> None:
    circ = Circuit(5)
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ.CX(0, 1)
    circ.CX(0, 3)
    circ.CX(2, 4)
    circ.CX(1, 4)
    circ.CX(0, 4)
    lp = LinePlacement(arc)
    lp.place(circ)
    out_circ = route(circ, arc)
    assert out_circ.valid_connectivity(arc, False)
    assert len(out_circ.get_commands()) == 6


def test_basic_routing_with_noise_map() -> None:
    circ = Circuit(5)
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ.CX(0, 1)
    circ.CX(0, 3)
    circ.CX(2, 4)
    circ.CX(1, 4)
    circ.CX(0, 4)

    oq_fids = [
        [Node(0), 0.999],
        [Node(1), 0.999],
        [Node(2), 0.999],
        [Node(3), 0.999],
        [Node(4), 0.999],
    ]
    tq_fids = [
        [[Node(0), Node(1)], 0.9],
        [[Node(1), Node(0)], 0.9],
        [[Node(1), Node(2)], 0.89],
        [[Node(2), Node(1)], 0.89],
        [[Node(2), Node(3)], 0.7],
        [[Node(3), Node(2)], 0.7],
        [[Node(3), Node(4)], 0.59],
        [[Node(4), Node(3)], 0.59],
    ]

    tq_errs_dict = {
        (Node(0), Node(1)): 0.1,
        (Node(1), Node(0)): 0.1,
        (Node(1), Node(2)): 0.11,
        (Node(2), Node(1)): 0.11,
        (Node(2), Node(3)): 0.3,
        (Node(3), Node(2)): 0.3,
        (Node(3), Node(4)): 0.41,
        (Node(4), Node(3)): 0.41,
    }
    oq_errs_dict = {node: 1.0 - value for node, value in oq_fids}

    nap = NoiseAwarePlacement(arc, oq_errs_dict, tq_errs_dict)
    nap.place(circ)
    out_circ = route(circ, arc)
    assert len(out_circ.get_commands()) == 6
    assert out_circ.valid_connectivity(arc, False)


def test_greedy_noise_route() -> None:
    circ = Circuit(5)
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ.CX(0, 1)
    circ.CX(0, 3)
    circ.CX(2, 4)
    circ.CX(1, 4)
    circ.CX(0, 4)

    oq_fids = [
        [Node(0), 0.999],
        [Node(1), 0.999],
        [Node(2), 0.999],
        [Node(3), 0.999],
        [Node(4), 0.999],
    ]

    tq_errs_dict = {
        (Node(0), Node(1)): 0.1,
        (Node(1), Node(0)): 0.1,
        (Node(1), Node(2)): 0.11,
        (Node(2), Node(1)): 0.11,
        (Node(2), Node(3)): 0.3,
        (Node(3), Node(2)): 0.3,
        (Node(3), Node(4)): 0.41,
        (Node(4), Node(3)): 0.41,
    }
    oq_errs_dict = {node: 1.0 - value for node, value in oq_fids}
    nap = NoiseAwarePlacement(arc, oq_errs_dict, tq_errs_dict)
    nap.place(circ)
    out_circ = route(circ, arc)

    assert len(out_circ.get_commands()) == 6
    assert out_circ.valid_connectivity(arc, False)


def test_decompose_swap_to_cx() -> None:
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

    out_circ = route(circ, arc)
    assert out_circ.valid_connectivity(arc, False)
    Transform.DecomposeSWAPtoCX(arc).apply(out_circ)
    assert len(out_circ.get_commands()) == 20
    Transform.DecomposeCXDirected(arc).apply(out_circ)
    assert out_circ.valid_connectivity(arc, True)
    assert len(out_circ.get_commands()) == 40


def test_commuting_sq_through_swap() -> None:
    circ = Circuit(5)
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ.H(0)
    circ.H(1)
    circ.H(2)
    circ.H(3)
    circ.H(4)
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

    out_circ = route(circ, arc, initial_mapping=init_map)
    assert out_circ.valid_connectivity(arc, False)
    # oq_fidelities = [
    #     [Node(0), OpType.H, 0.9],
    #     [Node(1), OpType.H, 0.3],
    #     [Node(2), OpType.H, 0.5],
    #     [Node(3), OpType.H, 0.67],
    #     [Node(4), OpType.H, 0.99999],
    # ]

    # _commute_single_gates_through_swaps(out_circ,arc,oq_fidelities) TODO: UN COMMENT WHEN THE DEVICE CLASS IS EXPOSED!!
    # Transform.CommuteSQThroughSWAP(devi).apply(out_circ)

    Transform.DecomposeSWAPtoCX(arc).apply(out_circ)
    Transform.DecomposeCXDirected(arc).apply(out_circ)
    assert out_circ.valid_connectivity(arc, True)


def test_noncontiguous_arc() -> None:
    arc = Architecture([[0, 2]])
    pass1 = DefaultMappingPass(arc)
    c = Circuit(2)
    pass1.apply(c)


def test_noncontiguous_arc_phase_poly() -> None:
    # testing non-contiguous ascending named nodes
    arc = Architecture([[0, 2]])
    pass1 = AASRouting(arc, lookahead=1)
    c = Circuit(2).H(0).H(1)
    pass1.apply(c)
    assert c.n_gates_of_type(OpType.H) == 2
    assert c.n_gates_of_type(OpType.CX) == 0
    assert c.n_gates_of_type(OpType.CX) == 0


def test_RoutingPass() -> None:
    arc = Architecture([[0, 2], [1, 3], [2, 3], [2, 4]])
    circ = Circuit(5)
    circ.CX(0, 1)
    circ.CX(0, 3)
    circ.CX(2, 4)
    circ.CX(1, 4)
    circ.CX(1, 3)
    circ.CX(1, 2)
    cu_0 = CompilationUnit(circ)
    cu_1 = CompilationUnit(circ)
    placer = GraphPlacement(arc)
    p_pass = PlacementPass(placer)
    r_pass_0 = RoutingPass(arc, swap_lookahead=10, bridge_interactions=10)
    r_pass_1 = RoutingPass(arc, swap_lookahead=10, bridge_interactions=0)
    p_pass.apply(cu_0)
    p_pass.apply(cu_1)
    r_pass_0.apply(cu_0)
    r_pass_1.apply(cu_1)
    out_circ_0 = cu_0.circuit
    out_circ_1 = cu_1.circuit
    # TODO Should we expect BRIDGE gates in out_circ_0? If not, replace with an example
    # where we would. See See https://github.com/CQCL-DEV/tket/pull/747.
    # assert out_circ_0.n_gates_of_type(OpType.BRIDGE) == 1
    assert out_circ_0.valid_connectivity(arc, False, True)
    assert out_circ_1.n_gates_of_type(OpType.BRIDGE) == 0
    assert out_circ_1.valid_connectivity(arc, False, True)


def test_FullMappingPass() -> None:
    arc = Architecture([[0, 2], [1, 3], [2, 3], [2, 4]])
    circ = Circuit(5)
    circ.CX(0, 1).CX(0, 3).CX(2, 4).CX(1, 4).CX(0, 4).CX(2, 1).CX(3, 0)
    cu_0 = CompilationUnit(circ)
    cu_1 = CompilationUnit(circ)
    gp_placer = GraphPlacement(arc)
    lp_placer = LinePlacement(arc)
    m_pass_0 = FullMappingPass(
        arc, gp_placer, swap_lookahead=10, bridge_interactions=10
    )
    m_pass_1 = FullMappingPass(arc, lp_placer)
    m_pass_0.apply(cu_0)
    m_pass_1.apply(cu_1)
    out_circ_0 = cu_0.circuit
    out_circ_1 = cu_1.circuit
    assert out_circ_0.n_gates < out_circ_1.n_gates
    assert out_circ_0.valid_connectivity(arc, False, True)
    assert out_circ_1.valid_connectivity(arc, False, True)


def test_AAS() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_2() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc)
    assert pass1.apply(circ)


def test_AAS_3() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_4() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc)
    assert pass1.apply(circ)


def test_AAS_5() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_6() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc)
    assert pass1.apply(circ)


def test_AAS_7() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_8() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.CX(0, 1)
    circ.H(0)
    circ.Z(1)
    circ.CX(0, 3)
    circ.Rx(1.5, 3)
    circ.CX(2, 4)
    circ.X(2)
    circ.CX(1, 4)
    circ.CX(0, 4)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_9() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]])
    circ = Circuit(9)
    circ.CX(0, 8).CX(8, 1).CX(1, 7).CX(7, 2).CX(2, 6).CX(6, 3).CX(3, 5).CX(5, 4)
    circ.Rz(0.5, 4)
    pass1 = AASRouting(arc, lookahead=2)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() < 56


def test_AAS_10() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=2)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() < 33


def test_AAS_11() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.SWAP)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 119


def test_AAS_12() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.HamPath)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 36


def test_AAS_13() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.Rec)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 28


def test_AAS_14() -> None:
    arc = Architecture([[0, 1], [1, 0], [1, 2], [2, 1]])
    circ = Circuit(3).CZ(0, 1)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.Rec)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 3


def test_AAS_15() -> None:
    arc = Architecture([[0, 1], [1, 0], [1, 2], [2, 1]])
    circ = Circuit(2).CZ(0, 1)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.Rec)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 3


def test_CXMappingPass() -> None:
    arc = Architecture([[0, 2], [1, 3], [2, 3], [2, 4]])
    circ = Circuit(5)
    circ.Y(4).CX(0, 1).S(3).CX(0, 3).H(0).CX(2, 4).CX(1, 4).Y(1).CX(0, 4).CX(2, 1).Z(
        2
    ).CX(3, 0).CX(2, 0).CX(1, 3)
    circ.measure_all()
    cu_0 = CompilationUnit(circ)
    cu_1 = CompilationUnit(circ)
    gp_placer = GraphPlacement(arc)
    lp_placer = LinePlacement(arc)
    m_pass_0 = CXMappingPass(
        arc, gp_placer, swap_lookahead=10, bridge_interactions=10, directed_cx=True
    )
    m_pass_1 = CXMappingPass(arc, lp_placer, delay_measures=False)
    m_pass_0.apply(cu_0)
    m_pass_1.apply(cu_1)
    out_circ_0 = cu_0.circuit
    out_circ_1 = cu_1.circuit

    measure_pred = NoMidMeasurePredicate()
    assert measure_pred.verify(cu_0.circuit) == True
    assert measure_pred.verify(cu_1.circuit) == False
    assert out_circ_0.valid_connectivity(arc, True)
    assert out_circ_1.valid_connectivity(arc, False)


def test_CXMappingPass_correctness() -> None:
    # TKET-1045
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    placer = NoiseAwarePlacement(arc)
    p = CXMappingPass(arc, placer, directed_cx=True, delay_measures=True)
    c = Circuit(3).CX(0, 1).CX(1, 2).CCX(2, 1, 0).CY(1, 0).CY(2, 1)
    cu = CompilationUnit(c)
    p.apply(cu)
    c1 = cu.circuit
    u1 = c1.get_unitary()
    assert all(np.isclose(abs(x), 0) or np.isclose(abs(x), 1) for x in u1.flatten())


def test_place_with_map_twice() -> None:
    # TKET-671
    c = Circuit(6).CX(0, 1).CX(2, 3).CX(4, 3).CX(2, 4).CX(3, 5).CX(0, 5)

    index_map = {0: 1, 1: 2, 2: 0, 3: 4, 4: 3}
    uid_map = {Qubit(i): Node(j) for i, j in index_map.items()}
    c_qbs = c.qubits
    assert uid_map[c_qbs[0]] == Node(1)
    assert uid_map[c_qbs[1]] == Node(2)
    assert uid_map[c_qbs[2]] == Node(0)
    assert uid_map[c_qbs[3]] == Node(4)
    assert uid_map[c_qbs[4]] == Node(3)

    assert all(qb.reg_name == "q" for qb in c.qubits)
    place_with_map(c, uid_map)
    assert all(qb.reg_name in ["node", "unplaced"] for qb in c.qubits)
    place_with_map(c, uid_map)
    assert all(qb.reg_name == "unplaced" for qb in c.qubits)


def test_big_placement() -> None:
    # TKET-1275
    c = circuit_from_qasm(
        Path(__file__).resolve().parent / "qasm_test_files" / "test14.qasm"
    )
    arc = Architecture(
        [
            [0, 1],
            [0, 14],
            [1, 0],
            [1, 2],
            [1, 13],
            [2, 1],
            [2, 3],
            [2, 12],
            [3, 2],
            [3, 4],
            [3, 11],
            [4, 3],
            [4, 5],
            [4, 10],
            [5, 4],
            [5, 6],
            [5, 9],
            [6, 5],
            [6, 8],
            [7, 8],
            [8, 6],
            [8, 7],
            [8, 9],
            [9, 5],
            [9, 8],
            [9, 10],
            [10, 4],
            [10, 9],
            [10, 11],
            [11, 3],
            [11, 10],
            [11, 12],
            [12, 2],
            [12, 11],
            [12, 13],
            [13, 1],
            [13, 12],
            [13, 14],
            [14, 0],
            [14, 13],
        ]
    )
    assert PauliSimp().apply(c)
    assert DefaultMappingPass(arc).apply(c)


def test_CXMappingPass_terminates() -> None:
    # TKET-1376
    c = circuit_from_qasm(
        Path(__file__).resolve().parent / "qasm_test_files" / "test13.qasm"
    )
    arc = Architecture(
        [
            [0, 1],
            [1, 0],
            [1, 2],
            [1, 4],
            [2, 1],
            [2, 3],
            [3, 2],
            [3, 5],
            [4, 1],
            [4, 7],
            [5, 3],
            [5, 8],
            [6, 7],
            [7, 4],
            [7, 6],
            [7, 10],
            [8, 5],
            [8, 9],
            [8, 11],
            [9, 8],
            [10, 7],
            [10, 12],
            [11, 8],
            [11, 14],
            [12, 10],
            [12, 13],
            [12, 15],
            [13, 12],
            [13, 14],
            [14, 11],
            [14, 13],
            [14, 16],
            [15, 12],
            [15, 18],
            [16, 14],
            [16, 19],
            [17, 18],
            [18, 15],
            [18, 17],
            [18, 21],
            [19, 16],
            [19, 20],
            [19, 22],
            [20, 19],
            [21, 18],
            [21, 23],
            [22, 19],
            [22, 25],
            [23, 21],
            [23, 24],
            [24, 23],
            [24, 25],
            [25, 22],
            [25, 24],
            [25, 26],
            [26, 25],
        ]
    )
    placer = NoiseAwarePlacement(arc)
    placer.modify_config(timeout=10000)
    p = CXMappingPass(arc, placer, directed_cx=False, delay_measures=False)
    assert p.apply(c)


if __name__ == "__main__":
    test_architectures()
    test_placements()
    test_placement_config()
    test_convert_index_mapping()
    test_basic_routing()
    test_basic_routing_with_line_map()
    test_commuting_sq_through_swap()
    test_decompose_swap_to_cx()
    test_greedy_noise_route()
    test_basic_routing_with_noise_map()
    test_noncontiguous_arc()
    test_noncontiguous_arc_phase_poly()
    test_RoutingPass()
    test_FullMappingPass()
    test_CXMappingPass()
    test_place_with_map_twice()
