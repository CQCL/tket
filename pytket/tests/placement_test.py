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

from pathlib import Path
from pytket import Circuit  # type: ignore
from pytket.circuit import Node, Qubit  # type: ignore
from pytket.architecture import Architecture, FullyConnected  # type: ignore
from pytket.placement import (  # type: ignore
    Placement,
    LinePlacement,
    GraphPlacement,
    NoiseAwarePlacement,
    place_with_map,
    place_fully_connected,
)
from pytket.passes import PauliSimp, DefaultMappingPass  # type: ignore
from pytket.mapping import MappingManager, LexiRouteRoutingMethod, LexiLabellingMethod  # type: ignore
from pytket.qasm import circuit_from_qasm  # type: ignore

import json


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

    mm = MappingManager(test_architecture)
    mm.route_circuit(base_placed, [LexiLabellingMethod(), LexiRouteRoutingMethod()])
    mm.route_circuit(line_placed, [LexiLabellingMethod(), LexiRouteRoutingMethod()])
    mm.route_circuit(graph_placed, [LexiLabellingMethod(), LexiRouteRoutingMethod()])

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

    mm = MappingManager(test_architecture)
    mm.route_circuit(circ1, [LexiLabellingMethod(), LexiRouteRoutingMethod()])
    mm.route_circuit(circ2, [LexiLabellingMethod(), LexiRouteRoutingMethod()])
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


def test_place_fully_connected() -> None:
    c = Circuit(5)
    fc5 = FullyConnected(5)
    place_fully_connected(c, fc5)
    qbs = c.qubits
    assert qbs[0].reg_name == "fcNode"
    assert qbs[1].reg_name == "fcNode"
    assert qbs[2].reg_name == "fcNode"
    assert qbs[3].reg_name == "fcNode"
    assert qbs[4].reg_name == "fcNode"


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


if __name__ == "__main__":
    test_placements()
    test_placements_serialization()
    test_placement_config()
    test_convert_index_mapping()
    test_place_with_map_twice()
    test_big_placement()
    test_place_fully_connected()
