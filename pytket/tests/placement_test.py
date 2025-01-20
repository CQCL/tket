# Copyright Quantinuum
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

import json
from pathlib import Path

import pytest

from pytket import Circuit
from pytket.architecture import Architecture, FullyConnected
from pytket.circuit import Node, Qubit
from pytket.mapping import LexiLabellingMethod, LexiRouteRoutingMethod, MappingManager
from pytket.passes import DefaultMappingPass, PauliSimp
from pytket.placement import (
    GraphPlacement,
    LinePlacement,
    NoiseAwarePlacement,
    Placement,
    place_fully_connected,
    place_with_map,
)
from pytket.qasm import circuit_from_qasm


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

    assert line_map[Qubit(0)] == Node(3)
    assert line_map[Qubit(1)] == Node(1)
    assert line_map[Qubit(2)] == Node("unplaced", 0)
    assert line_map[Qubit(3)] == Node("unplaced", 1)
    assert line_map[Qubit(4)] == Node(4)
    assert line_map[Qubit(5)] == Node(5)

    assert base_placed.qubits[0] == base_map[circ_qbs[0]]
    assert base_placed.qubits[1] == base_map[circ_qbs[1]]
    assert base_placed.qubits[2] == base_map[circ_qbs[2]]

    assert graph_map[Qubit(0)] == Node(2)
    assert graph_map[Qubit(1)] == Node(1)
    assert graph_map[Qubit(2)] == Node(3)
    assert graph_map[Qubit(3)] == Node(0)
    assert graph_map[Qubit(4)] == Node(4)
    assert graph_map[Qubit(5)] == Node(5)

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
        Path(__file__).resolve().parent / "json_test_files" / "placements.json"
    ) as f:
        d = json.load(f)
        base_pl_serial = d["base_placement"]
        line_pl_serial = d["line_placement"]
        graph_pl_serial = d["graph_placement"]
        noise_pl_serial = d["noise_placement"]

    assert Placement.from_dict(base_pl_serial).to_dict() == base_pl_serial
    assert LinePlacement.from_dict(line_pl_serial).to_dict() == line_pl_serial
    assert GraphPlacement.from_dict(graph_pl_serial).to_dict() == graph_pl_serial
    assert NoiseAwarePlacement.from_dict(noise_pl_serial).to_dict() == noise_pl_serial


def test_placement_config() -> None:
    test_coupling = [(0, 1), (1, 2), (2, 3)]
    test_architecture = Architecture(test_coupling)
    test_pl = GraphPlacement(test_architecture)
    c = Circuit(4).CX(0, 1).CX(1, 2).CX(2, 3)
    pm = test_pl.get_placement_map(c)
    assert test_pl.get_placement_map(c) == pm


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
    assert new_circ_qbs[5] == Qubit(5)

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
    assert all(qb.reg_name in ["node", "q"] for qb in c.qubits)
    place_with_map(c, uid_map)
    assert all(qb.reg_name in ["node", "q"] for qb in c.qubits)


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

    fc5 = FullyConnected(5, "fcNodetest")
    place_fully_connected(c, fc5)
    qbs = c.qubits
    assert qbs[0].reg_name == "fcNodetest"
    assert qbs[1].reg_name == "fcNodetest"
    assert qbs[2].reg_name == "fcNodetest"
    assert qbs[3].reg_name == "fcNodetest"
    assert qbs[4].reg_name == "fcNodetest"


def test_big_placement() -> None:
    # TKET-1275
    c = circuit_from_qasm(
        Path(__file__).resolve().parent / "qasm_test_files" / "test14.qasm"
    )
    arc = Architecture(
        [
            (0, 1),
            (0, 14),
            (1, 0),
            (1, 2),
            (1, 13),
            (2, 1),
            (2, 3),
            (2, 12),
            (3, 2),
            (3, 4),
            (3, 11),
            (4, 3),
            (4, 5),
            (4, 10),
            (5, 4),
            (5, 6),
            (5, 9),
            (6, 5),
            (6, 8),
            (7, 8),
            (8, 6),
            (8, 7),
            (8, 9),
            (9, 5),
            (9, 8),
            (9, 10),
            (10, 4),
            (10, 9),
            (10, 11),
            (11, 3),
            (11, 10),
            (11, 12),
            (12, 2),
            (12, 11),
            (12, 13),
            (13, 1),
            (13, 12),
            (13, 14),
            (14, 0),
            (14, 13),
        ]
    )
    assert PauliSimp().apply(c)
    assert DefaultMappingPass(arc).apply(c)


def test_large_error_rate_noise_aware() -> None:
    nodes = [Node(i) for i in range(3)]
    arc = Architecture([(nodes[0], nodes[1]), (nodes[1], nodes[2])])
    nap = NoiseAwarePlacement(
        arc,
        {},
        {
            (nodes[0], nodes[1]): 0.2,
            (nodes[1], nodes[0]): 0.5,
            (nodes[1], nodes[2]): 2,
            (nodes[2], nodes[1]): 1,
        },
    )

    c = Circuit(2).CX(0, 1)

    placement_map = nap.get_placement_map(c)
    assert placement_map[Qubit(0)] == Node(0)
    assert placement_map[Qubit(1)] == Node(1)


if __name__ == "__main__":
    test_placements()
    test_placements_serialization()
    test_convert_index_mapping()
    test_place_with_map_twice()
    test_big_placement()
    test_place_fully_connected()
    test_placement_config()
    test_large_error_rate_noise_aware()
