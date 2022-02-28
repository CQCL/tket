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

from pytket.circuit import Node, Op, OpType, Circuit, Qubit, PhasePolyBox  # type: ignore
from pytket.architecture import Architecture, SquareGrid, FullyConnected  # type: ignore
import numpy as np


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
    fc = FullyConnected(2)
    assert isinstance(fc, FullyConnected)
    sg = SquareGrid(2, 2, 2)
    assert isinstance(sg, SquareGrid)


def test_valid_operation() -> None:
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 5), (5, 6)]
    arc = Architecture(edges)

    assert not arc.valid_operation([Node(1), Node(3)])
    assert arc.valid_operation([Node(0)])
    assert arc.valid_operation([Node(0), Node(1)])
    assert not arc.valid_operation([Node(0), Node(1), Node(2)])
    assert not arc.valid_operation([Node(10)])
    assert not arc.valid_operation([Node(10), Node(11), Node(15)])
    assert not arc.valid_operation([Node(0), Node(1), Node(2), Node(3)])
    assert not arc.valid_operation([Node(0), Node(4)])
    assert not arc.valid_operation([Node(0), Node(1), Node(2)])
    assert not arc.valid_operation([Node(0), Node(1), Node(4)])


if __name__ == "__main__":
    test_architectures()
    test_architecture_eq()
    test_fully_connected()
    test_arch_types()
    test_valid_operation()
