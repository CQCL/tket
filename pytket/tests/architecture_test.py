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

from jsonschema import Draft7Validator  # type: ignore
from referencing import Registry
from referencing.jsonschema import DRAFT7

from pytket.architecture import Architecture, FullyConnected, RingArch, SquareGrid
from pytket.circuit import Node

curr_file_path = Path(__file__).resolve().parent
schema_dir = curr_file_path.parent.parent / "schemas"
with open(schema_dir / "circuit_v1.json") as f:
    circ_schema = json.load(f)
with open(schema_dir / "architecture_v1.json") as f:
    arch_schema = json.load(f)
with open(schema_dir / "fullyconnected_v1.json") as f:
    fc_schema = json.load(f)

schema_store = [
    (circ_schema["$id"], DRAFT7.create_resource(circ_schema)),
    (arch_schema["$id"], DRAFT7.create_resource(arch_schema)),
    (fc_schema["$id"], DRAFT7.create_resource(fc_schema)),
]
registry: Registry = Registry().with_resources(schema_store)
arch_validator = Draft7Validator(arch_schema, registry=registry)
fc_validator = Draft7Validator(fc_schema, registry=registry)


def check_arch_serialisation(arch: Architecture) -> None:
    serialised_arch = arch.to_dict()
    arch_validator.validate(serialised_arch)
    new_arch = Architecture.from_dict(serialised_arch)
    new_serialised_arch = new_arch.to_dict()
    assert new_serialised_arch == serialised_arch


def check_fc_serialisation(fc: FullyConnected) -> None:
    serialised_fc = fc.to_dict()
    fc_validator.validate(serialised_fc)
    new_fc = FullyConnected.from_dict(serialised_fc)
    new_serialised_fc = new_fc.to_dict()
    assert new_serialised_fc == serialised_fc


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
    g00, g01, g10, g11 = (
        Node("gridNode", [i, j, 0]) for i in range(2) for j in range(2)
    )
    sq_arc = Architecture([(g00, g01), (g01, g11), (g00, g10), (g10, g11)])
    assert sq_arc == SquareGrid(2, 2)
    assert sq_arc != Architecture([(g00, g01), (g01, g11), (g00, g10)])


def test_fully_connected() -> None:
    fc = FullyConnected(3)
    assert fc.nodes == [Node("fcNode", i) for i in range(3)]
    d = fc.to_dict()
    assert list(d.keys()) == ["nodes"]
    fc1 = FullyConnected.from_dict(d)
    assert fc == fc1


def test_arch_types() -> None:
    arch = Architecture([(0, 1)])
    assert isinstance(arch, Architecture)
    check_arch_serialisation(arch)
    fc = FullyConnected(2)
    assert fc.nodes[0].reg_name == "fcNode"
    assert isinstance(fc, FullyConnected)
    check_fc_serialisation(fc)
    sg = SquareGrid(2, 2, 2)
    assert sg.nodes[0].reg_name == "gridNode"
    assert isinstance(sg, SquareGrid)
    check_arch_serialisation(sg)
    ra = RingArch(2)
    check_arch_serialisation(ra)
    assert ra.nodes[0].reg_name == "ringNode"


def test_arch_names() -> None:
    fc = FullyConnected(2, "fc_test")
    assert fc.nodes[0].reg_name == "fc_test"
    sg = SquareGrid(2, 2, "sg_test")
    assert sg.nodes[0].reg_name == "sg_test"
    sg = SquareGrid(2, 2, 1, "sg_test")
    assert sg.nodes[0].reg_name == "sg_test"
    ra = RingArch(2, "ring_test")
    assert ra.nodes[0].reg_name == "ring_test"


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
    test_arch_names()
