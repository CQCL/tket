# Copyright 2019-2024 Cambridge Quantum Computing
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
from copy import copy, deepcopy

from pytket.unit_id import Bit, Node, Qubit


def test_copying_qubits() -> None:
    q = Qubit(0)

    q1 = copy(q)
    q2 = deepcopy(q)

    assert type(q) is Qubit
    assert type(q1) is Qubit
    assert type(q2) is Qubit


def test_copying_bits() -> None:
    b = Bit(0)

    b1 = copy(b)
    b2 = deepcopy(b)

    assert type(b) is Bit
    assert type(b1) is Bit
    assert type(b2) is Bit


def test_copying_nodes() -> None:
    n = Node(0)

    n1 = copy(n)
    n2 = deepcopy(n)

    assert type(n) is Node
    assert type(n1) is Node
    assert type(n2) is Node
