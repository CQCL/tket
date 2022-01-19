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

from pytket import Circuit, OpType
import pytest  # type: ignore


def gen_bell_state() -> Circuit:
    circ = Circuit(2)
    circ.H(0)
    circ.CX(0, 1)
    return circ


def test_direct_add() -> None:
    a = gen_bell_state()
    b = Circuit(1)
    b.add_gate(OpType.X, [0])
    a.add_circuit(b, [0])
    assert a.n_gates == 3
    assert a.n_qubits == 2
    assert a.depth() == 3


def test_swap_add() -> None:
    a = gen_bell_state()
    c = a.dagger()
    b = gen_bell_state()
    c.add_circuit(b, [1, 0])
    assert c.depth() == 3
    assert c.n_qubits == 2
    assert c.n_gates == 4


def test_append() -> None:
    a = gen_bell_state()
    b = gen_bell_state()
    a.append(b)
    assert a.n_gates == 4
    assert a.n_qubits == 2
    assert a.depth() == 4


def test_sequencing() -> None:
    c = Circuit(5, 5)
    c.H(0).Z(2).CX(0, 1).Rz(0.1, 2).CZ(2, 4).measure_all().get_commands()
    assert len(c.get_commands()) == 10


if __name__ == "__main__":
    test_direct_add()
    test_swap_add()
    test_append()
    test_sequencing()
