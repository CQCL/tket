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

from pytket.circuit import Circuit, OpType, Qubit, UnitID


def gen_bell_state(reset_start: bool = False) -> Circuit:
    circ = Circuit(2)
    if reset_start:
        circ = circ.Reset(0).Reset(1)
    circ.H(0)
    circ.CX(0, 1)
    return circ


def test_direct_add() -> None:
    a = gen_bell_state(True)
    b = Circuit(1)
    b.add_gate(OpType.X, [0])
    a.add_circuit(b, [0])
    assert a.n_gates == 5
    assert a.n_qubits == 2
    assert a.depth() == 4


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


def test_add_with_map() -> None:
    c = Circuit(3).X(0).Y(1).Z(2)
    c1 = Circuit(2, 2).CX(0, 1).measure_all()
    unit_map: dict[UnitID, UnitID] = {Qubit(0): Qubit(1), Qubit(1): Qubit(0)}
    c.add_circuit_with_map(c1, unit_map)
    assert c == Circuit(3, 2).X(0).Y(1).Z(2).CX(1, 0).Measure(0, 1).Measure(1, 0)
