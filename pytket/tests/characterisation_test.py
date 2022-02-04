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

from pytket.circuit import Circuit, OpType, Qubit  # type: ignore
from pytket.tailoring import (  # type: ignore
    FrameRandomisation,
    PauliFrameRandomisation,
    UniversalFrameRandomisation,
    apply_clifford_basis_change,
)
from pytket.pauli import Pauli, QubitPauliString  # type: ignore

import pytest  # type: ignore


def test_single_cycle_single_frame_randomisation() -> None:
    test_fr = FrameRandomisation(
        {OpType.CX, OpType.X},
        {OpType.H},
        {
            OpType.X: {(OpType.H,): (OpType.H,)},
            OpType.CX: {(OpType.H, OpType.H): (OpType.H, OpType.H)},
        },
    )
    circ = Circuit(2).CX(0, 1).S(0).S(1).CX(0, 1)
    fr_circs = test_fr.get_all_circuits(circ)
    assert len(fr_circs) == 1


def test_multi_cycle_multi_frame_randomisation() -> None:
    test_fr = FrameRandomisation(
        {OpType.CX, OpType.X},
        {OpType.H, OpType.Y},
        {
            OpType.X: {(OpType.H,): (OpType.Y,), (OpType.Y,): (OpType.H,)},
            OpType.CX: {
                (OpType.H, OpType.H): (OpType.Y, OpType.Y),
                (OpType.H, OpType.Y): (OpType.Y, OpType.H),
                (OpType.Y, OpType.H): (OpType.H, OpType.Y),
                (OpType.Y, OpType.Y): (OpType.H, OpType.H),
            },
        },
    )
    circ = Circuit(2).CX(0, 1).X(0).S(0).S(1).CX(0, 1)
    fr_circs = test_fr.get_all_circuits(circ)
    assert len(fr_circs) == 16
    coms_0 = fr_circs[0].get_commands()
    assert coms_0[0].op.type == OpType.Y
    assert coms_0[1].op.type == OpType.Y
    assert coms_0[6].op.type == OpType.Y
    assert coms_0[7].op.type == OpType.H
    assert coms_0[10].op.type == OpType.Y
    assert coms_0[11].op.type == OpType.Y
    assert coms_0[15].op.type == OpType.H
    assert coms_0[16].op.type == OpType.H


def test_pauli_frame_randomisation() -> None:
    test_pfr = PauliFrameRandomisation()
    circ = Circuit(2).CX(0, 1)
    pfr_circs = test_pfr.get_all_circuits(circ)
    assert len(pfr_circs) == 16
    coms_0 = pfr_circs[0].get_commands()
    coms_15 = pfr_circs[15].get_commands()
    assert coms_0[0].op.type == OpType.Z
    assert coms_0[1].op.type == OpType.Z
    assert coms_0[5].op.type == OpType.noop
    assert coms_0[6].op.type == OpType.Z
    assert coms_15[0].op.type == OpType.noop
    assert coms_15[1].op.type == OpType.noop
    assert coms_15[5].op.type == OpType.noop
    assert coms_15[6].op.type == OpType.noop


def test_universal_frame_randomisation() -> None:
    test_pfr = UniversalFrameRandomisation()
    circ = Circuit(2).Rz(0.2, 0).CX(0, 1)
    pfr_circs = test_pfr.get_all_circuits(circ)
    assert len(pfr_circs) == 16
    coms_0 = pfr_circs[0].get_commands()
    coms_5 = pfr_circs[5].get_commands()
    coms_15 = pfr_circs[15].get_commands()

    assert coms_0[0].op.type == OpType.Z
    assert coms_0[1].op.type == OpType.Z
    assert coms_0[3].op.params[0] == 0.2
    assert coms_0[6].op.type == OpType.noop
    assert coms_0[7].op.type == OpType.Z

    assert coms_5[0].op.type == OpType.X
    assert coms_5[1].op.type == OpType.X
    assert coms_5[3].op.params[0] == 3.8
    assert coms_5[6].op.type == OpType.X
    assert coms_5[7].op.type == OpType.noop

    assert coms_15[0].op.type == OpType.noop
    assert coms_15[1].op.type == OpType.noop
    assert coms_15[3].op.params[0] == 0.2
    assert coms_15[6].op.type == OpType.noop
    assert coms_15[7].op.type == OpType.noop


def test_apply_clifford_basis_change() -> None:
    circ_0 = Circuit(1).H(0)
    z_op = QubitPauliString(Qubit(0), Pauli.Z)
    x_op = QubitPauliString(Qubit(0), Pauli.X)
    assert apply_clifford_basis_change(z_op, circ_0) == x_op

    circ_1 = Circuit(2).CX(0, 1)
    zz_op = QubitPauliString([Qubit(0), Qubit(1)], [Pauli.Z, Pauli.Z])
    iz_op = QubitPauliString([Qubit(0), Qubit(1)], [Pauli.I, Pauli.Z])
    assert apply_clifford_basis_change(zz_op, circ_1) == iz_op

    circ_2 = Circuit(2).H(0).CX(0, 1).S(1).X(0)
    zz_op = QubitPauliString([Qubit(0), Qubit(1)], [Pauli.Z, Pauli.Z])
    assert apply_clifford_basis_change(zz_op, circ_2) == iz_op

    circ_3 = (
        Circuit(4).H(0).H(1).H(2).H(3).CX(0, 1).CX(1, 2).CX(2, 3).S(0).Y(1).Z(2).X(3)
    )
    yxzi_op = QubitPauliString(
        [Qubit(0), Qubit(1), Qubit(2), Qubit(3)], [Pauli.Y, Pauli.X, Pauli.Z, Pauli.I]
    )
    yxyi_op = QubitPauliString(
        [Qubit(0), Qubit(1), Qubit(2), Qubit(3)], [Pauli.Y, Pauli.X, Pauli.Y, Pauli.I]
    )
    assert apply_clifford_basis_change(yxzi_op, circ_3) == yxyi_op


if __name__ == "__main__":
    test_single_cycle_single_frame_randomisation()
    test_multi_cycle_multi_frame_randomisation()
    test_pauli_frame_randomisation()
    test_universal_frame_randomisation()
    test_apply_clifford_basis_change()
