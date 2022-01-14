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

import pytest  # type: ignore
from pytket.circuit import Circuit, OpType, Qubit  # type: ignore
from pytket.tableau import UnitaryTableau, UnitaryTableauBox  # type: ignore
from pytket.utils.results import compare_unitaries
import numpy as np


def test_tableau_box_from_gates() -> None:
    tab = UnitaryTableau(3)
    tab.apply_gate_at_end(OpType.H, [Qubit(0)])
    tab.apply_gate_at_end(OpType.CX, [Qubit(0), Qubit(1)])
    tab.apply_gate_at_end(OpType.V, [Qubit(1)])
    tab.apply_gate_at_end(OpType.CZ, [Qubit(2), Qubit(1)])
    tab.apply_gate_at_end(OpType.Vdg, [Qubit(2)])
    circ = Circuit(3)
    circ.H(0)
    circ.CX(0, 1)
    circ.V(1)
    circ.CZ(2, 1)
    circ.Vdg(2)
    circ = circ.dagger()
    circ.add_gate(UnitaryTableauBox(tab), [0, 1, 2])
    assert compare_unitaries(circ.get_unitary(), np.eye(8, dtype=complex))


def test_tableau_box_from_matrix() -> None:
    xx = np.asarray(
        [
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        dtype=bool,
    )
    xz = np.asarray(
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 1],
        ],
        dtype=bool,
    )
    xph = np.asarray([0, 0, 1], dtype=bool)
    zx = np.asarray(
        [
            [0, 0, 0],
            [0, 1, 0],
            [0, 0, 0],
        ],
        dtype=bool,
    )
    zz = np.asarray(
        [
            [1, 0, 0],
            [1, 1, 0],
            [0, 0, 1],
        ],
        dtype=bool,
    )
    zph = np.asarray([1, 0, 1], dtype=bool)
    box = UnitaryTableauBox(xx, xz, xph, zx, zz, zph)
    circ = Circuit(3)
    circ.add_gate(box, [0, 1, 2])
    circ.CX(0, 1)
    circ.X(0)
    circ.Vdg(1)
    circ.X(2)
    circ.S(2)
    assert compare_unitaries(circ.get_unitary(), np.eye(8, dtype=complex))
