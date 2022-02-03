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

import pytest  # type: ignore
from pytket import Circuit
from pytket.quipper import circuit_from_quipper
from pytket.transform import Transform  # type: ignore
from pytket.utils.results import compare_unitaries
import numpy as np
import os
from typing import Any

curr_file_path = Path(__file__).resolve().parent


def approx_equal_up_to_phase(s0: np.ndarray, s1: np.ndarray, eps: float) -> bool:
    d = np.vdot(s0, s1)  # type: ignore
    return bool(abs(abs(d) - 1) < eps)


def unitary_from_simulate_output(simout: Any, n: int) -> np.ndarray:
    """
    Convert text output by Quipper's `Simulate' program to a unitary matrix.
    `n' is the number of qubits. `simout' is the output of `Simulate'.
    The encoding is little-endian.
    """
    N = pow(2, n)
    fmt = "{0:0%db}" % n
    reps = ["|" + fmt.format(i) + ">" for i in range(N)]
    lines = simout.split("\n")
    n_lines = len(lines)
    pos = 0
    a = np.zeros((N, N), dtype=complex)
    for i in range(N):
        assert lines[pos].startswith(reps[i] + " ->")
        pos += 1
        for j in range(N):
            line = lines[pos]
            if line.startswith(" "):
                line0 = lines[pos].strip()
                if line0.startswith(reps[j] + " "):
                    valstr = (
                        line0[n + 3 :]
                        .replace(" ", "")
                        .replace("*", "")
                        .replace("i", "j")
                    )
                    a[j][i] = complex(valstr)
                    pos += 1
    return np.array(a)


def test_quipper_1() -> None:
    circ = circuit_from_quipper(
        str(curr_file_path / "quipper_test_files" / "test1.quip")
    )
    assert circ.n_qubits == 16


def test_quipper_2() -> None:
    circ = circuit_from_quipper(
        str(curr_file_path / "quipper_test_files" / "test2.quip")
    )
    assert circ.n_qubits == 4


def test_quipper_3() -> None:
    # This circuit should be the identity (up to phase).
    circ = circuit_from_quipper(
        str(curr_file_path / "quipper_test_files" / "test3.quip")
    )
    assert circ.n_qubits == 4
    Transform.OptimisePostRouting().apply(circ)
    circ0 = Circuit(4, 0)
    s = circ.get_statevector()
    s0 = circ0.get_statevector()
    assert approx_equal_up_to_phase(s, s0, 1e-10)


def test_quipper_4() -> None:
    # Check some test circuits against simulator output.
    for i in range(11):
        fname = "test4-%d.quip" % i
        fpath = str(curr_file_path / "quipper_test_files" / fname)
        circ = circuit_from_quipper(fpath)
        n_qubits = circ.n_qubits
        Transform.DecomposeBoxes().apply(circ)
        Transform.RebaseToTket().apply(circ)
        a = circ.get_unitary()
        with open(fpath + ".simout") as f:
            simout = f.read()
        b = unitary_from_simulate_output(simout, n_qubits)
        # The unitary matrices a and b should be equal up to a phase.
        assert compare_unitaries(a, b)


def test_quipper_5() -> None:
    # Invalid operation (CCH gate).
    with pytest.raises(NotImplementedError):
        circ = circuit_from_quipper(
            str(curr_file_path / "quipper_test_files" / "test5.quip")
        )
