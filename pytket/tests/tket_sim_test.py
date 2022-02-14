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
from pytket.circuit import Circuit  # type: ignore
from enum import Enum  # type: ignore
import numpy as np
import math  # type: ignore
from typing import Any, Tuple

# Note: of course, one could write many more (circuit -> unitary) tests.
# But the functions are just wrappers around Simulation functions in tket,
# which already has quite extensive C++ unit tests.
# We worry about and test for conversion errors between (C++, Python)
# and (Eigen, numpy) more than specific numerical values.

# Note: NumPy, and the way pybind11 converts between Eigen and NumPy,
# is annoying, hence why there are slightly strange-looking tests with matrices
# and columnvectors in many different formats.
# Since NumPy is standard, people probably just have to put up with row/column
# vectors not being handled sensibly. I tried various methods in main.cpp
# to make numpy behave, but they didn't work. If we really don't like the
# row/column vector behaviour, I think it needs an extra layer of Python code
# somehow to force the data into NumPy matrix formats.


class Behaviour(Enum):
    """Eigen <--> numpy conversion has some oddities.
    Although Eigen is pretty sensible,
    numpy can have surprising behaviour.
    Use this enum to document the actual results."""

    SUCCESS = 1
    # This may be surprising, depending on what you expect
    RESULT_NEEDS_RESHAPING = 2
    INVALID_MATRIX_PRODUCT = 3


def append_gates_sequence_0(circ: Circuit) -> None:
    """Append a fixed sequence of gates to the circuit."""
    circ.H(0)
    circ.XXPhase(0.1, 1, 0)
    circ.Rx(0.7, 1)


def append_gates_sequence_1(circ: Circuit) -> None:
    """Append a fixed sequence of gates to the circuit."""
    circ.YYPhase(0.2, 0, 1)
    circ.Ry(0.3, 1)
    circ.Rz(0.15, 0)


def get_circuit_triple() -> Tuple[Circuit, Circuit, Circuit]:
    """Returns 3 circuits, where the final circuit is the first
    followed by the second."""
    c0 = Circuit(2)
    append_gates_sequence_0(c0)
    c1 = Circuit(2)
    append_gates_sequence_1(c1)
    c2 = Circuit(2)
    append_gates_sequence_0(c2)
    append_gates_sequence_1(c2)
    return (c0, c1, c2)


def get_checked_unitary_matrix(circ: Circuit) -> np.ndarray:
    """Asserts that the matrix of a circuit really is almost unitary,
    of the correct size, and returns it."""
    expected_size = math.pow(2, circ.n_qubits)
    assert int(expected_size) == expected_size
    expected_size = int(expected_size)
    unitary = np.asarray(circ.get_unitary())
    assert len(unitary) == expected_size
    assert len(unitary[0]) == expected_size
    almost_identity = unitary @ unitary.transpose().conjugate()
    assert np.allclose(almost_identity, np.eye(expected_size, dtype=np.cdouble))
    return unitary


def check_matmul_failure_exception_string(exception_string: str) -> None:
    """When we deliberately multiply incorrectly sized matrices,
    we expect an exception."""
    assert "matmul: Input operand " in exception_string
    assert " has a mismatch in its core dimension " in exception_string


def check_that_premultiplication_fails(
    circ: Circuit, unitary: np.ndarray, matr: Any
) -> None:
    """Matrices deliberately have the wrong size.
    Check that an exception is generated, with the correct message.
    The type of "matr" is deliberately NOT specified, because we want
    to see what happens in normal Python code where we pass in objects
    which are implicitly converted to NumPy objects (or not)."""
    with pytest.raises(ValueError) as e1:
        product = unitary @ matr
    check_matmul_failure_exception_string(str(e1.value))

    with pytest.raises(RuntimeError) as e2:
        product = circ.get_unitary_times_other(matr)
    message = str(e2)
    # This error message came from the tket_sim C++
    assert "Error trying to simulate circuit" in message
    assert "premultiplying M with" in message
    assert "M has wrong number of" in message


def test_premultiplication() -> None:
    """Check that the "get_unitary_times_other" function returns the same
    result as calculating the unitary directly and multiplying."""

    # The precise object we want to use as a matrix,
    # followed by the expected result.
    initial_matrices_data = [
        (((1, 2), (3, 4), (5, 6), (7, 8)), Behaviour.SUCCESS),
        (((1, 2, -1), (3, 0, 4), (5, -9.72342, 6), (7, 2j, 8)), Behaviour.SUCCESS),
        # When passed into C++, arrays are silently converted
        # into an Eigen column vector, NOT a row vector.
        ([1, 2, 3, 4], Behaviour.RESULT_NEEDS_RESHAPING),
        (np.asarray([1, 2, 3, 4]), Behaviour.RESULT_NEEDS_RESHAPING),
        # Include some deliberately wrong sized matrices.
        # This really IS a row vector!
        (((1, 2, 3, 4),), Behaviour.INVALID_MATRIX_PRODUCT),
        (np.asarray([[1, 2, 3, 4]]), Behaviour.INVALID_MATRIX_PRODUCT),
        ([[1, 0], [0, 0]], Behaviour.INVALID_MATRIX_PRODUCT),
        ([1, 0, 0], Behaviour.INVALID_MATRIX_PRODUCT),
        ([1, 0, 0, -9, 2], Behaviour.INVALID_MATRIX_PRODUCT),
        ([[1, 0], [0, 2]], Behaviour.INVALID_MATRIX_PRODUCT),
    ]
    for circ in get_circuit_triple():
        unitary = get_checked_unitary_matrix(circ)
        for entry in initial_matrices_data:
            matr: Any = entry[0]
            expected_behaviour = entry[1]
            if expected_behaviour == Behaviour.INVALID_MATRIX_PRODUCT:
                check_that_premultiplication_fails(circ, unitary, matr)
                continue
            product = unitary @ matr
            product_again = circ.get_unitary_times_other(matr)
            assert len(product_again.shape) == 2
            assert product_again.shape[0] == 4
            if expected_behaviour == Behaviour.RESULT_NEEDS_RESHAPING:
                # Should only occur with row/column vectors
                assert product.shape == (4,)
                assert len(product) == 4
                product = np.asarray(product)
                product = product.reshape(4, 1)
            assert np.allclose(product, product_again)


def test_circuit_unitaries_homomorphism_property() -> None:
    """Converting circuits to unitaries converts circuit concatenation
    into matrix multiplication. Transposes would mess up the order,
    so extra transposes caused by C++/Python conversion
    should cause this to fail."""
    circuits = get_circuit_triple()
    assert len(circuits) == 3
    unitaries = [get_checked_unitary_matrix(c) for c in circuits]
    assert all(u.shape == (4, 4) for u in unitaries)
    # This checks that at least one matrix is changed by transposing,
    # so that we genuinely test for transpose conversion errors
    assert any(
        not np.allclose(np.asarray(u), np.asarray(u).transpose()) for u in unitaries
    )

    assert np.allclose(
        np.asarray(unitaries[2]), np.asarray(unitaries[1]) @ np.asarray(unitaries[0])
    )
    # They don't commute!
    assert not np.allclose(
        np.asarray(unitaries[2]), np.asarray(unitaries[0]) @ np.asarray(unitaries[1])
    )


def test_ry_matrix() -> None:
    """Check the specific matrix returned.
    We are worried about potential row/column order differences
    causing back/forth transpose conversion errors between Eigen
    and numpy, so we choose something which is different when transposed."""
    alpha = 0.1
    angle = math.pi * 0.5 * alpha
    cos = math.cos(angle)
    sin = math.sin(angle)

    expected_unitary = np.asarray([[cos, -sin], [sin, cos]])
    assert not np.allclose(expected_unitary, expected_unitary.transpose())

    circ = Circuit(1)
    circ.Ry(alpha, 0)

    calc_unitary = get_checked_unitary_matrix(circ)
    assert np.allclose(expected_unitary, calc_unitary)


def test_statevector() -> None:
    # The precise object we want to use as a column vector,
    # followed by the expected result.
    column_vectors_data = (
        ([[1], [0], [0], [0]], Behaviour.SUCCESS),
        (((1,), (0,), (0,), (0,)), Behaviour.SUCCESS),
        ((1, 0, 0, 0), Behaviour.RESULT_NEEDS_RESHAPING),
        ([1, 0, 0, 0], Behaviour.RESULT_NEEDS_RESHAPING),
        (np.asarray([[1], [0], [0], [0]]), Behaviour.SUCCESS),
        (np.asarray((1, 0, 0, 0)), Behaviour.RESULT_NEEDS_RESHAPING),
        (np.asarray([1, 0, 0, 0]), Behaviour.RESULT_NEEDS_RESHAPING),
        ([[1, 0], [0, 0]], Behaviour.INVALID_MATRIX_PRODUCT),
        ([1, 0, 0], Behaviour.INVALID_MATRIX_PRODUCT),
        ([1, 0, 0, 0, 0, 0], Behaviour.INVALID_MATRIX_PRODUCT),
    )

    for circ in get_circuit_triple():
        sv = circ.get_statevector()

        # This is necessary because of how numpy works;
        # see the pybind11 docs:
        #
        # "When returning an Eigen vector to numpy, the conversion is
        # ambiguous: a row vector of length 4 could be returned as either a
        # 1D array of length 4, or as a 2D array of size 1x4".
        assert sv.shape == (4,)
        sv = sv.reshape(4, 1)
        sv_with_premultiplication = circ.get_unitary_times_other([1, 0, 0, 0])
        assert np.allclose(sv, sv_with_premultiplication)

        for entry in column_vectors_data:
            column_vector: Any = entry[0]
            expected_behaviour = entry[1]
            try:
                unitary = get_checked_unitary_matrix(circ)
                recalc_sv = unitary @ column_vector
                recalc_sv_again = np.matmul(unitary, column_vector)
                assert np.allclose(recalc_sv, recalc_sv_again)

                if expected_behaviour == Behaviour.RESULT_NEEDS_RESHAPING:
                    # This is a numpy oddity; "matmul" and the
                    # @ operator do NOT always mean "matrix multiplication"
                    # in the pure mathematical sense. If they did, then the
                    # only possible valid result would have to be a
                    # column vector, so it would already have the correct shape.
                    assert recalc_sv.shape == (4,)
                    recalc_sv = recalc_sv.reshape(4, 1)

                assert np.allclose(sv, recalc_sv)

            except ValueError as e:
                assert expected_behaviour == Behaviour.INVALID_MATRIX_PRODUCT
                check_matmul_failure_exception_string(str(e))


if __name__ == "__main__":
    test_premultiplication()
    test_circuit_unitaries_homomorphism_property()
    test_ry_matrix()
    test_statevector()
