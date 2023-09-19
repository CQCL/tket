from typing import cast

import numpy as np

from pytket.circuit import Circuit, Qubit, PhasePolyBox
from pytket.passes import ComposePhasePolyBoxes, DecomposeBoxes
from pytket.utils import compare_unitaries
from useful_typedefs import ParamType  # type: ignore

PhasePolynomial = list[tuple[list[bool], ParamType]]


def phase_polynomials_are_equal(
    phase_poly_0: PhasePolynomial, phase_poly_1: PhasePolynomial
) -> bool:
    to_compare_0 = {tuple(pair[0]): pair[1] for pair in phase_poly_0}
    to_compare_1 = {tuple(pair[0]): pair[1] for pair in phase_poly_1}
    return to_compare_0 == to_compare_1


def test_phase_polybox() -> None:
    c = Circuit(1, 1)
    n_qb = 1
    qubit_indices = {Qubit(0): 0}
    phase_polynomial: PhasePolynomial = [([True], 0.1)]
    linear_transformation = np.array([[1]])
    p_box = PhasePolyBox(n_qb, qubit_indices, phase_polynomial, linear_transformation)

    b = p_box.get_circuit()

    p_box_ii = PhasePolyBox(b)

    c.add_phasepolybox(p_box, [0])
    c.add_phasepolybox(p_box, [Qubit(0)])

    c.add_phasepolybox(p_box_ii, [0])

    assert p_box.n_qubits == n_qb
    assert p_box_ii.n_qubits == n_qb
    assert p_box.qubit_indices == qubit_indices
    assert p_box_ii.qubit_indices == qubit_indices
    assert phase_polynomials_are_equal(
        cast(PhasePolynomial, p_box.phase_polynomial_as_list), phase_polynomial
    )
    assert phase_polynomials_are_equal(
        cast(PhasePolynomial, p_box_ii.phase_polynomial_as_list), phase_polynomial
    )
    assert np.array_equal(p_box.linear_transformation, linear_transformation)
    assert np.array_equal(p_box_ii.linear_transformation, linear_transformation)
    assert DecomposeBoxes().apply(c)


def test_phase_polybox_II() -> None:
    c = Circuit(1, 1)
    n_qb = 1
    qubit_indices = {Qubit(0): 0}
    phase_polynomial: PhasePolynomial = [([True], 0.1), ([True], 0.3)]
    linear_transformation = np.array([[1]])
    p_box = PhasePolyBox(n_qb, qubit_indices, phase_polynomial, linear_transformation)

    b = p_box.get_circuit()

    p_box_ii = PhasePolyBox(b)

    c.add_phasepolybox(p_box, [0])
    c.add_phasepolybox(p_box, [Qubit(0)])

    c.add_phasepolybox(p_box_ii, [0])

    assert p_box.n_qubits == n_qb
    assert p_box_ii.n_qubits == n_qb
    assert p_box.qubit_indices == qubit_indices
    assert p_box_ii.qubit_indices == qubit_indices
    assert phase_polynomials_are_equal(
        cast(PhasePolynomial, p_box.phase_polynomial_as_list), phase_polynomial
    )
    assert phase_polynomials_are_equal(
        cast(PhasePolynomial, p_box_ii.phase_polynomial_as_list), phase_polynomial
    )
    assert np.array_equal(p_box.linear_transformation, linear_transformation)
    assert np.array_equal(p_box_ii.linear_transformation, linear_transformation)
    assert DecomposeBoxes().apply(c)


def test_phase_polybox_big() -> None:
    c = Circuit(3, 3)
    n_qb = 3
    qubit_indices = {Qubit(0): 0, Qubit(1): 1, Qubit(2): 2}
    phase_polynomial: PhasePolynomial = [
        ([True, False, True], 0.333),
        ([False, False, True], 0.05),
        ([False, True, False], 1.05),
    ]
    linear_transformation = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    p_box = PhasePolyBox(n_qb, qubit_indices, phase_polynomial, linear_transformation)

    b = p_box.get_circuit()

    p_box_ii = PhasePolyBox(b)

    c.add_phasepolybox(p_box, [0, 1, 2])
    c.add_phasepolybox(p_box_ii, [0, 1, 2])
    assert p_box.n_qubits == n_qb
    assert p_box_ii.n_qubits == n_qb
    assert p_box.qubit_indices == qubit_indices
    assert p_box_ii.qubit_indices == qubit_indices
    assert phase_polynomials_are_equal(
        cast(PhasePolynomial, p_box.phase_polynomial_as_list), phase_polynomial
    )
    assert phase_polynomials_are_equal(
        cast(PhasePolynomial, p_box_ii.phase_polynomial_as_list), phase_polynomial
    )
    assert np.array_equal(p_box.linear_transformation, linear_transformation)
    assert np.array_equal(p_box_ii.linear_transformation, linear_transformation)
    assert DecomposeBoxes().apply(c)


def test_compose_phase_polybox_default_registers() -> None:
    c = Circuit(3)
    c.CX(0, 2).CX(1, 0)
    d = c.copy()
    ComposePhasePolyBoxes().apply(c)
    assert compare_unitaries(d.get_unitary(), c.get_unitary())


def test_compose_phase_polybox_one_registers() -> None:
    c = Circuit(3)
    areg = c.add_q_register("a", 3)
    c.CX(areg[0], areg[2]).CX(areg[1], areg[0])
    d = c.copy()
    ComposePhasePolyBoxes().apply(c)
    DecomposeBoxes().apply(c)
    assert compare_unitaries(d.get_unitary(), c.get_unitary())


def test_compose_phase_polybox_custom_registers() -> None:
    c = Circuit()
    areg = c.add_q_register("a", 2)
    breg = c.add_q_register("b", 1)
    c.CX(areg[0], breg[0]).CX(areg[1], areg[0])
    ComposePhasePolyBoxes().apply(c)


def test_add_and_decompose_phase_polybox_custom_registers() -> None:
    c = Circuit()
    areg = c.add_q_register("a", 2)
    breg = c.add_q_register("b", 1)
    c.CX(areg[0], breg[0]).CX(areg[1], areg[0])
    ppbox = PhasePolyBox(c)
    d = Circuit(3)
    d.add_phasepolybox(ppbox, [0, 1, 2])
    DecomposeBoxes().apply(d)
    assert compare_unitaries(d.get_unitary(), c.get_unitary())
