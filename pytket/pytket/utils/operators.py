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

import copy
from typing import Dict, TYPE_CHECKING, Union, List, Optional, Set, Any

import numpy
import numpy as np
from sympy import Symbol, sympify, Expr, re, im
from pytket.pauli import QubitPauliString, pauli_string_mult
from pytket.circuit import Qubit
from pytket.utils.serialization import complex_to_list, list_to_complex


CoeffTypeAccepted = Union[int, float, complex, Expr]
CoeffTypeConverted = Union[Expr]

if TYPE_CHECKING:
    from scipy.sparse import csc_matrix


def _coeff_convert(coeff: Union[CoeffTypeAccepted, str]) -> CoeffTypeConverted:
    sympy_val = sympify(coeff)  # type: ignore
    if not isinstance(sympy_val, Expr):
        raise ValueError("Unsupported value for QubitPauliString coefficient")
    return sympy_val


class QubitPauliOperator:
    """
    Generic data structure for generation of circuits and expectation
    value calculation. Contains a dictionary from QubitPauliString to
    sympy Expr. Capacity for symbolic expressions allows the operator
    to be used to generate ansätze for variational algorithms.

    Represents a mathematical object :math:`\\sum_j \\alpha_j P_j`,
    where each :math:`\\alpha_j` is a complex symbolic expression and
    :math:`P_j` is a Pauli string, i.e. :math:`P_j \\in \\{ I, X, Y,
    Z\\}^{\\otimes n}`.

    A prototypical example is a molecular Hamiltonian, for which one
    may wish to calculate the expectation value :math:`\\langle \\Psi
    | H | \\Psi \\rangle` by decomposing :math:`H` into individual
    Pauli measurements. Alternatively, one may wish to evolve a state
    by the operator :math:`e^{-iHt}` for digital quantum simulation.
    In this case, the whole operator must be decomposed into native
    operations.

    In both cases, :math:`H` may be represented by a
    QubitPauliOperator.
    """

    def __init__(
        self,
        dictionary: Optional[Dict[QubitPauliString, CoeffTypeAccepted]] = None,
    ) -> None:
        self._dict: Dict[QubitPauliString, CoeffTypeConverted] = dict()
        if dictionary:
            for key, value in dictionary.items():
                self._dict[key] = _coeff_convert(value)
        self._collect_qubits()

    def __repr__(self) -> str:
        return self._dict.__repr__()

    def __getitem__(self, key: QubitPauliString) -> CoeffTypeConverted:
        return self._dict[key]

    def get(
        self, key: QubitPauliString, default: CoeffTypeAccepted
    ) -> CoeffTypeConverted:
        return self._dict.get(key, _coeff_convert(default))

    def __setitem__(self, key: QubitPauliString, value: CoeffTypeAccepted) -> None:
        """Update value in dictionary ([]). Automatically converts value into sympy
        Expr.

        :param key: String to use as key
        :type key: QubitPauliString
        :param value: Associated coefficient
        :type value: Union[int, float, complex, Expr]
        """
        self._dict[key] = _coeff_convert(value)
        self._all_qubits.update(key.map.keys())

    def __getstate__(self) -> Dict[QubitPauliString, CoeffTypeConverted]:
        return self._dict

    def __setstate__(self, _dict: Dict[QubitPauliString, CoeffTypeConverted]) -> None:
        # values assumed to be already sympified
        self._dict = _dict
        self._collect_qubits()

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QubitPauliOperator):
            return self._dict == other._dict
        return False

    def __iadd__(self, addend: "QubitPauliOperator") -> "QubitPauliOperator":
        """In-place addition (+=) of QubitPauliOperators.

        :param addend: The operator to add
        :type addend: QubitPauliOperator
        :return: Updated operator (self)
        :rtype: QubitPauliOperator
        """
        if isinstance(addend, QubitPauliOperator):
            for key, value in addend._dict.items():
                self[key] = self.get(key, 0.0) + value
            self._all_qubits.update(addend._all_qubits)
        else:
            raise TypeError("Cannot add {} to QubitPauliOperator.".format(type(addend)))

        return self

    def __add__(self, addend: "QubitPauliOperator") -> "QubitPauliOperator":
        """Addition (+) of QubitPauliOperators.

        :param addend: The operator to add
        :type addend: QubitPauliOperator
        :return: Sum operator
        :rtype: QubitPauliOperator
        """
        summand = copy.deepcopy(self)
        summand += addend
        return summand

    def __imul__(
        self, multiplier: Union[float, Expr, "QubitPauliOperator"]
    ) -> "QubitPauliOperator":
        """In-place multiplication (*=) with QubitPauliOperator or scalar.
        Multiply coefficients and terms.

        :param multiplier: The operator or scalar to multiply
        :type multiplier: Union[QubitPauliOperator, int, float, complex, Expr]
        :return: Updated operator (self)
        :rtype: QubitPauliOperator
        """

        # Handle operator of the same type
        if isinstance(multiplier, QubitPauliOperator):
            result_terms: Dict = dict()
            for left_key, left_value in self._dict.items():
                for right_key, right_value in multiplier._dict.items():
                    new_term, bonus_coeff = pauli_string_mult(left_key, right_key)
                    new_coefficient = bonus_coeff * left_value * right_value

                    # Update result dict.
                    if new_term in result_terms:
                        result_terms[new_term] += new_coefficient
                    else:
                        result_terms[new_term] = new_coefficient
            self._dict = result_terms
            self._all_qubits.update(multiplier._all_qubits)
            return self

        # Handle scalars.
        elif isinstance(multiplier, (float, Expr)):
            for key in self._dict:
                self[key] *= multiplier
            return self

        # Invalid multiplier type
        else:
            raise TypeError(
                "Cannot multiply QubitPauliOperator with {}".format(type(multiplier))
            )

    def __mul__(
        self, multiplier: Union[float, Expr, "QubitPauliOperator"]
    ) -> "QubitPauliOperator":
        """Multiplication (*) by QubitPauliOperator or scalar.

        :param multiplier: The scalar to multiply by
        :type multiplier: Union[int, float, complex, Expr, QubitPauliOperator]
        :return: Product operator
        :rtype: QubitPauliOperator
        """
        product = copy.deepcopy(self)
        product *= multiplier
        return product

    def __rmul__(self, multiplier: CoeffTypeAccepted) -> "QubitPauliOperator":
        """Multiplication (*) by a scalar.
        We only define __rmul__ for scalars because left multiply is
        queried as default behaviour, and is used for
        QubitPauliOperator*QubitPauliOperator.

        :param multiplier: The scalar to multiply by
        :type multiplier: Union[int, float, complex, Expr]
        :return: Product operator
        :rtype: QubitPauliOperator
        """
        return self * _coeff_convert(multiplier)

    @property
    def all_qubits(self) -> Set[Qubit]:
        """
        :return: The set of all qubits the operator ranges over (including qubits
            that were provided explicitly as identities)
        :rtype: Set[Qubit]
        """
        return self._all_qubits

    def subs(self, symbol_dict: Dict[Symbol, complex]) -> None:
        """Substitutes any matching symbols in the QubitPauliOperator.

        :param symbol_dict: A dictionary of symbols to fixed values.
        :type symbol_dict: Dict[Symbol, complex]
        """
        for key, value in self._dict.items():
            self._dict[key] = value.subs(symbol_dict)  # type: ignore

    def to_list(self) -> List[Dict[str, Any]]:
        """Generate a list serialized representation of QubitPauliOperator,
         suitable for writing to JSON.

        :return: JSON serializable list of dictionaries.
        :rtype: List[Dict[str, Any]]
        """
        ret: List[Dict[str, Any]] = []
        for k, v in self._dict.items():
            try:
                coeff = complex_to_list(complex(v))
            except TypeError:
                assert type(Expr(v)) == Expr  # type: ignore
                coeff = str(v)
            ret.append(
                {
                    "string": k.to_list(),
                    "coefficient": coeff,
                }
            )
        return ret

    @classmethod
    def from_list(cls, pauli_list: List[Dict[str, Any]]) -> "QubitPauliOperator":
        """Construct a QubitPauliOperator from a serializable JSON list format,
        as returned by QubitPauliOperator.to_list()

        :return: New QubitPauliOperator instance.
        :rtype: QubitPauliOperator
        """

        def get_qps(obj: Dict[str, Any]) -> QubitPauliString:
            return QubitPauliString.from_list(obj["string"])

        def get_coeff(obj: Dict[str, Any]) -> CoeffTypeConverted:
            coeff = obj["coefficient"]
            if type(coeff) is str:
                return _coeff_convert(coeff)
            else:
                return _coeff_convert(list_to_complex(coeff))

        return QubitPauliOperator({get_qps(obj): get_coeff(obj) for obj in pauli_list})

    def to_sparse_matrix(
        self, qubits: Union[List[Qubit], int, None] = None
    ) -> "csc_matrix":
        """Represents the sparse operator as a dense operator under the ordering
        scheme specified by ``qubits``, and generates the corresponding matrix.

        - When ``qubits`` is an explicit list, the qubits are ordered with
          ``qubits[0]`` as the most significant qubit for indexing into the matrix.
        - If ``None``, then no padding qubits are introduced and we use the ILO-BE
          convention, e.g. ``Qubit("a", 0)`` is more significant than
          ``Qubit("a", 1)`` or ``Qubit("b")``.
        - Giving a number specifies the number of qubits to use in the final
          operator, treated as sequentially indexed from 0 in the default register
          (padding with identities as necessary) and ordered by ILO-BE so
          ``Qubit(0)`` is the most significant.

        :param qubits: Sequencing of qubits in the matrix, either as an explicit
            list, number of qubits to pad to, or infer from the operator.
            Defaults to None
        :type qubits: Union[List[Qubit], int, None], optional
        :return: A sparse matrix representation of the operator.
        :rtype: csc_matrix
        """
        if qubits is None:
            qubits_ = sorted(list(self._all_qubits))
            return sum(
                complex(coeff) * pauli.to_sparse_matrix(qubits_)
                for pauli, coeff in self._dict.items()
            )
        return sum(
            complex(coeff) * pauli.to_sparse_matrix(qubits)
            for pauli, coeff in self._dict.items()
        )

    def dot_state(
        self, state: np.ndarray, qubits: Optional[List[Qubit]] = None
    ) -> np.ndarray:
        """Applies the operator to the given state, mapping qubits to indexes
        according to ``qubits``.

        - When ``qubits`` is an explicit list, the qubits are ordered with
          ``qubits[0]`` as the most significant qubit for indexing into ``state``.
        - If ``None``, qubits sequentially indexed from 0 in the default register
          and ordered by ILO-BE so ``Qubit(0)`` is the most significant.

        :param state: The initial statevector
        :type state: numpy.ndarray
        :param qubits: Sequencing of qubits in ``state``, if not mapped to the
            default register. Defaults to None
        :type qubits: Union[List[Qubit], None], optional
        :return: The dot product of the operator with the statevector
        :rtype: numpy.ndarray
        """
        if qubits:
            product_sum = sum(
                complex(coeff) * pauli.dot_state(state, qubits)
                for pauli, coeff in self._dict.items()
            )
        else:
            product_sum = sum(
                complex(coeff) * pauli.dot_state(state)
                for pauli, coeff in self._dict.items()
            )
        return product_sum if isinstance(product_sum, numpy.ndarray) else state

    def state_expectation(
        self, state: np.ndarray, qubits: Optional[List[Qubit]] = None
    ) -> complex:
        """Calculates the expectation value of the given statevector with respect
        to the operator, mapping qubits to indexes according to ``qubits``.

        - When ``qubits`` is an explicit list, the qubits are ordered with
          ``qubits[0]`` as the most significant qubit for indexing into ``state``.
        - If ``None``, qubits sequentially indexed from 0 in the default register
          and ordered by ILO-BE so ``Qubit(0)`` is the most significant.

        :param state: The initial statevector
        :type state: numpy.ndarray
        :param qubits: Sequencing of qubits in ``state``, if not mapped to the
            default register. Defaults to None
        :type qubits: Union[List[Qubit], None], optional
        :return: The expectation value of the statevector and operator
        :rtype: complex
        """
        if qubits:
            return sum(
                complex(coeff) * pauli.state_expectation(state, qubits)
                for pauli, coeff in self._dict.items()
            )
        return sum(
            complex(coeff) * pauli.state_expectation(state)
            for pauli, coeff in self._dict.items()
        )

    def compress(self, abs_tol: float = 1e-10) -> None:
        """Substitutes all free symbols in the QubitPauliOperator with
        1, and then removes imaginary and real components which have
        magnitudes below the tolerance. If the resulting expression is
        0, the term is removed entirely.

        Warning: This methods assumes significant expression structure
        is known a priori, and is best suited to operators which have
        simple product expressions, such as excitation operators for
        VQE ansätze and digital quantum simulation. Otherwise, it may
        remove terms relevant to computation. Each expression is of
        the form :math:`f(a_1,a_2,\\ldots,a_n)` for some symbols
        :math:`a_i`. :math:`|f(a_1,a_2,\\ldots,a_n)|` is assumed to
        monotonically increase in both real and imaginary components
        for all :math:`a_i \\in [0, 1]`.

        :param abs_tol: The threshold below which to remove values.
        :type abs_tol: float
        """

        to_delete = []
        for key, value in self._dict.items():
            placeholder = value.subs({s: 1 for s in value.free_symbols})  # type: ignore
            if abs(re(placeholder)) <= abs_tol:
                if abs(im(placeholder)) <= abs_tol:
                    to_delete.append(key)
                else:
                    self._dict[key] = im(value) * 1j
            elif abs(im(placeholder)) <= abs_tol:
                self._dict[key] = re(value)

        for key in to_delete:
            del self._dict[key]

    def _collect_qubits(self) -> None:
        self._all_qubits: set[Qubit] = set()
        for key in self._dict.keys():
            for q in key.map.keys():
                self._all_qubits.add(q)
