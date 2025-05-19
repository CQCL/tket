from collections.abc import Mapping, Sequence
import enum
from typing import Annotated, overload

from numpy.typing import NDArray
import scipy.sparse

import pytket._tket.unit_id


class Pauli(enum.Enum):
    I = 0

    X = 1

    Y = 2

    Z = 3

I: Pauli = Pauli.I

X: Pauli = Pauli.X

Y: Pauli = Pauli.Y

Z: Pauli = Pauli.Z

class QubitPauliString:
    """
    A string of Pauli letters from the alphabet {I, X, Y, Z}, implemented as a sparse list, indexed by qubit.
    """

    @overload
    def __init__(self) -> None:
        """Constructs an empty QubitPauliString."""

    @overload
    def __init__(self, qubit: pytket._tket.unit_id.Qubit, pauli: Pauli) -> None:
        """Constructs a QubitPauliString with a single Pauli term."""

    @overload
    def __init__(self, qubits: Sequence[pytket._tket.unit_id.Qubit], paulis: Sequence[Pauli]) -> None:
        """
        Constructs a QubitPauliString from two matching lists of Qubits and Paulis.
        """

    @overload
    def __init__(self, map: Mapping[pytket._tket.unit_id.Qubit, Pauli]) -> None:
        """
        Construct a QubitPauliString from a dictionary mapping :py:class:`Qubit` to :py:class:`Pauli`.
        """

    def __hash__(self) -> int: ...

    def __repr__(self) -> str: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    def __lt__(self, arg: QubitPauliString, /) -> bool: ...

    def __getitem__(self, arg: pytket._tket.unit_id.Qubit, /) -> Pauli: ...

    def __setitem__(self, arg0: pytket._tket.unit_id.Qubit, arg1: Pauli, /) -> None: ...

    @property
    def map(self) -> dict[pytket._tket.unit_id.Qubit, Pauli]:
        """
        The QubitPauliString's underlying dict mapping :py:class:`Qubit` to :py:class:`Pauli`
        """

    def to_list(self) -> list:
        """
        A JSON-serializable representation of the QubitPauliString.

        :return: a list of :py:class:`Qubit`-to-:py:class:`Pauli` entries, represented as dicts.
        """

    @staticmethod
    def from_list(arg: list, /) -> QubitPauliString:
        """
        Construct a new QubitPauliString instance from a JSON serializable list representation.
        """

    def compress(self) -> None:
        """Removes I terms to compress the sparse representation."""

    def commutes_with(self, other: QubitPauliString) -> bool:
        """:return: True if the two strings commute, else False"""

    @overload
    def to_sparse_matrix(self) -> scipy.sparse.csc_matrix[complex]:
        """
        Represents the sparse string as a dense string (without padding for extra qubits) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit("a", 0)`` is more significant that ``Qubit("a", 1)`` and ``Qubit("b")`` for indexing into the matrix.

        :return: a sparse matrix corresponding to the operator
        """

    @overload
    def to_sparse_matrix(self, n_qubits: int) -> scipy.sparse.csc_matrix[complex]:
        """
        Represents the sparse string as a dense string over `n_qubits` qubits (sequentially indexed from 0 in the default register) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit(0)`` is the most significant bit for indexing into the matrix.

        :param n_qubits: the number of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """

    @overload
    def to_sparse_matrix(self, qubits: Sequence[pytket._tket.unit_id.Qubit]) -> scipy.sparse.csc_matrix[complex]:
        """
        Represents the sparse string as a dense string and generates the matrix for the tensor. Orders qubits according to `qubits` (padding with identities if they are not in the sparse string), so ``qubits[0]`` is the most significant bit for indexing into the matrix.

        :param qubits: the ordered list of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """

    @overload
    def dot_state(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]) -> Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]:
        """
        Performs the dot product of the state with the pauli string. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.

        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: dot product of operator with state
        """

    @overload
    def dot_state(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')], qubits: Sequence[pytket._tket.unit_id.Qubit]) -> Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]:
        """
        Performs the dot product of the state with the pauli string. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.

        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: dot product of operator with state
        """

    @overload
    def state_expectation(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]) -> complex:
        """
        Calculates the expectation value of the state with the pauli string. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.

        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: expectation value with respect to state
        """

    @overload
    def state_expectation(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')], qubits: Sequence[pytket._tket.unit_id.Qubit]) -> complex:
        """
        Calculates the expectation value of the state with the pauli string. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.

        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: expectation value with respect to state
        """

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...

def pauli_string_mult(qubitpaulistring1: QubitPauliString, qubitpaulistring2: QubitPauliString) -> tuple[QubitPauliString, complex]:
    """
    :return: the product of two QubitPauliString objects as a pair (QubitPauliString, complex)
    """

class PauliStabiliser:
    """
    A string of Pauli letters from the alphabet {I, X, Y, Z} with a +/- 1 coefficient.
    """

    @overload
    def __init__(self) -> None:
        """Constructs an empty PauliStabiliser."""

    @overload
    def __init__(self, string: Sequence[Pauli], coeff: int) -> None:
        """Constructs a PauliStabiliser with a list of Pauli terms."""

    @property
    def coeff(self) -> int:
        """The coefficient of the stabiliser"""

    @property
    def string(self) -> list[Pauli]:
        """The list of Pauli terms"""

    def __eq__(self, arg: object, /) -> bool: ...

    def __hash__(self) -> int:
        """
        Hashing is not implemented for this class, attempting to hash an object will raise a type error
        """

    def __ne__(self, arg: object, /) -> bool: ...

class QubitPauliTensor:
    """
    A tensor formed by Pauli terms, consisting of a sparse map from :py:class:`Qubit` to :py:class:`Pauli` (implemented as a :py:class:`QubitPauliString`) and a complex coefficient.
    """

    @overload
    def __init__(self, coeff: complex = 1.0) -> None:
        """Constructs an empty QubitPauliTensor, representing the identity."""

    @overload
    def __init__(self, qubit: pytket._tket.unit_id.Qubit, pauli: Pauli, coeff: complex = 1.0) -> None:
        """Constructs a QubitPauliTensor with a single Pauli term."""

    @overload
    def __init__(self, qubits: Sequence[pytket._tket.unit_id.Qubit], paulis: Sequence[Pauli], coeff: complex = 1.0) -> None:
        """
        Constructs a QubitPauliTensor from two matching lists of Qubits and Paulis.
        """

    @overload
    def __init__(self, map: Mapping[pytket._tket.unit_id.Qubit, Pauli], coeff: complex = 1.0) -> None:
        """
        Construct a QubitPauliTensor from a dictionary mapping :py:class:`Qubit` to :py:class:`Pauli`.
        """

    @overload
    def __init__(self, string: QubitPauliString, coeff: complex = 1.0) -> None:
        """Construct a QubitPauliTensor from a QubitPauliString."""

    def __hash__(self) -> int: ...

    def __repr__(self) -> str: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    def __lt__(self, arg: QubitPauliTensor, /) -> bool: ...

    def __getitem__(self, arg: pytket._tket.unit_id.Qubit, /) -> Pauli: ...

    def __setitem__(self, arg0: pytket._tket.unit_id.Qubit, arg1: Pauli, /) -> None: ...

    def __mul__(self, arg: QubitPauliTensor, /) -> QubitPauliTensor: ...

    def __rmul__(self, arg: complex, /) -> QubitPauliTensor: ...

    @property
    def string(self) -> QubitPauliString:
        """The QubitPauliTensor's underlying :py:class:`QubitPauliString`"""

    @string.setter
    def string(self, arg: QubitPauliString, /) -> None: ...

    @property
    def coeff(self) -> complex:
        """The global coefficient of the tensor"""

    @coeff.setter
    def coeff(self, arg: complex, /) -> None: ...

    def compress(self) -> None:
        """Removes I terms to compress the sparse representation."""

    def commutes_with(self, other: QubitPauliTensor) -> bool:
        """:return: True if the two tensors commute, else False"""

    @overload
    def to_sparse_matrix(self) -> scipy.sparse.csc_matrix[complex]:
        """
        Represents the sparse string as a dense string (without padding for extra qubits) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit("a", 0)`` is more significant that ``Qubit("a", 1)`` and ``Qubit("b")`` for indexing into the matrix.

        :return: a sparse matrix corresponding to the tensor
        """

    @overload
    def to_sparse_matrix(self, n_qubits: int) -> scipy.sparse.csc_matrix[complex]:
        """
        Represents the sparse string as a dense string over `n_qubits` qubits (sequentially indexed from 0 in the default register) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit(0)`` is the most significant bit for indexing into the matrix.

        :param n_qubits: the number of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """

    @overload
    def to_sparse_matrix(self, qubits: Sequence[pytket._tket.unit_id.Qubit]) -> scipy.sparse.csc_matrix[complex]:
        """
        Represents the sparse string as a dense string and generates the matrix for the tensor. Orders qubits according to `qubits` (padding with identities if they are not in the sparse string), so ``qubits[0]`` is the most significant bit for indexing into the matrix.

        :param qubits: the ordered list of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """

    @overload
    def dot_state(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]) -> Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]:
        """
        Performs the dot product of the state with the pauli tensor. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.

        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: dot product of operator with state
        """

    @overload
    def dot_state(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')], qubits: Sequence[pytket._tket.unit_id.Qubit]) -> Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]:
        """
        Performs the dot product of the state with the pauli tensor. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.

        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: dot product of operator with state
        """

    @overload
    def state_expectation(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')]) -> complex:
        """
        Calculates the expectation value of the state with the pauli operator. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.

        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: expectation value with respect to state
        """

    @overload
    def state_expectation(self, state: Annotated[NDArray, dict(dtype='complex128', shape=(None), order='C')], qubits: Sequence[pytket._tket.unit_id.Qubit]) -> complex:
        """
        Calculates the expectation value of the state with the pauli operator. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.

        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: expectation value with respect to state
        """

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...
