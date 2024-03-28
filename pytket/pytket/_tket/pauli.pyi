from numpy.typing import NDArray
from __future__ import annotations
import numpy
import pytket._tket.unit_id
import scipy.sparse
import typing
__all__ = ['I', 'Pauli', 'PauliStabiliser', 'QubitPauliString', 'QubitPauliTensor', 'X', 'Y', 'Z', 'pauli_string_mult']
class Pauli:
    """
    Members:
    
      I
    
      X
    
      Y
    
      Z
    """
    I: typing.ClassVar[Pauli]  # value = <Pauli.I: 0>
    X: typing.ClassVar[Pauli]  # value = <Pauli.X: 1>
    Y: typing.ClassVar[Pauli]  # value = <Pauli.Y: 2>
    Z: typing.ClassVar[Pauli]  # value = <Pauli.Z: 3>
    __members__: typing.ClassVar[dict[str, Pauli]]  # value = {'I': <Pauli.I: 0>, 'X': <Pauli.X: 1>, 'Y': <Pauli.Y: 2>, 'Z': <Pauli.Z: 3>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class PauliStabiliser:
    """
    A string of Pauli letters from the alphabet {I, X, Y, Z} with a +/- 1 coefficient.
    """
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __hash__(self) -> int:
        """
        Hashing is not implemented for this class, attempting to hash an object will raise a type error
        """
    @typing.overload
    def __init__(self) -> None:
        """
        Constructs an empty PauliStabiliser.
        """
    @typing.overload
    def __init__(self, string: typing.Sequence[Pauli], coeff: int) -> None:
        """
        Constructs a PauliStabiliser with a list of Pauli terms.
        """
    def __ne__(self, arg0: typing.Any) -> bool:
        ...
    @property
    def coeff(self) -> int:
        """
        The coefficient of the stabiliser
        """
    @property
    def string(self) -> list[Pauli]:
        """
        The list of Pauli terms
        """
class QubitPauliString:
    """
    A string of Pauli letters from the alphabet {I, X, Y, Z}, implemented as a sparse list, indexed by qubit.
    """
    @staticmethod
    def from_list(arg0: list) -> QubitPauliString:
        """
        Construct a new QubitPauliString instance from a JSON serializable list representation.
        """
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __getitem__(self, arg0: pytket._tket.unit_id.Qubit) -> Pauli:
        ...
    def __getstate__(self) -> tuple:
        ...
    def __hash__(self) -> int:
        ...
    @typing.overload
    def __init__(self) -> None:
        """
        Constructs an empty QubitPauliString.
        """
    @typing.overload
    def __init__(self, qubit: pytket._tket.unit_id.Qubit, pauli: Pauli) -> None:
        """
        Constructs a QubitPauliString with a single Pauli term.
        """
    @typing.overload
    def __init__(self, qubits: typing.Sequence[pytket._tket.unit_id.Qubit], paulis: typing.Sequence[Pauli]) -> None:
        """
        Constructs a QubitPauliString from two matching lists of Qubits and Paulis.
        """
    @typing.overload
    def __init__(self, map: dict[pytket._tket.unit_id.Qubit, Pauli]) -> None:
        """
        Construct a QubitPauliString from a dictionary mapping :py:class:`Qubit` to :py:class:`Pauli`.
        """
    def __lt__(self, arg0: QubitPauliString) -> bool:
        ...
    def __ne__(self, arg0: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setitem__(self, arg0: pytket._tket.unit_id.Qubit, arg1: Pauli) -> None:
        ...
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def commutes_with(self, other: QubitPauliString) -> bool:
        """
        :return: True if the two strings commute, else False
        """
    def compress(self) -> None:
        """
        Removes I terms to compress the sparse representation.
        """
    @typing.overload
    def dot_state(self, state: NDArray[numpy.complex128]) -> NDArray[numpy.complex128]:
        """
        Performs the dot product of the state with the pauli string. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.
        
        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: dot product of operator with state
        """
    @typing.overload
    def dot_state(self, state: NDArray[numpy.complex128], qubits: typing.Sequence[pytket._tket.unit_id.Qubit]) -> NDArray[numpy.complex128]:
        """
        Performs the dot product of the state with the pauli string. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.
        
        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: dot product of operator with state
        """
    @typing.overload
    def state_expectation(self, state: NDArray[numpy.complex128]) -> complex:
        """
        Calculates the expectation value of the state with the pauli string. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.
        
        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: expectation value with respect to state
        """
    @typing.overload
    def state_expectation(self, state: NDArray[numpy.complex128], qubits: typing.Sequence[pytket._tket.unit_id.Qubit]) -> complex:
        """
        Calculates the expectation value of the state with the pauli string. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.
        
        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: expectation value with respect to state
        """
    def to_list(self) -> list:
        """
        A JSON-serializable representation of the QubitPauliString.
        
        :return: a list of :py:class:`Qubit`-to-:py:class:`Pauli` entries, represented as dicts.
        """
    @typing.overload
    def to_sparse_matrix(self) -> scipy.sparse.csc_matrix:
        """
        Represents the sparse string as a dense string (without padding for extra qubits) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit("a", 0)`` is more significant that ``Qubit("a", 1)`` and ``Qubit("b")`` for indexing into the matrix.
        
        :return: a sparse matrix corresponding to the operator
        """
    @typing.overload
    def to_sparse_matrix(self, n_qubits: int) -> scipy.sparse.csc_matrix:
        """
        Represents the sparse string as a dense string over `n_qubits` qubits (sequentially indexed from 0 in the default register) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit(0)`` is the most significant bit for indexing into the matrix.
        
        :param n_qubits: the number of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """
    @typing.overload
    def to_sparse_matrix(self, qubits: typing.Sequence[pytket._tket.unit_id.Qubit]) -> scipy.sparse.csc_matrix:
        """
        Represents the sparse string as a dense string and generates the matrix for the tensor. Orders qubits according to `qubits` (padding with identities if they are not in the sparse string), so ``qubits[0]`` is the most significant bit for indexing into the matrix.
        
        :param qubits: the ordered list of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """
    @property
    def map(self) -> dict[pytket._tket.unit_id.Qubit, Pauli]:
        """
        :return: the QubitPauliString's underlying dict mapping :py:class:`Qubit` to :py:class:`Pauli`
        """
class QubitPauliTensor:
    """
    A tensor formed by Pauli terms, consisting of a sparse map from :py:class:`Qubit` to :py:class:`Pauli` (implemented as a :py:class:`QubitPauliString`) and a complex coefficient.
    """
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __getitem__(self, arg0: pytket._tket.unit_id.Qubit) -> Pauli:
        ...
    def __getstate__(self) -> tuple:
        ...
    def __hash__(self) -> int:
        ...
    @typing.overload
    def __init__(self, coeff: complex = 1.0) -> None:
        """
        Constructs an empty QubitPauliTensor, representing the identity.
        """
    @typing.overload
    def __init__(self, qubit: pytket._tket.unit_id.Qubit, pauli: Pauli, coeff: complex = 1.0) -> None:
        """
        Constructs a QubitPauliTensor with a single Pauli term.
        """
    @typing.overload
    def __init__(self, qubits: typing.Sequence[pytket._tket.unit_id.Qubit], paulis: typing.Sequence[Pauli], coeff: complex = 1.0) -> None:
        """
        Constructs a QubitPauliTensor from two matching lists of Qubits and Paulis.
        """
    @typing.overload
    def __init__(self, map: dict[pytket._tket.unit_id.Qubit, Pauli], coeff: complex = 1.0) -> None:
        """
        Construct a QubitPauliTensor from a dictionary mapping :py:class:`Qubit` to :py:class:`Pauli`.
        """
    @typing.overload
    def __init__(self, string: QubitPauliString, coeff: complex = 1.0) -> None:
        """
        Construct a QubitPauliTensor from a QubitPauliString.
        """
    def __lt__(self, arg0: QubitPauliTensor) -> bool:
        ...
    def __mul__(self, arg0: QubitPauliTensor) -> QubitPauliTensor:
        ...
    def __ne__(self, arg0: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __rmul__(self, arg0: complex) -> QubitPauliTensor:
        ...
    def __setitem__(self, arg0: pytket._tket.unit_id.Qubit, arg1: Pauli) -> None:
        ...
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def commutes_with(self, other: QubitPauliTensor) -> bool:
        """
        :return: True if the two tensors commute, else False
        """
    def compress(self) -> None:
        """
        Removes I terms to compress the sparse representation.
        """
    @typing.overload
    def dot_state(self, state: NDArray[numpy.complex128]) -> NDArray[numpy.complex128]:
        """
        Performs the dot product of the state with the pauli tensor. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.
        
        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: dot product of operator with state
        """
    @typing.overload
    def dot_state(self, state: NDArray[numpy.complex128], qubits: typing.Sequence[pytket._tket.unit_id.Qubit]) -> NDArray[numpy.complex128]:
        """
        Performs the dot product of the state with the pauli tensor. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.
        
        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: dot product of operator with state
        """
    @typing.overload
    def state_expectation(self, state: NDArray[numpy.complex128]) -> complex:
        """
        Calculates the expectation value of the state with the pauli operator. Maps the qubits of the statevector with sequentially-indexed qubits in the default register, with ``Qubit(0)`` being the most significant qubit.
        
        :param state: statevector for qubits ``Qubit(0)`` to ``Qubit(n-1)``
        :return: expectation value with respect to state
        """
    @typing.overload
    def state_expectation(self, state: NDArray[numpy.complex128], qubits: typing.Sequence[pytket._tket.unit_id.Qubit]) -> complex:
        """
        Calculates the expectation value of the state with the pauli operator. Maps the qubits of the statevector according to the ordered list `qubits`, with ``qubits[0]`` being the most significant qubit.
        
        :param state: statevector
        :param qubits: order of qubits in `state` from most to least significant
        :return: expectation value with respect to state
        """
    @typing.overload
    def to_sparse_matrix(self) -> scipy.sparse.csc_matrix:
        """
        Represents the sparse string as a dense string (without padding for extra qubits) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit("a", 0)`` is more significant that ``Qubit("a", 1)`` and ``Qubit("b")`` for indexing into the matrix.
        
        :return: a sparse matrix corresponding to the tensor
        """
    @typing.overload
    def to_sparse_matrix(self, n_qubits: int) -> scipy.sparse.csc_matrix:
        """
        Represents the sparse string as a dense string over `n_qubits` qubits (sequentially indexed from 0 in the default register) and generates the matrix for the tensor. Uses the ILO-BE convention, so ``Qubit(0)`` is the most significant bit for indexing into the matrix.
        
        :param n_qubits: the number of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """
    @typing.overload
    def to_sparse_matrix(self, qubits: typing.Sequence[pytket._tket.unit_id.Qubit]) -> scipy.sparse.csc_matrix:
        """
        Represents the sparse string as a dense string and generates the matrix for the tensor. Orders qubits according to `qubits` (padding with identities if they are not in the sparse string), so ``qubits[0]`` is the most significant bit for indexing into the matrix.
        
        :param qubits: the ordered list of qubits in the full operator
        :return: a sparse matrix corresponding to the operator
        """
    @property
    def coeff(self) -> complex:
        """
        The global coefficient of the tensor
        """
    @coeff.setter
    def coeff(self, arg0: complex) -> None:
        ...
    @property
    def string(self) -> QubitPauliString:
        """
        The QubitPauliTensor's underlying :py:class:`QubitPauliString`
        """
    @string.setter
    def string(self, arg1: QubitPauliString) -> None:
        ...
def pauli_string_mult(qubitpaulistring1: QubitPauliString, qubitpaulistring2: QubitPauliString) -> tuple[QubitPauliString, complex]:
    """
    :return: the product of two QubitPauliString objects as a pair (QubitPauliString, complex)
    """
I: Pauli  # value = <Pauli.I: 0>
X: Pauli  # value = <Pauli.X: 1>
Y: Pauli  # value = <Pauli.Y: 2>
Z: Pauli  # value = <Pauli.Z: 3>
