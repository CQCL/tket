from numpy.typing import NDArray
from __future__ import annotations
import numpy
import pytket._tket.circuit
import pytket._tket.pauli
import pytket._tket.unit_id
import typing
__all__ = ['UnitaryTableau', 'UnitaryTableauBox']
class UnitaryTableau:
    """
    Stabilizer tableau for a unitary in the style of Aaronson&Gottesman "Improved Simulation of Stabilizer Circuits": rows indicate the action at the output corresponding to either an X or a Z on a single input.
    """
    @typing.overload
    def __init__(self, nqb: int) -> None:
        """
        Constructs a :py:class:`UnitaryTableau` representing the identity operation over some number of qubits. Qubits will be indexed sequentially in the default register.
        
        :param nqb: The number of qubits in the unitary.
        """
    @typing.overload
    def __init__(self, xx: NDArray[numpy.bool_], xz: NDArray[numpy.bool_], xph: NDArray[numpy.bool_], zx: NDArray[numpy.bool_], zz: NDArray[numpy.bool_], zph: NDArray[numpy.bool_]) -> None:
        """
        Constructs a :py:class:`UnitaryTableau` from the binary tables of its components.
        
        :param xx: The X component of the X rows.
        :param xz: The Z component of the X rows.
        :param xph: The phases of the X rows.
        :param zx: The X component of the Z rows.
        :param zz: The Z component of the Z rows.
        :param zph: The phases of the Z rows.
        """
    @typing.overload
    def __init__(self, arg0: pytket._tket.circuit.Circuit) -> None:
        """
        Constructs a :py:class:`UnitaryTableau` from a unitary :py:class:`Circuit`. Throws an exception if the input contains non-unitary operations.
        
        :param circ: The unitary circuit to convert to a tableau.
        """
    def __repr__(self) -> str:
        ...
    def apply_gate_at_end(self, type: pytket._tket.circuit.OpType, qbs: typing.Sequence[pytket._tket.unit_id.Qubit]) -> None:
        """
        Update the tableau according to adding a Clifford gate after the current unitary, i.e. updates :math:`U` to :math:`GU` for a gate :math:`G`.
        
        :param type: The :py:class:`OpType` of the gate to add. Must be an unparameterised Clifford gate type.
        :param qbs: The qubits to apply the gate to. Length must match the arity of the given gate type.
        """
    def apply_gate_at_front(self, type: pytket._tket.circuit.OpType, qbs: typing.Sequence[pytket._tket.unit_id.Qubit]) -> None:
        """
        Update the tableau according to adding a Clifford gate before the current unitary, i.e. updates :math:`U` to :math:`UG` for a gate :math:`G`.
        
        :param type: The :py:class:`OpType` of the gate to add. Must be an unparameterised Clifford gate type.
        :param qbs: The qubits to apply the gate to. Length must match the arity of the given gate type.
        """
    def get_row_product(self, paulis: pytket._tket.pauli.QubitPauliTensor) -> pytket._tket.pauli.QubitPauliTensor:
        """
        Combine rows to yield the effect of a given Pauli string.
        
        :param paulis: The Pauli string :math:`P` to consider at the input.
        :return: The Pauli string :math:`Q` such that :math:`QU=UP`.
        """
    def get_xrow(self, qb: pytket._tket.unit_id.Qubit) -> pytket._tket.pauli.QubitPauliTensor:
        """
        Read off an X row as a Pauli string.
        
        :param qb: The qubits whose X row to read off.
        :return: The Pauli string :math:`P` such that :math:`PU=UX_{qb}`.
        """
    def get_zrow(self, qb: pytket._tket.unit_id.Qubit) -> pytket._tket.pauli.QubitPauliTensor:
        """
        Read off an Z row as a Pauli string.
        
        :param qb: The qubits whose Z row to read off.
        :return: The Pauli string :math:`P` such that :math:`PU=UZ_{qb}`.
        """
    def to_circuit(self) -> pytket._tket.circuit.Circuit:
        """
        Synthesises a unitary :py:class:`Circuit` realising the same unitary as the tableau. Uses the method from Aaronson & Gottesman: "Improved Simulation of Stabilizer Circuits", Theorem 8. This is not optimised for gate count, so is not recommended for performance-sensitive usage.
        """
class UnitaryTableauBox(pytket._tket.circuit.Op):
    """
    A Clifford unitary specified by its actions on Paulis.
    """
    @typing.overload
    def __init__(self, tab: UnitaryTableau) -> None:
        """
        Construct from a given tableau.
        
        :param tab: The :py:class:`UnitaryTableau` representing the desired unitary.
        """
    @typing.overload
    def __init__(self, xx: NDArray[numpy.bool_], xz: NDArray[numpy.bool_], xph: NDArray[numpy.bool_], zx: NDArray[numpy.bool_], zz: NDArray[numpy.bool_], zph: NDArray[numpy.bool_]) -> None:
        """
        Construct the tableau from the binary tables of its components.
        
        :param xx: The X component of the X rows.
        :param xz: The Z component of the X rows.
        :param xph: The phases of the X rows.
        :param zx: The X component of the Z rows.
        :param zz: The Z component of the Z rows.
        :param zph: The phases of the Z rows.
        """
    def get_circuit(self) -> pytket._tket.circuit.Circuit:
        """
        :return: The :py:class:`Circuit` described by the box.
        """
    def get_tableau(self) -> UnitaryTableau:
        """
        :return: The tableau representing the unitary operation.
        """
