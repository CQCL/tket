from typing import Any
from __future__ import annotations
import pytket._tket.circuit
import pytket._tket.pauli
import typing
__all__ = ['GraphColourMethod', 'MeasurementBitMap', 'MeasurementSetup', 'PauliPartitionStrat', 'measurement_reduction', 'term_sequence']
class GraphColourMethod:
    """
    Enum for available methods to perform graph colouring.
    
    Members:
    
      Lazy : Does not build the graph before performing the colouring; partitions while iterating through the Pauli tensors in the input order.
    
      LargestFirst : Builds the graph and then greedily colours by iterating through the vertices, with the highest degree first.
    
      Exhaustive : Builds the graph and then systematically checks all possibilities until it finds a colouring with the minimum possible number of colours. Such colourings need not be unique. Exponential time in the worst case, but often runs much faster.
    """
    Exhaustive: typing.ClassVar[GraphColourMethod]  # value = <GraphColourMethod.Exhaustive: 2>
    LargestFirst: typing.ClassVar[GraphColourMethod]  # value = <GraphColourMethod.LargestFirst: 1>
    Lazy: typing.ClassVar[GraphColourMethod]  # value = <GraphColourMethod.Lazy: 0>
    __members__: typing.ClassVar[dict[str, GraphColourMethod]]  # value = {'Lazy': <GraphColourMethod.Lazy: 0>, 'LargestFirst': <GraphColourMethod.LargestFirst: 1>, 'Exhaustive': <GraphColourMethod.Exhaustive: 2>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):  # type: ignore
        ...
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
class MeasurementBitMap:
    """
    Maps Pauli tensors to Clifford circuit indices and bits required for measurement. A MeasurementBitMap belongs to a MeasurementSetup object, and dictates which bits are to be included in the measurement. As Clifford circuits may flip the parity of the corresponding Pauli tensor, the MeasurementBitMap optionally inverts the result.
    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):  # type: ignore
        ...
    @staticmethod
    def from_dict(arg0: dict) -> MeasurementBitMap:
        """
        Construct MeasurementBitMap instance from dict representation.
        """
    def __init__(self, circ_index: int, bits: typing.Sequence[int], invert: bool = False) -> None:
        """
        Constructs a MeasurementBitMap for some Clifford circuit index and bits, with an option to invert the result.
        
        :param circ_index: which measurement circuit the measurement map refers to
        :param bits: which bits are included in the measurement
        :param invert: whether to flip the parity of the result
        """
    def __repr__(self) -> str:
        ...
    def to_dict(self) -> dict:
        """
        JSON-serializable dict representation of the MeasurementBitMap.
        
        :return: dict representation of the MeasurementBitMap
        """
    @property
    def bits(self) -> list[int]:
        """
        Bits to measure
        """
    @property
    def circ_index(self) -> int:
        """
        Clifford circuit index
        """
    @property
    def invert(self) -> bool:
        """
        Whether result is inverted or not
        """
class MeasurementSetup:
    """
    Encapsulates an experiment in which the expectation value of an operator is to be measured via decomposition into QubitPauliStrings. Each tensor expectation value can be measured using shots. These values are then summed together with some weights to retrieve the desired operator expctation value.
    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):  # type: ignore
        ...
    @staticmethod
    def from_dict(arg0: dict) -> MeasurementSetup:
        """
        Construct MeasurementSetup instance from dict representation.
        """
    def __init__(self) -> None:
        """
        Constructs an empty MeasurementSetup object
        """
    def __repr__(self) -> str:
        ...
    def add_measurement_circuit(self, circ: pytket._tket.circuit.Circuit) -> None:
        """
        Add a Clifford circuit that rotates into some Pauli basis
        """
    def add_result_for_term(self, term: pytket._tket.pauli.QubitPauliString, result: MeasurementBitMap) -> None:
        """
        Add a new Pauli string with a corresponding BitMap
        """
    def to_dict(self) -> dict:
        """
        JSON-serializable dict representation of the MeasurementSetup.
        
        :return: dict representation of the MeasurementSetup
        """
    def verify(self) -> bool:
        """
        Checks that the strings to be measured correspond to the correct strings generated by the measurement circs. Checks for parity by comparing to the `invert` flag.
        
        :return: True or False
        """
    @property
    def measurement_circs(self) -> list[pytket._tket.circuit.Circuit]:
        """
        Clifford measurement circuits.
        """
    @property
    def results(self) -> dict[pytket._tket.pauli.QubitPauliString, list[MeasurementBitMap]]:
        """
        Map from Pauli strings to MeasurementBitMaps
        """
class PauliPartitionStrat:
    """
    Enum for available strategies to partition Pauli tensors.
    
    Members:
    
      NonConflictingSets : Build sets of Pauli tensors in which each qubit has the same Pauli or Pauli.I. Requires no additional CX gates for diagonalisation.
    
      CommutingSets : Build sets of mutually commuting Pauli tensors. Requires O(n^2) CX gates to diagonalise.
    """
    CommutingSets: typing.ClassVar[PauliPartitionStrat]  # value = <PauliPartitionStrat.CommutingSets: 1>
    NonConflictingSets: typing.ClassVar[PauliPartitionStrat]  # value = <PauliPartitionStrat.NonConflictingSets: 0>
    __members__: typing.ClassVar[dict[str, PauliPartitionStrat]]  # value = {'NonConflictingSets': <PauliPartitionStrat.NonConflictingSets: 0>, 'CommutingSets': <PauliPartitionStrat.CommutingSets: 1>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):  # type: ignore
        ...
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
def measurement_reduction(strings: typing.Sequence[pytket._tket.pauli.QubitPauliString], strat: PauliPartitionStrat, method: GraphColourMethod = GraphColourMethod.Lazy, cx_config: pytket._tket.circuit.CXConfigType = pytket._tket.circuit.CXConfigType.Snake) -> MeasurementSetup:
    """
    Automatically performs graph colouring and diagonalisation to reduce measurements required for Pauli strings.
    
    :param strings: A list of `QubitPauliString` objects to be partitioned.
    :param strat: The `PauliPartitionStrat` to use.
    :param method: The `GraphColourMethod` to use.
    :param cx_config: Whenever diagonalisation is required, use this configuration of CX gates
    :return: a :py:class:`MeasurementSetup` object
    """
def term_sequence(strings: typing.Sequence[pytket._tket.pauli.QubitPauliString], strat: PauliPartitionStrat = PauliPartitionStrat.CommutingSets, method: GraphColourMethod = GraphColourMethod.Lazy) -> list[list[pytket._tket.pauli.QubitPauliString]]:
    """
    Takes in a list of QubitPauliString objects and partitions them into mutually commuting sets according to some PauliPartitionStrat, then sequences in an arbitrary order.
    
    :param tensors: A list of `QubitPauliString` objects to be sequenced. Assumes that each Pauli tensor is unique, and does not combine equivalent tensors.
    :param strat: The `PauliPartitionStrat` to use. Defaults to `CommutingSets`.
    :param method: The `GraphColourMethod` to use.
    :return: a list of lists of :py:class:`QubitPauliString` s
    """
