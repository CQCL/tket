from typing import ClassVar, Dict, List

import pytket._tket.circuit
import pytket._tket.pauli
from .type_helpers import json

class GraphColourMethod:
    __members__: ClassVar[dict] = ...  # read-only
    Exhaustive: ClassVar[GraphColourMethod] = ...
    LargestFirst: ClassVar[GraphColourMethod] = ...
    Lazy: ClassVar[GraphColourMethod] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MeasurementBitMap:
    def __init__(self, circ_index: int, bits: List[int], invert: bool = ...) -> None: ...
    @classmethod
    def from_dict(cls, arg0: json) -> MeasurementBitMap: ...
    def to_dict(self) -> json: ...
    @property
    def bits(self) -> List[int]: ...
    @property
    def circ_index(self) -> int: ...
    @property
    def invert(self) -> bool: ...

class MeasurementSetup:
    def __init__(self) -> None: ...
    def add_measurement_circuit(self, circ: pytket._tket.circuit.Circuit) -> None: ...
    def add_result_for_term(self, term: pytket._tket.pauli.QubitPauliString, result: MeasurementBitMap) -> None: ...
    @classmethod
    def from_dict(cls, arg0: json) -> MeasurementSetup: ...
    def to_dict(self) -> json: ...
    def verify(self) -> bool: ...
    @property
    def measurement_circs(self) -> List[pytket._tket.circuit.Circuit]: ...
    @property
    def results(self) -> Dict[pytket._tket.pauli.QubitPauliString,List[MeasurementBitMap]]: ...

class PauliPartitionStrat:
    __members__: ClassVar[dict] = ...  # read-only
    CommutingSets: ClassVar[PauliPartitionStrat] = ...
    NonConflictingSets: ClassVar[PauliPartitionStrat] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

def measurement_reduction(strings: List[pytket._tket.pauli.QubitPauliString], strat: PauliPartitionStrat, method: GraphColourMethod = ..., cx_config: pytket._tket.circuit.CXConfigType = ...) -> MeasurementSetup: ...
def term_sequence(strings: List[pytket._tket.pauli.QubitPauliString], strat: PauliPartitionStrat = ..., method: GraphColourMethod = ...) -> List[List[pytket._tket.pauli.QubitPauliString]]: ...
