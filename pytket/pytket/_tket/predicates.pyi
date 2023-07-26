from typing import Callable, Dict, List, Set, Any

from typing import overload
import pytket._tket.architecture
import pytket._tket.circuit
from .type_helpers import json

class CliffordCircuitPredicate(Predicate):
    def __init__(self) -> None: ...

class CommutableMeasuresPredicate(Predicate):
    def __init__(self) -> None: ...

class CompilationUnit:
    @overload
    def __init__(self, circuit: pytket._tket.circuit.Circuit) -> None: ...
    @overload
    def __init__(self, circuit: pytket._tket.circuit.Circuit, predicates: List[Predicate]) -> None: ...
    def check_all_predicates(self) -> bool: ...
    @property
    def circuit(self) -> pytket._tket.circuit.Circuit: ...
    @property
    def final_map(self) -> Dict[pytket._tket.circuit.UnitID,pytket._tket.circuit.UnitID]: ...
    @property
    def initial_map(self) -> Dict[pytket._tket.circuit.UnitID,pytket._tket.circuit.UnitID]: ...

class ConnectivityPredicate(Predicate):
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None: ...

class DefaultRegisterPredicate(Predicate):
    def __init__(self) -> None: ...

class DirectednessPredicate(Predicate):
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None: ...

class GateSetPredicate(Predicate):
    def __init__(self, allowed_types: Set[pytket._tket.circuit.OpType]) -> None: ...
    @property
    def gate_set(self) -> Set[pytket._tket.circuit.OpType]: ...

class MaxNClRegPredicate(Predicate):
    def __init__(self, arg0: int) -> None: ...

class MaxNQubitsPredicate(Predicate):
    def __init__(self, arg0: int) -> None: ...

class MaxTwoQubitGatesPredicate(Predicate):
    def __init__(self) -> None: ...

class NoBarriersPredicate(Predicate):
    def __init__(self) -> None: ...

class NoClassicalBitsPredicate(Predicate):
    def __init__(self) -> None: ...

class NoClassicalControlPredicate(Predicate):
    def __init__(self) -> None: ...

class NoFastFeedforwardPredicate(Predicate):
    def __init__(self) -> None: ...

class NoMidMeasurePredicate(Predicate):
    def __init__(self) -> None: ...

class NoSymbolsPredicate(Predicate):
    def __init__(self) -> None: ...

class NoWireSwapsPredicate(Predicate):
    def __init__(self) -> None: ...

class NormalisedTK2Predicate(Predicate):
    def __init__(self) -> None: ...

class PlacementPredicate(Predicate):
    @overload
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None: ...
    @overload
    def __init__(self, nodes: Set[pytket._tket.circuit.Node]) -> None: ...

class Predicate:
    def __init__(self, *args, **kwargs: Any) -> None: ...
    @classmethod
    def from_dict(cls, arg0: json) -> Predicate: ...
    def implies(self, other: Predicate) -> bool: ...
    def to_dict(self) -> json: ...
    def verify(self, circuit: pytket._tket.circuit.Circuit) -> bool: ...

class UserDefinedPredicate(Predicate):
    def __init__(self, check_function: Callable[[pytket._tket.circuit.Circuit],bool]) -> None: ...
