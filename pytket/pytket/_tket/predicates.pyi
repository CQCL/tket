from __future__ import annotations
import pytket._tket.architecture
import pytket._tket.circuit
import pytket._tket.unit_id
import typing
__all__ = ['CliffordCircuitPredicate', 'CommutableMeasuresPredicate', 'CompilationUnit', 'ConnectivityPredicate', 'DefaultRegisterPredicate', 'DirectednessPredicate', 'GateSetPredicate', 'MaxNClRegPredicate', 'MaxNQubitsPredicate', 'MaxTwoQubitGatesPredicate', 'NoBarriersPredicate', 'NoClassicalBitsPredicate', 'NoClassicalControlPredicate', 'NoFastFeedforwardPredicate', 'NoMidMeasurePredicate', 'NoSymbolsPredicate', 'NoWireSwapsPredicate', 'NormalisedTK2Predicate', 'PlacementPredicate', 'Predicate', 'UserDefinedPredicate']
class CliffordCircuitPredicate(Predicate):
    """
    Predicate asserting that a circuit has only Clifford gates.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class CommutableMeasuresPredicate(Predicate):
    """
    Predicate asserting that all measurements can be delayed to the end of the circuit.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class CompilationUnit:
    """
    This class comprises a circuit and the predicates that the circuit is required to satisfy, for example to run on a backend.
    """
    @typing.overload
    def __init__(self, circuit: pytket._tket.circuit.Circuit) -> None:
        """
        Construct from a circuit, with no predicates.
        """
    @typing.overload
    def __init__(self, circuit: pytket._tket.circuit.Circuit, predicates: typing.Sequence[Predicate]) -> None:
        """
        Construct from a circuit and some required predicates.
        """
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    def check_all_predicates(self) -> bool:
        """
        :return: True if all predicates are satisfied, else False
        """
    @property
    def circuit(self) -> pytket._tket.circuit.Circuit:
        """
        Return a copy of the circuit.
        """
    @property
    def final_map(self) -> dict[pytket._tket.unit_id.UnitID, pytket._tket.unit_id.UnitID]:
        """
        Returns the map from the original qubits to their corresponding qubits at the end of the current circuit.
        """
    @property
    def initial_map(self) -> dict[pytket._tket.unit_id.UnitID, pytket._tket.unit_id.UnitID]:
        """
        Returns the map from the original qubits to the corresponding qubits at the start of the current circuit.
        """
class ConnectivityPredicate(Predicate):
    """
    Predicate asserting that a circuit satisfies a given connectivity graph. The graph is always considered to be undirected.
    """
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None:
        """
        Construct from an :py:class:`Architecture`.
        """
class DefaultRegisterPredicate(Predicate):
    """
    Predicate asserting that a circuit only uses the default quantum and classical registers.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class DirectednessPredicate(Predicate):
    """
    Predicate asserting that a circuit satisfies a given connectivity graph. The graph is always considered to be directed.
    """
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None:
        """
        Construct from an :py:class:`Architecture`.
        """
class GateSetPredicate(Predicate):
    """
    Predicate asserting that all operations are in the specified set of types.
    
    Note that the following are always permitted and do not need to be included in the specified set:
    
    - 'meta' operations (inputs, outputs, barriers);
    - ``OpType.Phase`` gates (which have no input or output wires).
    
    Classically conditioned operations are permitted provided that the conditioned operation is of a permitted type.
    """
    def __init__(self, allowed_types: set[pytket._tket.circuit.OpType]) -> None:
        """
        Construct from a set of gate types.
        """
    @property
    def gate_set(self) -> set[pytket._tket.circuit.OpType]:
        ...
class MaxNClRegPredicate(Predicate):
    """
    Predicate asserting that a circuit has at most n classical registers.
    """
    def __init__(self, arg0: int) -> None:
        """
        Constructor.
        """
class MaxNQubitsPredicate(Predicate):
    """
    Predicate asserting that a circuit has at most n qubits.
    """
    def __init__(self, arg0: int) -> None:
        """
        Constructor.
        """
class MaxTwoQubitGatesPredicate(Predicate):
    """
    Predicate asserting that a circuit has no gates with more than two input wires.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoBarriersPredicate(Predicate):
    """
    Predicate asserting that a circuit contains no Barrier operations.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoClassicalBitsPredicate(Predicate):
    """
    Predicate asserting that a circuit has no classical wires.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoClassicalControlPredicate(Predicate):
    """
    Predicate asserting that a circuit has no classical controls.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoFastFeedforwardPredicate(Predicate):
    """
    Predicate asserting that a circuit has no fast feedforward.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoMidMeasurePredicate(Predicate):
    """
    Predicate asserting that all measurements occur at the end of the circuit.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoSymbolsPredicate(Predicate):
    """
    Predicate asserting that no gates in the circuit have symbolic parameters.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NoWireSwapsPredicate(Predicate):
    """
    Predicate asserting that a circuit has no wire swaps.
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class NormalisedTK2Predicate(Predicate):
    """
    Asserts that all TK2 gates are normalised
    
    A gate TK2(a, b, c) is considered normalised if
    
     - If all expressions are non symbolic, then it must hold `0.5 ≥ a ≥ b ≥ |c|`.
     - In the ordering (a, b, c), any symbolic expression must appear before non-symbolic ones. The remaining non-symbolic expressions must still be ordered in non-increasing order and must be in the interval [0, 1/2], with the exception of the last one that may be in [-1/2, 1/2].
    """
    def __init__(self) -> None:
        """
        Constructor.
        """
class PlacementPredicate(Predicate):
    """
    Predicate asserting that a circuit has been acted on by some Placement object.
    """
    @typing.overload
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None:
        """
        Construct from an :py:class:`Architecture`.
        """
    @typing.overload
    def __init__(self, nodes: set[pytket._tket.unit_id.Node]) -> None:
        """
        Construct from a set of Node.
        """
class Predicate:
    """
    A predicate that may be satisfied by a circuit.
    """
    @staticmethod
    def from_dict(arg0: dict) -> Predicate:
        """
        Construct Predicate instance from JSON serializable dict representation of the Predicate.
        """
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    def implies(self, other: Predicate) -> bool:
        """
        :return: True if predicate implies another one, else False
        """
    def to_dict(self) -> dict:
        """
        Return a JSON serializable dict representation of the Predicate.
        
        :return: dict representation of the Predicate.
        """
    def verify(self, circuit: pytket._tket.circuit.Circuit) -> bool:
        """
        :return: True if circuit satisfies predicate, else False
        """
class UserDefinedPredicate(Predicate):
    """
    User-defined predicate.
    """
    def __init__(self, check_function: typing.Callable[[pytket._tket.circuit.Circuit], bool]) -> None:
        """
        Construct from a user-defined function from :py:class:`Circuit` to `bool`.
        """
