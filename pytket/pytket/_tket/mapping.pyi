from __future__ import annotations
import pytket._tket.architecture
import pytket._tket.circuit
import pytket._tket.unit_id
import typing
__all__ = ['AASLabellingMethod', 'AASRouteRoutingMethod', 'BoxDecompositionRoutingMethod', 'LexiLabellingMethod', 'LexiRouteRoutingMethod', 'MappingManager', 'MultiGateReorderRoutingMethod', 'RoutingMethod', 'RoutingMethodCircuit']
class AASLabellingMethod(RoutingMethod):
    """
    Defines a Labeling Method for aas for labelling all unplaced qubits in a circuit
    """
    def __init__(self) -> None:
        """
        AASLabellingMethod constructor.
        """
class AASRouteRoutingMethod(RoutingMethod):
    """
    Defines a RoutingMethod object for mapping circuits that uses the architecture aware synthesis method implemented in tket.
    """
    def __init__(self, aaslookahead: int) -> None:
        """
        AASRouteRoutingMethod constructor.
        
        :param aaslookahead: recursive interation depth of the architecture aware synthesis.method.
        """
class BoxDecompositionRoutingMethod(RoutingMethod):
    """
    Defines a RoutingMethod object for decomposing boxes.
    """
    def __init__(self) -> None:
        """
        BoxDecompositionRoutingMethod constructor.
        """
class LexiLabellingMethod(RoutingMethod):
    """
    Defines a RoutingMethod for labelling Qubits that uses the Lexicographical Comparison approach outlined in arXiv:1902.08091.
    """
    def __init__(self) -> None:
        """
        LexiLabellingMethod constructor.
        """
class LexiRouteRoutingMethod(RoutingMethod):
    """
    Defines a RoutingMethod object for mapping circuits that uses the Lexicographical Comparison approach outlined in arXiv:1902.08091.Only supports 1-qubit, 2-qubit and barrier gates.
    """
    def __init__(self, lookahead: int = 10) -> None:
        """
        LexiRoute constructor.
        
        :param lookahead: Maximum depth of lookahead employed when picking SWAP for purpose of logical to physical mapping.
        """
class MappingManager:
    """
    Defined by a pytket Architecture object, maps Circuit logical qubits to physically permitted Architecture qubits. Mapping is completed by sequential routing (full or partial) of subcircuits. A custom method for routing (full or partial) of subcircuits can be defined in Python.
    """
    def __init__(self, architecture: pytket._tket.architecture.Architecture) -> None:
        """
        MappingManager constructor.
        
        :param architecture: pytket Architecture object.
        """
    def route_circuit(self, circuit: pytket._tket.circuit.Circuit, routing_methods: typing.Sequence[RoutingMethod]) -> bool:
        """
        Maps from given logical circuit to physical circuit. Modification defined by route_subcircuit, but typically this proceeds by insertion of SWAP gates that permute logical qubits on physical qubits.
        
        :param circuit: pytket circuit to be mapped
        :param routing_methods: Ranked methods to use for routing subcircuits. In given order, each method is sequentially checked for viability, with the first viable method being used.
        """
class MultiGateReorderRoutingMethod(RoutingMethod):
    """
    Defines a RoutingMethod object for commuting physically permitted multi-qubit gates to the front of the subcircuit.
    """
    def __init__(self, max_depth: int = 10, max_size: int = 10) -> None:
        """
        MultiGateReorderRoutingMethod constructor.
        
        :param max_depth: Maximum number of layers of gates checked for simultaneous commutation. 
        :param max_size: Maximum number of gates checked for simultaneous commutation.
        """
class RoutingMethod:
    """
    Parent class for RoutingMethod, for inheritance purposes only, not for usage.
    """
    def __init__(self) -> None:
        ...
class RoutingMethodCircuit(RoutingMethod):
    """
    The RoutingMethod class captures a method for partially mapping logical subcircuits to physical operations as permitted by some architecture. Ranked RoutingMethod objects are used by the MappingManager to route whole circuits.
    """
    def __init__(self, route_subcircuit: typing.Callable[[pytket._tket.circuit.Circuit, pytket._tket.architecture.Architecture], tuple[bool, pytket._tket.circuit.Circuit, dict[pytket._tket.unit_id.UnitID, pytket._tket.unit_id.UnitID], dict[pytket._tket.unit_id.UnitID, pytket._tket.unit_id.UnitID]]], max_size: int, max_depth: int) -> None:
        """
        Constructor for a routing method defined by partially routing subcircuits.
        
        :param route_subcircuit: A function declaration that given a Circuit and Architecture object, returns a tuple containing a bool informing MappingManager whether to substitute the returned circuit into the circuit being routed, a new modified circuit, the initial logical to physical qubit mapping of the modified circuit and the permutation of logical to physical qubit mapping given operations in the modified circuit
        :param max_size: The maximum number of gates permitted in a subcircuit
        :param max_depth: The maximum permitted depth of a subcircuit.
        """
