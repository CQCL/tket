from typing import Any
from __future__ import annotations
import pytket._tket.architecture
import pytket._tket.circuit
import pytket._tket.unit_id
import typing
__all__ = ['GraphPlacement', 'LinePlacement', 'NoiseAwarePlacement', 'Placement', 'place_fully_connected', 'place_with_map']
class GraphPlacement(Placement):
    """
    The GraphPlacement class, contains methods for getting maps between Circuit Qubits and Architecture Nodes and for relabelling Circuit Qubits.
    """
    def __init__(self, arc: pytket._tket.architecture.Architecture, maximum_matches: int = 1000, timeout: int = 1000, maximum_pattern_gates: int = 100, maximum_pattern_depth: int = 100) -> None:
        """
        The constructor for a GraphPlacement object. The Architecture object describes the connectivity between qubits. To find a qubit to node assignment, this method constructs a pattern graph where vertices are Circuit qubits and edges mean a pair of qubits have an interaction in the circuit, and then tries to find a weighted subgraph monomorphsim to the architecture connectivity, or target, graph. Edges in the pattern graph are weighted by the circuit depth at which the interaction between a pair of qubit occurs. The number of edges added to the pattern graph is effected by the maximum_pattern_gates and maximum_pattern_depth arguments. If no subgraph monomorphism can be found, lower edge weights are removed from the pattern graph, are more edges are added to the target graph. Edges added to the pattern graph are weighted lower to reflect what the distance between the Nodes they  are added between was on the original target graph. 
        
        :param arc: An Architecture object.
        :param maximum_matches: The total number of weighted subgraph monomorphisms that can be found before matches are returned.
        :param timeout: Total time in milliseconds before stopping search for monomorphisms.
        :param maximum_pattern_gates: The upper bound on the number of circuit gates used to construct the pattern graph for finding subgraph monomorphisms.
        :param maximum_pattern_depth: The upper bound on the circuit depth gates are added to the pattern graph to for finding subgraph monomorphisms.
        """
    def __repr__(self: Placement) -> str:
        ...
    def modify_config(self, **kwargs: Any) -> None:
        """
        Deprecated and no longer modifies parameters for finding solutions. Please create a new GraphPlacement object instead
        """
class LinePlacement(Placement):
    """
    The LinePlacement class, contains methods for getting maps between Circuit Qubits and Architecture Nodes and for relabelling Circuit Qubits.
    """
    def __init__(self, arc: pytket._tket.architecture.Architecture, maximum_line_gates: int = 100, maximum_line_depth: int = 100) -> None:
        """
        The constructor for a LinePlacement object. The Architecture object describes the connectivity between qubits. In this class, a reduced qubit interaction subgraph is constructed where each node has maximum outdegree 2 and does not construct a circle (i.e. lines). To place the Circuit, a Hamiltonian Path is found in the Architecture and this subgraph of lines is assigned along it.
        
        :param arc: An Architecture object.
        :param maximum_line_gates: maximum number of gates in the circuit considered when constructing lines for assigning to the graph
        :param maximum_line_depth: maximum depth of circuit considered when constructing lines for assigning to the graph
        """
    def __repr__(self: Placement) -> str:
        ...
class NoiseAwarePlacement(Placement):
    """
    The NoiseAwarePlacement class, contains methods for getting maps between Circuit Qubits and Architecture Nodes and for relabelling Circuit Qubits. It uses gate error rates and readout errors to find the best placement map.
    """
    def __init__(self, arc: pytket._tket.architecture.Architecture, node_errors: dict[pytket._tket.unit_id.Node, float] = {}, link_errors: dict[tuple[pytket._tket.unit_id.Node, pytket._tket.unit_id.Node], float] = {}, readout_errors: dict[pytket._tket.unit_id.Node, float] = {}, maximum_matches: int = 1000, timeout: int = 1000, maximum_pattern_gates: int = 100, maximum_pattern_depth: int = 100) -> None:
        """
        The constructor for a NoiseAwarePlacement object. The Architecture object describes the connectivity between qubits. The dictionaries passed as parameters indicate the average gate errors for single- and two-qubit gates as well as readouterrors.  If no error is given for a given node or pair of nodes,the fidelity is assumed to be 1.
        
        :param arc: An Architecture object
        :param node_errors: a dictionary mapping nodes in the architecture to average single-qubit gate errors
        :param link_errors: a dictionary mapping pairs of nodes in the architecture to average two-qubit gate errors
        :param readout_errors: a dictionary mapping nodes in the architecture to average measurement readout errors.
        :param maximum_matches: The total number of weighted subgraph monomorphisms that can be found before matches are returned.
        :param timeout: Total time in milliseconds before stopping search for monomorphisms.
        :param maximum_pattern_gates: The upper bound on the number of circuit gates used to construct the pattern graph for finding subgraph monomorphisms.
        :param maximum_pattern_depth: The upper bound on the circuit depth gates are added to the pattern graph to for finding subgraph monomorphisms.
        """
    def __repr__(self: Placement) -> str:
        ...
    def modify_config(self, **kwargs: Any) -> None:
        """
        Deprecated and no longer modifies paramters for finding solutions. Please create a new NoiseAwarePlacement object instead
        """
class Placement:
    """
    The base Placement class, contains methods for getting maps between Circuit Qubits and Architecture Nodes and for relabelling Circuit Qubits.
    """
    @staticmethod
    def from_dict(arg0: dict) -> Placement:
        """
        Construct Placement instance from JSON serializable dict representation of the Placement.
        """
    @staticmethod
    def place_with_map(circuit: pytket._tket.circuit.Circuit, qmap: dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node]) -> bool:
        """
        Relabels Circuit Qubits to Architecture Nodes using given map. 
        
        :param circuit: The circuit being relabelled
        :param qmap: The map from logical to physical qubits to apply.
        """
    def __init__(self, arc: pytket._tket.architecture.Architecture) -> None:
        """
        The constructor for a Placement object. The Architecture object describes the connectivity between qubits.
        
        :param arc: An Architecture object.
        """
    def __repr__(self) -> str:
        ...
    def get_placement_map(self, circuit: pytket._tket.circuit.Circuit) -> dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node]:
        """
        Returns a map from logical to physical qubits that is Architecture appropriate for the given Circuit. 
        
        :param circuit: The circuit a map is designed for.
        :return: dictionary mapping :py:class:`Qubit` s to :py:class:`Node` s
        """
    def get_placement_maps(self, circuit: pytket._tket.circuit.Circuit, matches: int = 100) -> list[dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node]]:
        """
        Returns a list of maps from logical to physical qubits that are Architecture appropriate for the given Circuit. Each map is estimated to given a similar SWAP overheard after routing. 
        
        :param circuit: The circuit the maps are designed for.
        :param matches: The maximum number of maps returned by the method.
        :return: list of dictionaries mapping :py:class:`Qubit` s to :py:class:`Node` s
        """
    def place(self, circuit: pytket._tket.circuit.Circuit) -> bool:
        """
        Relabels Circuit Qubits to Architecture Nodes and 'unplaced'. For base Placement, all Qubits and labelled 'unplaced'. 
        
        :param circuit: The Circuit being relabelled.
        """
    def to_dict(self) -> typing.Any:
        """
        Return a JSON serializable dict representation of the Placement.
        
        :return: dict representing the Placement.
        """
def place_fully_connected(circuit: pytket._tket.circuit.Circuit, fully_connected: pytket._tket.architecture.FullyConnected) -> None:
    """
    Relabels all Circuit Qubits to the Node objects of a FullyConnected object. 
    
    :param circuit: The Circuit being relabelled. 
    :param fully_connected: FullyConnected object Qubits being relabelled to match.
    """
def place_with_map(circuit: pytket._tket.circuit.Circuit, qmap: dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node]) -> None:
    """
    Relabels Circuit Qubits according to a map. If provided map is partial, remaining Circuit Qubits are left 'unplaced'. 
    
    :param circuit: The Circuit being relabelled. 
    :param qmap: The map from logical to physical qubits to apply.
    """
