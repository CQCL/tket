from __future__ import annotations
import pytket._tket.unit_id
import typing
__all__ = ['Architecture', 'FullyConnected', 'RingArch', 'SquareGrid']
class Architecture:
    """
    Class describing the connectivity of qubits on a general device.
    """
    @staticmethod
    def from_dict(arg0: dict) -> Architecture:
        """
        Construct Architecture instance from JSON serializable dict representation of the Architecture.
        """
    def __deepcopy__(self, arg0: dict) -> Architecture:
        ...
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __hash__(self) -> int:
        """
        Hashing is not implemented for this class, attempting to hash an object will raise a type error
        """
    @typing.overload
    def __init__(self) -> None:
        """
        Produces an empty architecture
        """
    @typing.overload
    def __init__(self, connections: typing.Sequence[tuple[int, int]]) -> None:
        """
        The constructor for an architecture with connectivity between qubits.
        
        :param connections: A list of pairs representing qubit indices that can perform two-qubit operations
        """
    @typing.overload
    def __init__(self, connections: typing.Sequence[tuple[pytket._tket.unit_id.Node, pytket._tket.unit_id.Node]]) -> None:
        """
        The constructor for an architecture with connectivity between qubits.
        
        :param connections: A list of pairs representing Nodes that can perform two-qubit operations
        """
    def __repr__(self) -> str:
        ...
    def get_adjacent_nodes(self, node: pytket._tket.unit_id.Node) -> set[pytket._tket.unit_id.Node]:
        """
        given a node, returns adjacent nodes in Architecture.
        """
    def get_distance(self, node_0: pytket._tket.unit_id.Node, node_1: pytket._tket.unit_id.Node) -> int:
        """
        given two nodes in Architecture, returns distance between them
        """
    def to_dict(self) -> dict:
        """
        Return a JSON serializable dict representation of the Architecture.
        
        :return: dict containing nodes and links.
        """
    def valid_operation(self, uids: typing.Sequence[pytket._tket.unit_id.Node]) -> bool:
        """
        nodes can be executed on the Architecture connectivity graph.
        
        :param uids: list of UnitIDs validity is being checked for
        """
    @property
    def coupling(self) -> list[tuple[pytket._tket.unit_id.Node, pytket._tket.unit_id.Node]]:
        """
        Returns the coupling map of the Architecture as UnitIDs. 
        """
    @property
    def nodes(self) -> list[pytket._tket.unit_id.Node]:
        """
        Returns all nodes of architecture as Node objects.
        """
class FullyConnected:
    """
    A specialised non-Architecture object emulating an architecture with all qubits connected. Not compatible with Routing or Placement methods.
    """
    @staticmethod
    def from_dict(arg0: dict) -> FullyConnected:
        """
        Construct FullyConnected instance from dict representation.
        """
    def __deepcopy__(self, arg0: dict) -> FullyConnected:
        ...
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __hash__(self) -> int:
        """
        Hashing is not implemented for this class, attempting to hash an object will raise a type error
        """
    def __init__(self, n: int, label: str = 'fcNode') -> None:
        """
        Construct a fully-connected architecture.
        
        :param n: number of qubits
        :param label: Name for Node in FullyConnected Architecture
        """
    def __repr__(self) -> str:
        ...
    def to_dict(self) -> dict:
        """
        JSON-serializable dict representation of the architecture.
        
        :return: dict containing nodes
        """
    @property
    def nodes(self) -> list[pytket._tket.unit_id.Node]:
        """
        All nodes of the architecture as :py:class:`Node` objects.
        """
class RingArch(Architecture):
    """
    Inherited Architecture class for number of qubits arranged in a ring.
    """
    def __deepcopy__(self, arg0: dict) -> RingArch:
        ...
    def __init__(self, nodes: int, label: str = 'ringNode') -> None:
        """
        The constructor for a RingArchitecture with some undirected connectivity between qubits.
        
        :param number of qubits:
        :param label: Name for Node in RingArch Architecture
        """
    def __repr__(self) -> str:
        ...
class SquareGrid(Architecture):
    """
    Inherited Architecture class for qubits arranged in a square lattice of given number of rows and columns. Qubits are arranged with qubits values increasing first along rows then along columns i.e. for a 3 x 3 grid:
    
     0 1 2
    
     3 4 5
    
     6 7 8
    """
    def __deepcopy__(self, arg0: dict) -> SquareGrid:
        ...
    @typing.overload
    def __init__(self, n_rows: int, n_columns: int, label: str = 'gridNode') -> None:
        """
        The constructor for a Square Grid architecture with some undirected connectivity between qubits.
        
        :param n_rows: The number of rows in the grid
        :param n_columns: The number of columns in the grid
        :param label: Name for Node in SquareGrid Architecture
        """
    @typing.overload
    def __init__(self, n_rows: int, n_columns: int, n_layers: int = 1, label: str = 'gridNode') -> None:
        """
        The constructor for  a Square Grid architecture with some undirected connectivity between qubits.
        
        :param n_rows: The number of rows in the grid
        :param n_columns: The number of columns in the grid
        :param n_layers: The number of layers of grids
        :param label: Name for Node in SquareGrid Architecture
        """
    def __repr__(self) -> str:
        ...
    def qind_to_squind(self, index: int) -> tuple[int, int]:
        """
        Converts a single qubit index to a (row,column) index for a square grid.
        
        :param index: The global qubit index
        :return: the corresponding grid index as a pair (row,column)
        """
    def squind_to_qind(self, row: int, column: int) -> int:
        """
        Converts a (row,column) index for a square grid to a single qubit index
        
        :param row: The given row index
        :param column: The given column index
        :return: the corresponding global qubit index
        """
