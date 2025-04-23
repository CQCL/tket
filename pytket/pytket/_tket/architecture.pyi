from collections.abc import Sequence
from typing import overload

import pytket._tket.unit_id


class Architecture:
    """Class describing the connectivity of qubits on a general device."""

    @overload
    def __init__(self) -> None:
        """Produces an empty architecture"""

    @overload
    def __init__(self, connections: Sequence[tuple[int, int]]) -> None:
        """
        The constructor for an architecture with connectivity between qubits.

        :param connections: A list of pairs representing qubit indices that can perform two-qubit operations
        """

    @overload
    def __init__(self, connections: Sequence[tuple[pytket._tket.unit_id.Node, pytket._tket.unit_id.Node]]) -> None:
        """
        The constructor for an architecture with connectivity between qubits.

        :param connections: A list of pairs representing Nodes that can perform two-qubit operations
        """

    def __repr__(self) -> str: ...

    def get_distance(self, node_0: pytket._tket.unit_id.Node, node_1: pytket._tket.unit_id.Node) -> int:
        """given two nodes in Architecture, returns distance between them"""

    def valid_operation(self, uids: Sequence[pytket._tket.unit_id.Node], bidirectional: bool = True) -> bool:
        """
        Checks if a given operation on the given nodes can be executed on the Architecture's connectivity graph.
        The operation is considered valid if:

        - The operation acts on a single node that belongs to the Architecture.
        - The operation acts on two nodes, and either:

          - `bidirectional` is True and an edge exists between the two nodes in either direction.
          - `bidirectional` is False and an edge exists from uids[0] to uids[1].

        The function always returns False if the number of nodes exceeds 2.

        :param uids: list of UnitIDs validity is being checked for.
        :param bidirectional: If True, treats edges in the coupling graph as bidirectional. Defaults to True.
        """

    def get_adjacent_nodes(self, node: pytket._tket.unit_id.Node) -> set[pytket._tket.unit_id.Node]:
        """given a node, returns adjacent nodes in Architecture."""

    @property
    def nodes(self) -> list[pytket._tket.unit_id.Node]:
        """Returns all nodes of architecture as Node objects."""

    @property
    def coupling(self) -> list[tuple[pytket._tket.unit_id.Node, pytket._tket.unit_id.Node]]:
        """Returns the coupling map of the Architecture as UnitIDs."""

    def to_dict(self) -> dict:
        """
        Return a JSON serializable dict representation of the Architecture.

        :return: dict containing nodes and links.
        """

    @staticmethod
    def from_dict(arg: dict, /) -> Architecture:
        """
        Construct Architecture instance from JSON serializable dict representation of the Architecture.
        """

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...

    def __deepcopy__(self, arg: dict, /) -> Architecture: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __hash__(self) -> int:
        """
        Hashing is not implemented for this class, attempting to hash an object will raise a type error
        """

class SquareGrid(Architecture):
    """
    Inherited Architecture class for qubits arranged in a square lattice of given number of rows and columns. Qubits are arranged with qubits values increasing first along rows then along columns i.e. for a 3 x 3 grid:

     0 1 2

     3 4 5

     6 7 8
    """

    @overload
    def __init__(self, n_rows: int, n_columns: int, label: str = 'gridNode') -> None:
        """
        The constructor for a Square Grid architecture with some undirected connectivity between qubits.

        :param n_rows: The number of rows in the grid
        :param n_columns: The number of columns in the grid
        :param label: Name for Node in SquareGrid Architecture
        """

    @overload
    def __init__(self, n_rows: int, n_columns: int, n_layers: int = 1, label: str = 'gridNode') -> None:
        """
        The constructor for  a Square Grid architecture with some undirected connectivity between qubits.

        :param n_rows: The number of rows in the grid
        :param n_columns: The number of columns in the grid
        :param n_layers: The number of layers of grids
        :param label: Name for Node in SquareGrid Architecture
        """

    def squind_to_qind(self, row: int, column: int) -> int:
        """
        Converts a (row,column) index for a square grid to a single qubit index

        :param row: The given row index
        :param column: The given column index
        :return: the corresponding global qubit index
        """

    def qind_to_squind(self, index: int) -> tuple[int, int]:
        """
        Converts a single qubit index to a (row,column) index for a square grid.

        :param index: The global qubit index
        :return: the corresponding grid index as a pair (row,column)
        """

    def __deepcopy__(self, arg: dict, /) -> SquareGrid: ...

    def __repr__(self) -> str: ...

class RingArch(Architecture):
    """Inherited Architecture class for number of qubits arranged in a ring."""

    def __init__(self, nodes: int, label: str = 'ringNode') -> None:
        """
        The constructor for a RingArchitecture with some undirected connectivity between qubits.

        :param number of qubits:
        :param label: Name for Node in RingArch Architecture
        """

    def __deepcopy__(self, arg: dict, /) -> RingArch: ...

    def __repr__(self) -> str: ...

class FullyConnected:
    """
    A specialised non-Architecture object emulating an architecture with all qubits connected. Not compatible with Routing or Placement methods.
    """

    def __init__(self, n: int, label: str = 'fcNode') -> None:
        """
        Construct a fully-connected architecture.

        :param n: number of qubits
        :param label: Name for Node in FullyConnected Architecture
        """

    def __deepcopy__(self, arg: dict, /) -> FullyConnected: ...

    def __repr__(self) -> str: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __hash__(self) -> int:
        """
        Hashing is not implemented for this class, attempting to hash an object will raise a type error
        """

    @property
    def nodes(self) -> list[pytket._tket.unit_id.Node]:
        """All nodes of the architecture as :py:class:`Node` objects."""

    def to_dict(self) -> dict:
        """
        JSON-serializable dict representation of the architecture.

        :return: dict containing nodes
        """

    @staticmethod
    def from_dict(arg: dict, /) -> FullyConnected:
        """Construct FullyConnected instance from dict representation."""

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...
