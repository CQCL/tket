from collections.abc import Sequence
import enum
from typing import Union, overload

import pytket.circuit.logic_exp


_TEMP_REG_SIZE: int = 64

_TEMP_BIT_NAME: str = 'tk_SCRATCH_BIT'

_TEMP_BIT_REG_BASE: str = 'tk_SCRATCH_BITREG'

_DEBUG_ONE_REG_PREFIX: str = 'tk_DEBUG_ONE_REG'

_DEBUG_ZERO_REG_PREFIX: str = 'tk_DEBUG_ZERO_REG'

class UnitType(enum.Enum):
    """Enum for data types of units in circuits (e.g. Qubits vs Bits)."""

    qubit = 0
    """A single Qubit"""

    wasmstate = 2
    """A single WasmState"""

    bit = 1
    """A single classical Bit"""

class UnitID:
    """A handle to a computational unit (e.g. qubit, bit)"""

    def __init__(self) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __lt__(self, arg: UnitID, /) -> bool: ...

    def __repr__(self) -> str: ...

    def __hash__(self) -> int: ...

    def __copy__(self) -> UnitID: ...

    def __deepcopy__(self, arg: dict, /) -> UnitID: ...

    @property
    def reg_name(self) -> str:
        """Readable name of register"""

    @property
    def index(self) -> list[int]:
        """
        Index vector describing position in the register. The length of this vector is the dimension of the register
        """

    @property
    def type(self) -> UnitType:
        """
        Type of unit, either ``UnitType.qubit`` or ``UnitType.bit`` or ``UnitType.wasmstate``
        """

class Qubit(UnitID):
    """A handle to a qubit"""

    @overload
    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default qubit register

        :param index: The index in the register
        """

    @overload
    def __init__(self, name: str) -> None:
        """
        Constructs a named id (i.e. corresponding to a singleton register)

        :param name: The readable name for the id
        """

    @overload
    def __init__(self, name: str, index: int) -> None:
        """
        Constructs an indexed id (i.e. corresponding to an element in a linear register)

        :param name: The readable name for the register
        :param index: The numerical index
        """

    @overload
    def __init__(self, name: str, row: int, col: int) -> None:
        """
        Constructs a doubly-indexed id (i.e. corresponding to an element in a grid register)

        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        """

    @overload
    def __init__(self, name: str, index: Sequence[int]) -> None:
        """
        Constructs an id with an arbitrary-dimensional index

        :param name: The readable name for the register
        :param index: The index vector
        """

    def __copy__(self) -> Qubit: ...

    def __deepcopy__(self, arg: dict, /) -> Qubit: ...

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...

    def to_list(self) -> list:
        """:return: a JSON serializable list representation of the Qubit"""

    @staticmethod
    def from_list(arg: list, /) -> Qubit:
        """
        Construct Qubit instance from JSON serializable list representation of the Qubit.
        """

class Bit(UnitID):
    """A handle to a bit"""

    @overload
    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default classical register

        :param index: The index in the register
        """

    @overload
    def __init__(self, name: str) -> None:
        """
        Constructs a named id (i.e. corresponding to a singleton register)

        :param name: The readable name for the id
        """

    @overload
    def __init__(self, name: str, index: int) -> None:
        """
        Constructs an indexed id (i.e. corresponding to an element in a linear register)

        :param name: The readable name for the register
        :param index: The numerical index
        """

    @overload
    def __init__(self, name: str, row: int, col: int) -> None:
        """
        Constructs a doubly-indexed id (i.e. corresponding to an element in a grid register)

        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        """

    @overload
    def __init__(self, name: str, index: Sequence[int]) -> None:
        """
        Constructs an id with an arbitrary-dimensional index

        :param name: The readable name for the register
        :param index: The index vector
        """

    def __copy__(self) -> Bit: ...

    def __deepcopy__(self, arg: dict, /) -> Bit: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __hash__(self) -> int: ...

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...

    def to_list(self) -> list:
        """
        Return a JSON serializable list representation of the Bit.

        :return: list containing register name and index
        """

    @staticmethod
    def from_list(arg: list, /) -> Bit:
        """
        Construct Bit instance from JSON serializable list representation of the Bit.
        """

    def __and__(self: Union[pytket.circuit.logic_exp.LogicExp, Bit, int], other: Union[pytket.circuit.logic_exp.LogicExp, Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp: ...

    def __rand__(self: Union[pytket.circuit.logic_exp.LogicExp, Bit, int], other: Union[pytket.circuit.logic_exp.LogicExp, Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp: ...

    def __or__(self: Union[pytket.circuit.logic_exp.LogicExp, Bit, int], other: Union[pytket.circuit.logic_exp.LogicExp, Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp: ...

    def __ror__(self: Union[pytket.circuit.logic_exp.LogicExp, Bit, int], other: Union[pytket.circuit.logic_exp.LogicExp, Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp: ...

    def __xor__(self: Union[pytket.circuit.logic_exp.LogicExp, Bit, int], other: Union[pytket.circuit.logic_exp.LogicExp, Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp: ...

    def __rxor__(self: Union[pytket.circuit.logic_exp.LogicExp, Bit, int], other: Union[pytket.circuit.logic_exp.LogicExp, Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp: ...

class WasmState(UnitID):
    """A handle to a wasmstate"""

    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default wasm register

        :param index: The index in the register
        """

    def __copy__(self) -> WasmState: ...

    def __deepcopy__(self, arg: dict, /) -> WasmState: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __hash__(self) -> int: ...

    def __getstate__(self) -> tuple: ...

    def __setstate__(self, arg: tuple, /) -> None: ...

    def to_list(self) -> list:
        """
        Return a JSON serializable list representation of the WasmState.

        :return: list containing register name and index
        """

    @staticmethod
    def from_list(arg: list, /) -> WasmState:
        """
        Construct WasmState instance from JSON serializable list representation of the WasmState.
        """

class Node(Qubit):
    """A handle to a device node"""

    @overload
    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default physical register

        :param index: The index in the register
        """

    @overload
    def __init__(self, name: str, index: int) -> None:
        """
        Constructs an indexed id (i.e. corresponding to an element in a linear register)

        :param name: The readable name for the register
        :param index: The numerical index
        """

    @overload
    def __init__(self, name: str, row: int, col: int) -> None:
        """
        Constructs a doubly-indexed id (i.e. corresponding to an element in a grid register)

        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        """

    @overload
    def __init__(self, name: str, row: int, col: int, layer: int) -> None:
        """
        Constructs a triply-indexed id (i.e. corresponding to an element in a 3D grid register)

        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        :param layer: The layer index
        """

    @overload
    def __init__(self, name: str, index: Sequence[int]) -> None:
        """
        Constructs an id with an arbitrary-dimensional index

        :param name: The readable name for the register
        :param index: The index vector
        """

    def __copy__(self) -> Node: ...

    def __deepcopy__(self, arg: dict, /) -> Node: ...

    def to_list(self) -> list:
        """:return: a JSON serializable list representation of the Node"""

    @staticmethod
    def from_list(arg: list, /) -> Node:
        """
        Construct Node instance from JSON serializable list representation of the Node.
        """

class BitRegister:
    """Linear register of UnitID types."""

    def __init__(self, name: str, size: int) -> None:
        """
        Construct a new BitRegister.

        :param name: Name of the register.
        :param size: Size of register.
        """

    def __getitem__(self, arg: int, /) -> Bit: ...

    def __lt__(self, arg: BitRegister, /) -> bool: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __contains__(self, arg: Bit, /) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def __iter__(self) -> BitRegister: ...

    @property
    def name(self) -> str:
        """Name of register."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def size(self) -> int:
        """Size of register."""

    @size.setter
    def size(self, arg: int, /) -> None: ...

    @property
    def _current(self) -> int:
        """Internal property for iteration."""

    @_current.setter
    def _current(self, arg: int, /) -> None: ...

    def to_list(self) -> list[Bit]: ...

    def __hash__(self) -> int: ...

    def __copy__(self) -> BitRegister: ...

    def __deepcopy__(self, arg: dict, /) -> BitRegister: ...

    def __next__(self: BitRegister) -> Bit: ...

    def __and__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __rand__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __or__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __ror__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __xor__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __rxor__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __add__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __sub__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __mul__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __floordiv__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __pow__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __lshift__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

    def __rshift__(self: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int], other: Union[pytket.circuit.logic_exp.LogicExp, BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp: ...

class QubitRegister:
    """Linear register of UnitID types."""

    def __init__(self, name: str, size: int) -> None:
        """
        Construct a new QubitRegister.

        :param name: Name of the register.
        :param size: Size of register.
        """

    def __getitem__(self, arg: int, /) -> Qubit: ...

    def __lt__(self, arg: QubitRegister, /) -> bool: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __contains__(self, arg: Qubit, /) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def __iter__(self) -> QubitRegister: ...

    @property
    def name(self) -> str:
        """Name of register."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def size(self) -> int:
        """Size of register."""

    @size.setter
    def size(self, arg: int, /) -> None: ...

    @property
    def _current(self) -> int:
        """Internal property for iteration."""

    @_current.setter
    def _current(self, arg: int, /) -> None: ...

    def to_list(self) -> list[Qubit]: ...

    def __hash__(self) -> int: ...

    def __copy__(self) -> QubitRegister: ...

    def __deepcopy__(self, arg: dict, /) -> QubitRegister: ...

    def __next__(self: QubitRegister) -> Qubit: ...
