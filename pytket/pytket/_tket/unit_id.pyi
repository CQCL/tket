from __future__ import annotations
import pytket.circuit.logic_exp
import typing
__all__ = ['Bit', 'BitRegister', 'Node', 'Qubit', 'QubitRegister', 'UnitID', 'UnitType']
class Bit(UnitID):
    """
    A handle to a bit
    """
    @staticmethod
    def from_list(arg0: list) -> Bit:
        """
        Construct Bit instance from JSON serializable list representation of the Bit.
        """
    def __and__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp:
        ...
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __hash__(self) -> int:
        ...
    @typing.overload
    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default classical register
        
        :param index: The index in the register
        """
    @typing.overload
    def __init__(self, name: str) -> None:
        """
        Constructs a named id (i.e. corresponding to a singleton register)
        
        :param name: The readable name for the id
        """
    @typing.overload
    def __init__(self, name: str, index: int) -> None:
        """
        Constructs an indexed id (i.e. corresponding to an element in a linear register)
        
        :param name: The readable name for the register
        :param index: The numerical index
        """
    @typing.overload
    def __init__(self, name: str, row: int, col: int) -> None:
        """
        Constructs a doubly-indexed id (i.e. corresponding to an element in a grid register)
        
        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        """
    @typing.overload
    def __init__(self, name: str, index: typing.Sequence[int]) -> None:
        """
        Constructs an id with an arbitrary-dimensional index
        
        :param name: The readable name for the register
        :param index: The index vector
        """
    def __or__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp:
        ...
    def __rand__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp:
        ...
    def __ror__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp:
        ...
    def __rxor__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp:
        ...
    def __xor__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.Bit, int]) -> pytket.circuit.logic_exp.BitLogicExp:
        ...
    def to_list(self) -> list:
        """
        Return a JSON serializable list representation of the Bit.
        
        :return: list containing register name and index
        """
class BitRegister:
    """
    Linear register of UnitID types.
    """
    def __add__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __and__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __contains__(self, arg0: Bit) -> bool:
        ...
    def __copy__(self) -> BitRegister:
        ...
    def __deepcopy__(self, arg0: dict) -> BitRegister:
        ...
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __floordiv__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __getitem__(self, arg0: int) -> Bit:
        ...
    def __hash__(self) -> int:
        ...
    def __init__(self, name: str, size: int) -> None:
        """
        Construct a new BitRegister.
        
        :param name: Name of the register.
        :param size: Size of register.
        """
    def __iter__(self) -> BitRegister:
        ...
    def __len__(self) -> int:
        ...
    def __lshift__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __lt__(self, arg0: BitRegister) -> bool:
        ...
    def __mul__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __next__(self) -> Bit:
        ...
    def __or__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __pow__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __rand__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __repr__(self) -> str:
        ...
    def __ror__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __rshift__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __rxor__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __str__(self) -> str:
        ...
    def __sub__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def __xor__(self, other: typing.Union[pytket.circuit.logic_exp.LogicExp, pytket._tket.unit_id.BitRegister, int]) -> pytket.circuit.logic_exp.RegLogicExp:
        ...
    def to_list(self) -> list[Bit]:
        ...
    @property
    def _current(self) -> int:
        """
        Internal property for iteration.
        """
    @_current.setter
    def _current(self, arg1: int) -> None:
        ...
    @property
    def name(self) -> str:
        """
        Name of register.
        """
    @name.setter
    def name(self, arg1: str) -> None:
        ...
    @property
    def size(self) -> int:
        """
        Size of register.
        """
    @size.setter
    def size(self, arg1: int) -> None:
        ...
class Node(Qubit):
    """
    A handle to a device node
    """
    @staticmethod
    def from_list(arg0: list) -> Node:
        """
        Construct Node instance from JSON serializable list representation of the Node.
        """
    @typing.overload
    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default physical register
        
        :param index: The index in the register
        """
    @typing.overload
    def __init__(self, name: str, index: int) -> None:
        """
        Constructs an indexed id (i.e. corresponding to an element in a linear register)
        
        :param name: The readable name for the register
        :param index: The numerical index
        """
    @typing.overload
    def __init__(self, name: str, row: int, col: int) -> None:
        """
        Constructs a doubly-indexed id (i.e. corresponding to an element in a grid register)
        
        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        """
    @typing.overload
    def __init__(self, name: str, row: int, col: int, layer: int) -> None:
        """
        Constructs a triply-indexed id (i.e. corresponding to an element in a 3D grid register)
        
        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        :param layer: The layer index
        """
    @typing.overload
    def __init__(self, name: str, index: typing.Sequence[int]) -> None:
        """
        Constructs an id with an arbitrary-dimensional index
        
        :param name: The readable name for the register
        :param index: The index vector
        """
    def to_list(self) -> list:
        """
        :return: a JSON serializable list representation of the Node
        """
class Qubit(UnitID):
    """
    A handle to a qubit
    """
    @staticmethod
    def from_list(arg0: list) -> Qubit:
        """
        Construct Qubit instance from JSON serializable list representation of the Qubit.
        """
    def __getstate__(self) -> tuple:
        ...
    @typing.overload
    def __init__(self, index: int) -> None:
        """
        Constructs an id for some index in the default qubit register
        
        :param index: The index in the register
        """
    @typing.overload
    def __init__(self, name: str) -> None:
        """
        Constructs a named id (i.e. corresponding to a singleton register)
        
        :param name: The readable name for the id
        """
    @typing.overload
    def __init__(self, name: str, index: int) -> None:
        """
        Constructs an indexed id (i.e. corresponding to an element in a linear register)
        
        :param name: The readable name for the register
        :param index: The numerical index
        """
    @typing.overload
    def __init__(self, name: str, row: int, col: int) -> None:
        """
        Constructs a doubly-indexed id (i.e. corresponding to an element in a grid register)
        
        :param name: The readable name for the register
        :param row: The row index
        :param col: The column index
        """
    @typing.overload
    def __init__(self, name: str, index: typing.Sequence[int]) -> None:
        """
        Constructs an id with an arbitrary-dimensional index
        
        :param name: The readable name for the register
        :param index: The index vector
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def to_list(self) -> list:
        """
        :return: a JSON serializable list representation of the Qubit
        """
class QubitRegister:
    """
    Linear register of UnitID types.
    """
    def __contains__(self, arg0: Qubit) -> bool:
        ...
    def __copy__(self) -> QubitRegister:
        ...
    def __deepcopy__(self, arg0: dict) -> QubitRegister:
        ...
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __getitem__(self, arg0: int) -> Qubit:
        ...
    def __hash__(self) -> int:
        ...
    def __init__(self, name: str, size: int) -> None:
        """
        Construct a new QubitRegister.
        
        :param name: Name of the register.
        :param size: Size of register.
        """
    def __iter__(self) -> QubitRegister:
        ...
    def __len__(self) -> int:
        ...
    def __lt__(self, arg0: QubitRegister) -> bool:
        ...
    def __next__(self) -> Qubit:
        ...
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    def to_list(self) -> list[Qubit]:
        ...
    @property
    def _current(self) -> int:
        """
        Internal property for iteration.
        """
    @_current.setter
    def _current(self, arg1: int) -> None:
        ...
    @property
    def name(self) -> str:
        """
        Name of register.
        """
    @name.setter
    def name(self, arg1: str) -> None:
        ...
    @property
    def size(self) -> int:
        """
        Size of register.
        """
    @size.setter
    def size(self, arg1: int) -> None:
        ...
class UnitID:
    """
    A handle to a computational unit (e.g. qubit, bit)
    """
    def __copy__(self) -> UnitID:
        ...
    def __deepcopy__(self, arg0: dict) -> UnitID:
        ...
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __lt__(self, arg0: UnitID) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def index(self) -> list[int]:
        """
        Index vector describing position in the register. The length of this vector is the dimension of the register
        """
    @property
    def reg_name(self) -> str:
        """
        Readable name of register
        """
    @property
    def type(self) -> UnitType:
        """
        Type of unit, either ``UnitType.qubit`` or ``UnitType.bit``
        """
class UnitType:
    """
    Enum for data types of units in circuits (e.g. Qubits vs Bits).
    
    Members:
    
      qubit : A single Qubit
    
      bit : A single classical Bit
    """
    __members__: typing.ClassVar[dict[str, UnitType]]  # value = {'qubit': <UnitType.qubit: 0>, 'bit': <UnitType.bit: 1>}
    bit: typing.ClassVar[UnitType]  # value = <UnitType.bit: 1>
    qubit: typing.ClassVar[UnitType]  # value = <UnitType.qubit: 0>
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
_DEBUG_ONE_REG_PREFIX: str = 'tk_DEBUG_ONE_REG'
_DEBUG_ZERO_REG_PREFIX: str = 'tk_DEBUG_ZERO_REG'
_TEMP_BIT_NAME: str = 'tk_SCRATCH_BIT'
_TEMP_BIT_REG_BASE: str = 'tk_SCRATCH_BITREG'
_TEMP_REG_SIZE: int = 32
