from typing import ClassVar, List

from typing import overload
_DEBUG_ONE_REG_PREFIX: str
_DEBUG_ZERO_REG_PREFIX: str
_TEMP_BIT_NAME: str
_TEMP_BIT_REG_BASE: str
_TEMP_REG_SIZE: int

class Bit(UnitID):
    __and__: ClassVar[function] = ...
    __or__: ClassVar[function] = ...
    __rand__: ClassVar[function] = ...
    __ror__: ClassVar[function] = ...
    __rxor__: ClassVar[function] = ...
    __xor__: ClassVar[function] = ...
    @overload
    def __init__(self, index: int) -> None: ...
    @overload
    def __init__(self, name: str) -> None: ...
    @overload
    def __init__(self, name: str, index: int) -> None: ...
    @overload
    def __init__(self, name: str, row: int, col: int) -> None: ...
    @overload
    def __init__(self, name: str, index: List[int]) -> None: ...
    @classmethod
    def from_list(cls, arg0: list) -> Bit: ...
    def to_list(self) -> list: ...
    def __eq__(self, arg0: object) -> bool: ...
    def __hash__(self) -> int: ...

class BitRegister:
    __add__: ClassVar[function] = ...
    __and__: ClassVar[function] = ...
    __floordiv__: ClassVar[function] = ...
    __lshift__: ClassVar[function] = ...
    __mul__: ClassVar[function] = ...
    __or__: ClassVar[function] = ...
    __pow__: ClassVar[function] = ...
    __rand__: ClassVar[function] = ...
    __ror__: ClassVar[function] = ...
    __rshift__: ClassVar[function] = ...
    __rxor__: ClassVar[function] = ...
    __sub__: ClassVar[function] = ...
    __xor__: ClassVar[function] = ...
    name: str
    size: int
    def __init__(self, name: str, size: int) -> None: ...
    def __contains__(self, arg0: Bit) -> bool: ...
    def __copy__(self) -> BitRegister: ...
    def __deepcopy__(self, arg0: dict) -> BitRegister: ...
    def __eq__(self, arg0: object) -> bool: ...
    def __getitem__(self, arg0: int) -> Bit: ...
    def __hash__(self) -> int: ...
    def __len__(self) -> int: ...
    def __lt__(self, arg0: BitRegister) -> bool: ...

class Node(Qubit):
    @overload
    def __init__(self, index: int) -> None: ...
    @overload
    def __init__(self, name: str, index: int) -> None: ...
    @overload
    def __init__(self, name: str, row: int, col: int) -> None: ...
    @overload
    def __init__(self, name: str, row: int, col: int, layer: int) -> None: ...
    @overload
    def __init__(self, name: str, index: List[int]) -> None: ...
    @classmethod
    def from_list(cls, arg0: list) -> Node: ...
    def to_list(self) -> list: ...

class Qubit(UnitID):
    @overload
    def __init__(self, index: int) -> None: ...
    @overload
    def __init__(self, name: str) -> None: ...
    @overload
    def __init__(self, name: str, index: int) -> None: ...
    @overload
    def __init__(self, name: str, row: int, col: int) -> None: ...
    @overload
    def __init__(self, name: str, index: List[int]) -> None: ...
    @classmethod
    def from_list(cls, arg0: list) -> Qubit: ...
    def to_list(self) -> list: ...
    def __getstate__(self) -> tuple: ...
    def __setstate__(self, arg0: tuple) -> None: ...

class QubitRegister:
    name: str
    size: int
    def __init__(self, name: str, size: int) -> None: ...
    def __contains__(self, arg0: Qubit) -> bool: ...
    def __copy__(self) -> QubitRegister: ...
    def __deepcopy__(self, arg0: dict) -> QubitRegister: ...
    def __eq__(self, arg0: object) -> bool: ...
    def __getitem__(self, arg0: int) -> Qubit: ...
    def __hash__(self) -> int: ...
    def __len__(self) -> int: ...
    def __lt__(self, arg0: QubitRegister) -> bool: ...

class UnitID:
    def __init__(self) -> None: ...
    def __copy__(self) -> UnitID: ...
    def __deepcopy__(self, arg0: dict) -> UnitID: ...
    def __eq__(self, arg0: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __lt__(self, arg0: UnitID) -> bool: ...
    @property
    def index(self) -> List[int]: ...
    @property
    def reg_name(self) -> str: ...
    @property
    def type(self) -> UnitType: ...

class UnitType:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    bit: ClassVar[UnitType] = ...
    qubit: ClassVar[UnitType] = ...
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
