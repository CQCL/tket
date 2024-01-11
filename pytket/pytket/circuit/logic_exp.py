# Copyright 2019-2024 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Classes and functions for constructing logical
 expressions over Bit and BitRegister."""
from typing import (
    Any,
    Iterable,
    Set,
    Tuple,
    Type,
    List,
    Union,
    Dict,
    ClassVar,
    Iterator,
    TypeVar,
    cast,
    Sequence,
)
from enum import Enum
from dataclasses import dataclass

from pytket.circuit import Bit, BitRegister

T = TypeVar("T")


def filter_by_type(seq: Iterable, var_type: Type[T]) -> Iterator[Tuple[int, T]]:
    """Return enumeration of seq, with only elements of type var_type."""
    return filter(lambda x: isinstance(x[1], var_type), enumerate(seq))


class BitWiseOp(Enum):
    """Enum for operations on Bit."""

    AND = "&"
    OR = "|"
    XOR = "^"
    EQ = "=="
    NEQ = "!="
    NOT = "~"


class RegWiseOp(Enum):
    """Enum for operations on BitRegister."""

    AND = "&"
    OR = "|"
    XOR = "^"
    EQ = "=="
    NEQ = "!="
    LT = "<"
    GT = ">"
    LEQ = "<="
    GEQ = ">="
    ADD = "+"
    SUB = "-"
    MUL = "*"
    DIV = "/"
    POW = "**"
    LSH = "<<"
    RSH = ">>"
    NOT = "~"
    NEG = "-"


Ops = Union[BitWiseOp, RegWiseOp]  # all op enum types


Constant = int  # constants in expression
Variable = Union[Bit, BitRegister]  # variables in expression
ArgType = Union["LogicExp", Variable, Constant]  # all possible arguments in expression


@dataclass(init=False)
class LogicExp:
    """Logical expressions over Bit or BitRegister.
    Encoded as a tree of expressions"""

    op: Ops  # enum for operation encoded by this node
    args: List[ArgType]  # arguments of operation
    # class level dictionary mapping enum to class
    op_cls_dict: ClassVar[Dict[Ops, Type["LogicExp"]]] = dict()

    @classmethod
    def factory(cls, op: Ops) -> Type["LogicExp"]:
        """Return matching operation class for enum."""
        # RegNeg cannot be initialised this way as "-" clashes with SUB
        if op == BitWiseOp.AND:
            return BitAnd
        if op == BitWiseOp.OR:
            return BitOr
        if op == BitWiseOp.XOR:
            return BitXor
        if op == BitWiseOp.NOT:
            return BitNot
        if op == BitWiseOp.EQ:
            return BitEq
        if op == BitWiseOp.NEQ:
            return BitNeq
        if op == RegWiseOp.AND:
            return RegAnd
        if op == RegWiseOp.OR:
            return RegOr
        if op == RegWiseOp.XOR:
            return RegXor
        if op == RegWiseOp.ADD:
            return RegAdd
        if op == RegWiseOp.SUB:
            return RegSub
        if op == RegWiseOp.MUL:
            return RegMul
        if op == RegWiseOp.DIV:
            return RegDiv
        if op == RegWiseOp.POW:
            return RegPow
        if op == RegWiseOp.LSH:
            return RegLsh
        if op == RegWiseOp.RSH:
            return RegRsh
        if op == RegWiseOp.EQ:
            return RegEq
        if op == RegWiseOp.NEQ:
            return RegNeq
        if op == RegWiseOp.LT:
            return RegLt
        if op == RegWiseOp.GT:
            return RegGt
        if op == RegWiseOp.LEQ:
            return RegLeq
        if op == RegWiseOp.GEQ:
            return RegGeq
        if op == RegWiseOp.NOT:
            return RegNot
        raise ValueError("op type not supported")

    def set_value(self, var: Variable, val: Constant) -> None:
        """Set value of var to val recursively."""
        for i, arg in enumerate(self.args):
            if isinstance(arg, (Bit, BitRegister)):
                if arg == var:
                    self.args[i] = val
            elif isinstance(arg, LogicExp):
                arg.set_value(var, val)

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        """Evaluate expression given constant values for all args."""
        raise NotImplementedError

    def eval_vals(self) -> ArgType:
        """Attempt to evaluate all sub-expressions; simple constant folding."""
        rval: ArgType = self
        for i, arg in filter_by_type(self.args, LogicExp):
            self.args[i] = arg.eval_vals()
        if all(isinstance(a, Constant) for a in self.args):
            try:
                rval = self._const_eval(cast(List[Constant], self.args))
            except NotImplementedError:
                pass
        return rval

    def all_inputs(self) -> Set[Variable]:
        """
        :return: All variables involved in expression.
        :rtype: Set[Variable]
        """
        outset: Set[Variable] = set()

        for arg in self.args:
            if isinstance(arg, LogicExp):
                outset.update(arg.all_inputs())
                continue
            if isinstance(self, BitLogicExp):
                if isinstance(arg, Bit):
                    outset.add(arg)
            else:
                if isinstance(arg, BitRegister):
                    outset.add(arg)
        return outset

    def all_inputs_ordered(self) -> list[Variable]:
        """
        :return: All variables involved in expression, in order of first appearance.
        :rtype: list[Variable]
        """
        # use dict[Variable, None] instead of set[Variable] to preserve order
        outset: dict[Variable, None] = {}

        for arg in self.args:
            if isinstance(arg, LogicExp):
                outset.update(dict.fromkeys(arg.all_inputs_ordered()))
                continue
            if isinstance(self, BitLogicExp):
                if isinstance(arg, Bit):
                    outset[arg] = None
            else:
                if isinstance(arg, BitRegister):
                    outset[arg] = None
        return list(outset)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LogicExp):
            return False
        return (self.op == other.op) and (self.args == other.args)

    def to_dict(self) -> Dict[str, Any]:
        """Output JSON serializable nested dictionary."""
        out: Dict[str, Any] = {"op": str(self.op)}
        args_ser: List[Union[Dict, Constant, List[Union[str, int]]]] = []

        for arg in self.args:
            if isinstance(arg, LogicExp):
                args_ser.append(arg.to_dict())
            elif isinstance(arg, Constant):
                args_ser.append(arg)
            elif isinstance(arg, Bit):
                args_ser.append(arg.to_list())
            elif isinstance(arg, BitRegister):
                args_ser.append({"name": arg.name, "size": arg.size})

        out["args"] = args_ser
        return out

    @classmethod
    def from_dict(cls, dic: Dict[str, Any]) -> "LogicExp":
        """Load from JSON serializable nested dictionary."""
        opset_name, op_name = dic["op"].split(".", 2)
        opset = BitWiseOp if opset_name == "BitWiseOp" else RegWiseOp
        op = cast(
            Union[BitWiseOp, RegWiseOp], next(o for o in opset if o.name == op_name)
        )
        args: List[ArgType] = []
        for arg_ser in dic["args"]:
            if isinstance(arg_ser, Constant):
                args.append(arg_ser)
            elif isinstance(arg_ser, List):
                args.append(Bit(arg_ser[0], arg_ser[1]))
            elif isinstance(arg_ser, Dict):
                if "op" in arg_ser:
                    args.append(LogicExp.from_dict(arg_ser))
                else:
                    args.append(BitRegister(arg_ser["name"], arg_ser["size"]))
        return create_logic_exp(op, args)

    def _rename_args_recursive(
        self, cmap: Dict[Bit, Bit], renamed_regs: Set[str]
    ) -> bool:
        success = False
        for i, arg in enumerate(self.args):
            if isinstance(arg, Bit):
                if arg in cmap:
                    self.args[i] = cmap[arg]
                    success = True
            elif isinstance(arg, BitRegister):
                if arg.name in renamed_regs:
                    raise ValueError(
                        f"""Can't rename bits in {arg.__repr__()} """
                        """because the register is being used """
                        """in a register-wise logic expression."""
                    )
            elif isinstance(arg, LogicExp):
                success |= arg._rename_args_recursive(cmap, renamed_regs)
        return success

    def rename_args(self, cmap: Dict[Bit, Bit]) -> bool:
        """Rename the Bits according to a Bit map. Raise ValueError if
        a bit is being used in a register-wise expression.
        """
        renamed_regs = set([key.reg_name for key in cmap.keys()])
        return self._rename_args_recursive(cmap, renamed_regs)


BitArgType = Union[LogicExp, Bit, Constant]
RegArgType = Union[LogicExp, BitRegister, Constant]


class BitLogicExp(LogicExp):
    """Expression acting only on Bit or Constant types."""

    def __and__(self, other: BitArgType) -> "BitAnd":
        return BitAnd(self, other)

    def __rand__(self, other: BitArgType) -> "BitAnd":
        return BitAnd(self, other)

    def __or__(self, other: BitArgType) -> "BitOr":
        return BitOr(self, other)

    def __ror__(self, other: BitArgType) -> "BitOr":
        return BitOr(self, other)

    def __xor__(self, other: BitArgType) -> "BitXor":
        return BitXor(self, other)

    def __rxor__(self, other: BitArgType) -> "BitXor":
        return BitXor(self, other)


class RegLogicExp(LogicExp):
    """Expression acting only on BitRegister or Constant types."""

    def __and__(self, other: RegArgType) -> "RegAnd":
        return RegAnd(self, other)

    def __rand__(self, other: RegArgType) -> "RegAnd":
        return RegAnd(self, other)

    def __or__(self, other: RegArgType) -> "RegOr":
        return RegOr(self, other)

    def __ror__(self, other: RegArgType) -> "RegOr":
        return RegOr(self, other)

    def __xor__(self, other: RegArgType) -> "RegXor":
        return RegXor(self, other)

    def __rxor__(self, other: RegArgType) -> "RegXor":
        return RegXor(self, other)

    def __add__(self, other: RegArgType) -> "RegAdd":
        return RegAdd(self, other)

    def __sub__(self, other: RegArgType) -> "RegSub":
        return RegSub(self, other)

    def __mul__(self, other: RegArgType) -> "RegMul":
        return RegMul(self, other)

    def __floordiv__(self, other: RegArgType) -> "RegDiv":
        return RegDiv(self, other)

    def __pow__(self, other: RegArgType) -> "RegPow":
        return RegPow(self, other)

    def __lshift__(self, other: RegArgType) -> "RegLsh":
        return RegLsh(self, other)

    def __rshift__(self, other: RegArgType) -> "RegRsh":
        return RegRsh(self, other)


class BinaryOp(LogicExp):
    """Expresion for operation on two arguments."""

    def __str__(self) -> str:
        return f"({self.args[0]} {self.op.value} {self.args[1]})"


class UnaryOp(LogicExp):
    """Expression for operation on one argument."""

    def __str__(self) -> str:
        return f"({self.op.value} {self.args[0]})"


class And(BinaryOp):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] & args[1]

    def eval_vals(self) -> ArgType:
        rval: ArgType = super().eval_vals()
        if 0 in self.args:
            return 0
        return rval


class Or(BinaryOp):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] | args[1]


class Xor(BinaryOp):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] ^ args[1]


class BitAnd(And, BitLogicExp):
    def __init__(self, arg1: BitArgType, arg2: BitArgType) -> None:
        self.op = BitWiseOp.AND
        self.args = [arg1, arg2]


class BitOr(Or, BitLogicExp):
    def __init__(self, arg1: BitArgType, arg2: BitArgType) -> None:
        self.op = BitWiseOp.OR
        self.args = [arg1, arg2]

    def eval_vals(self) -> ArgType:
        rval: ArgType = super().eval_vals()
        if 1 in self.args:
            return 1
        return rval


class BitXor(Xor, BitLogicExp):
    def __init__(self, arg1: BitArgType, arg2: BitArgType) -> None:
        self.op = BitWiseOp.XOR
        self.args = [arg1, arg2]


class BitNot(UnaryOp, BitLogicExp):
    def __init__(self, arg1: BitArgType) -> None:
        self.op = BitWiseOp.NOT
        self.args = [arg1]

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return 1 - args[0]


class RegAnd(And, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.AND
        self.args = [arg1, arg2]


class RegOr(Or, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.OR
        self.args = [arg1, arg2]


class RegXor(Xor, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.XOR
        self.args = [arg1, arg2]


class RegAdd(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.ADD
        self.args = [arg1, arg2]


class RegSub(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.SUB
        self.args = [arg1, arg2]


class RegMul(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.MUL
        self.args = [arg1, arg2]


class RegDiv(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.DIV
        self.args = [arg1, arg2]


class RegPow(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.POW
        self.args = [arg1, arg2]


class RegLsh(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.LSH
        self.args = [arg1, arg2]


class RegNeg(UnaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType) -> None:
        self.op = RegWiseOp.NEG
        self.args = [arg1]


class RegNot(UnaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType) -> None:
        self.op = RegWiseOp.NOT
        self.args = [arg1]


class RegRsh(BinaryOp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.RSH
        self.args = [arg1, arg2]


class PredicateExp(BinaryOp):
    """
    A binary predicate where the arguments are either
    Bits, BitRegisters, or Constants.
    """


class Eq(PredicateExp):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] == args[1]


class Neq(PredicateExp):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return 1 - Eq._const_eval(args)


class BitEq(Eq, BitLogicExp):
    def __init__(self, arg1: BitArgType, arg2: BitArgType) -> None:
        self.op = BitWiseOp.EQ
        self.args = [arg1, arg2]


class BitNeq(Neq, BitLogicExp):
    def __init__(self, arg1: BitArgType, arg2: BitArgType) -> None:
        self.op = BitWiseOp.NEQ
        self.args = [arg1, arg2]


class RegEq(Eq, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.EQ
        self.args = [arg1, arg2]


class RegNeq(Neq, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.NEQ
        self.args = [arg1, arg2]


class RegLt(PredicateExp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.LT
        self.args = [arg1, arg2]

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] < args[1]


class RegGt(PredicateExp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.GT
        self.args = [arg1, arg2]

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] > args[1]


class RegLeq(PredicateExp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.LEQ
        self.args = [arg1, arg2]

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] <= args[1]


class RegGeq(PredicateExp, RegLogicExp):
    def __init__(self, arg1: RegArgType, arg2: RegArgType) -> None:
        self.op = RegWiseOp.GEQ
        self.args = [arg1, arg2]

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] >= args[1]


def reg_eq(register: Union[RegLogicExp, BitRegister], value: Constant) -> RegLogicExp:
    """Function to express a BitRegister equality predicate, i.e.
    for a register ``r``, ``(r == 5)`` is expressed as ``reg_eq(r, 5)``"""
    return RegEq(register, value)


def reg_neq(register: Union[RegLogicExp, BitRegister], value: Constant) -> RegLogicExp:
    """Function to express a BitRegister inequality predicate, i.e.
    for a register ``r``, ``(r != 5)`` is expressed as ``reg_neq(r, 5)``"""
    return RegNeq(register, value)


def reg_lt(register: Union[RegLogicExp, BitRegister], value: Constant) -> RegLogicExp:
    """Function to express a BitRegister less than predicate, i.e.
    for a register ``r``, ``(r < 5)`` is expressed as ``reg_lt(r, 5)``"""
    return RegLt(register, value)


def reg_gt(register: Union[RegLogicExp, BitRegister], value: Constant) -> RegLogicExp:
    """Function to express a BitRegister greater than predicate, i.e.
    for a register ``r``, ``(r > 5)`` is expressed as ``reg_gt(r, 5)``"""
    return RegGt(register, value)


def reg_leq(register: Union[RegLogicExp, BitRegister], value: Constant) -> RegLogicExp:
    """Function to express a BitRegister less than or equal to predicate,
    i.e. for a register ``r``, ``(r <= 5)`` is expressed as ``reg_leq(r, 5)``"""
    return RegLeq(register, value)


def reg_geq(register: Union[RegLogicExp, BitRegister], value: Constant) -> RegLogicExp:
    """Function to express a BitRegister greater than or equal to
    predicate, i.e. for a register ``r``, ``(r >= 5)`` is expressed as
    ``reg_geq(r, 5)``"""
    return RegGeq(register, value)


def if_bit(bit: Union[Bit, BitLogicExp]) -> PredicateExp:
    """Equivalent of ``if bit:``."""
    return BitEq(bit, 1)


def if_not_bit(bit: Union[Bit, BitLogicExp]) -> PredicateExp:
    """Equivalent of ``if not bit:``."""
    return BitEq(bit, 0)


def create_bit_logic_exp(op: BitWiseOp, args: Sequence[BitArgType]) -> BitLogicExp:
    if op == BitWiseOp.AND:
        assert len(args) == 2
        return BitAnd(args[0], args[1])
    if op == BitWiseOp.OR:
        assert len(args) == 2
        return BitOr(args[0], args[1])
    if op == BitWiseOp.XOR:
        assert len(args) == 2
        return BitXor(args[0], args[1])
    if op == BitWiseOp.NOT:
        assert len(args) == 1
        return BitNot(args[0])
    if op == BitWiseOp.EQ:
        assert len(args) == 2
        return BitEq(args[0], args[1])
    if op == BitWiseOp.NEQ:
        assert len(args) == 2
        return BitNeq(args[0], args[1])


def create_reg_logic_exp(op: RegWiseOp, args: Sequence[RegArgType]) -> RegLogicExp:
    if op == RegWiseOp.AND:
        assert len(args) == 2
        return RegAnd(args[0], args[1])
    if op == RegWiseOp.OR:
        assert len(args) == 2
        return RegOr(args[0], args[1])
    if op == RegWiseOp.XOR:
        assert len(args) == 2
        return RegXor(args[0], args[1])
    if op == RegWiseOp.ADD:
        assert len(args) == 2
        return RegAdd(args[0], args[1])
    if op == RegWiseOp.SUB:
        if len(args) == 2:
            return RegSub(args[0], args[1])
        if len(args) == 1:
            return RegNeg(args[0])
    if op == RegWiseOp.NEG:
        assert len(args) == 1
        return RegNeg(args[0])
    if op == RegWiseOp.MUL:
        assert len(args) == 2
        return RegMul(args[0], args[1])
    if op == RegWiseOp.DIV:
        assert len(args) == 2
        return RegDiv(args[0], args[1])
    if op == RegWiseOp.POW:
        assert len(args) == 2
        return RegPow(args[0], args[1])
    if op == RegWiseOp.LSH:
        assert len(args) == 2
        return RegLsh(args[0], args[1])
    if op == RegWiseOp.RSH:
        assert len(args) == 2
        return RegRsh(args[0], args[1])
    if op == RegWiseOp.EQ:
        assert len(args) == 2
        return RegEq(args[0], args[1])
    if op == RegWiseOp.NEQ:
        assert len(args) == 2
        return RegNeq(args[0], args[1])
    if op == RegWiseOp.LT:
        assert len(args) == 2
        return RegLt(args[0], args[1])
    if op == RegWiseOp.GT:
        assert len(args) == 2
        return RegGt(args[0], args[1])
    if op == RegWiseOp.LEQ:
        assert len(args) == 2
        return RegLeq(args[0], args[1])
    if op == RegWiseOp.GEQ:
        assert len(args) == 2
        return RegGeq(args[0], args[1])
    if op == RegWiseOp.NOT:
        assert len(args) == 1
        return RegNot(args[0])
    raise ValueError("op type not supported")


def create_logic_exp(op: Ops, args: Sequence[ArgType]) -> LogicExp:
    if isinstance(op, BitWiseOp):
        bit_args = []
        for arg in args:
            assert isinstance(arg, (BitLogicExp, Bit, Constant))
            bit_args.append(arg)
        return create_bit_logic_exp(op, bit_args)
    else:
        assert isinstance(op, RegWiseOp)
        reg_args = []
        for arg in args:
            assert isinstance(arg, (RegLogicExp, BitRegister, Constant))
            reg_args.append(arg)
        return create_reg_logic_exp(op, reg_args)


def create_predicate_exp(op: Ops, args: Sequence[ArgType]) -> PredicateExp:
    if op == BitWiseOp.EQ:
        assert len(args) == 2
        assert isinstance(args[0], (BitLogicExp, Bit, int))
        assert isinstance(args[1], (BitLogicExp, Bit, int))
        return BitEq(args[0], args[1])
    if op == BitWiseOp.NEQ:
        assert len(args) == 2
        assert isinstance(args[0], (BitLogicExp, Bit, int))
        assert isinstance(args[1], (BitLogicExp, Bit, int))
        return BitNeq(args[0], args[1])
    if op == RegWiseOp.EQ:
        assert len(args) == 2
        assert isinstance(args[0], (RegLogicExp, BitRegister, int))
        assert isinstance(args[1], (RegLogicExp, BitRegister, int))
        return RegEq(args[0], args[1])
    if op == RegWiseOp.NEQ:
        assert len(args) == 2
        assert isinstance(args[0], (RegLogicExp, BitRegister, int))
        assert isinstance(args[1], (RegLogicExp, BitRegister, int))
        return RegNeq(args[0], args[1])
    if op == RegWiseOp.LT:
        assert len(args) == 2
        assert isinstance(args[0], (RegLogicExp, BitRegister, int))
        assert isinstance(args[1], (RegLogicExp, BitRegister, int))
        return RegLt(args[0], args[1])
    if op == RegWiseOp.GT:
        assert len(args) == 2
        assert isinstance(args[0], (RegLogicExp, BitRegister, int))
        assert isinstance(args[1], (RegLogicExp, BitRegister, int))
        return RegGt(args[0], args[1])
    if op == RegWiseOp.LEQ:
        assert len(args) == 2
        assert isinstance(args[0], (RegLogicExp, BitRegister, int))
        assert isinstance(args[1], (RegLogicExp, BitRegister, int))
        return RegLeq(args[0], args[1])
    if op == RegWiseOp.GEQ:
        assert len(args) == 2
        assert isinstance(args[0], (RegLogicExp, BitRegister, int))
        assert isinstance(args[1], (RegLogicExp, BitRegister, int))
        return RegGeq(args[0], args[1])
    raise ValueError("op type not supported")


#  ConstPredicate is deprecated in favour of PredicateExp
ConstPredicate = PredicateExp
