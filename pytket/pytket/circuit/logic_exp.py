# Copyright 2019-2021 Cambridge Quantum Computing
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
    Callable,
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


Ops = Union[BitWiseOp, RegWiseOp]  # all op enum types

Constant = int  # constants in expression
Variable = Union[Bit, BitRegister]  # variables in expression
ArgType = Union["LogicExp", Variable, Constant]  # all possible arguments in expression


@dataclass
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
        if not cls.op_cls_dict:
            cls.op_cls_dict = {
                cl.op: cl  # type:ignore
                for cl in (
                    BitAnd,
                    BitOr,
                    BitXor,
                    BitNot,
                    BitEq,
                    RegAnd,
                    RegOr,
                    RegXor,
                    RegEq,
                    RegNeq,
                    RegLt,
                    RegGt,
                    RegLeq,
                    RegGeq,
                )
            }
        return cls.op_cls_dict[op]

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
            rval = self._const_eval(cast(List[Constant], self.args))
        return rval

    def all_inputs(self) -> Set[Variable]:
        """
        :return: All variables involved in expression.
        :rtype: Set[Variable]
        """
        outset: Set[Variable] = set()

        var_type = Bit if isinstance(self, BitLogicExp) else BitRegister
        for arg in self.args:
            if isinstance(arg, var_type):
                outset.add(arg)
            elif isinstance(arg, LogicExp):
                outset.update(arg.all_inputs())
        return outset

    def __eq__(self, other: ArgType) -> bool:
        if not isinstance(other, LogicExp):
            return False
        return (self.op == other.op) and (self.args == other.args)

    def to_dict(self) -> Dict[str, Any]:
        """Output JSON serializable nested dictionary."""
        out: Dict[str, Any] = {"op": str(self.op)}
        args_ser: List[Union[Dict[str, Any], Constant, List[Union[str, int]]]] = []

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
        op = next(o for o in opset if o.name == op_name)
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

        return cls.factory(op)(*args)  # type: ignore


class BitLogicExp(LogicExp):
    """Expression acting only on Bit or Constant types."""

    def __init__(self, op: Ops, args: List[ArgType]):
        assert all(isinstance(a, (Bit, BitLogicExp, Constant)) for a in args)
        assert all(a in (0, 1) for a in args if isinstance(a, Constant))
        super().__init__(op, args)


class RegLogicExp(LogicExp):
    """Expression acting only on BitRegister or Constant types."""

    def __init__(self, op: Ops, args: List[ArgType]):
        assert all(isinstance(a, (BitRegister, RegLogicExp, Constant)) for a in args)
        super().__init__(op, args)


class BinaryOp(LogicExp):
    """Expresion for operation on two arguments."""

    def __init__(self, arg1: ArgType, arg2: ArgType):
        super().__init__(self.op, [arg1, arg2])

    def __str__(self) -> str:
        return f"({self.args[0]} {self.op.value} {self.args[1]})"


class UnaryOp(LogicExp):
    """Expression for operation on one argument."""

    def __init__(self, arg: ArgType):
        super().__init__(self.op, [arg])

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
    op = BitWiseOp.AND


class BitOr(Or, BitLogicExp):
    op = BitWiseOp.OR

    def eval_vals(self) -> ArgType:
        rval: ArgType = super().eval_vals()
        if 1 in self.args:
            return 1
        return rval


class BitXor(Xor, BitLogicExp):
    op = BitWiseOp.XOR


class BitNot(UnaryOp, BitLogicExp):
    op = BitWiseOp.NOT

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return 1 - args[0]


class RegAnd(And, RegLogicExp):
    op = RegWiseOp.AND


class RegOr(Or, RegLogicExp):
    op = RegWiseOp.OR


class RegXor(Xor, RegLogicExp):
    op = RegWiseOp.XOR


class ConstPredicate(BinaryOp):
    """A binary predicate where at least one of the arguments is constant."""

    def __init__(self, exp: Union[LogicExp, Variable], const: Constant) -> None:
        super().__init__(exp, const)
        assert isinstance(exp, (LogicExp, Bit, BitRegister))
        assert isinstance(const, Constant)


class Eq(ConstPredicate):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] == args[1]


class Neq(ConstPredicate):
    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return 1 - Eq._const_eval(args)


class BitEq(Eq, BitLogicExp):
    op = BitWiseOp.EQ


class RegEq(Eq, RegLogicExp):
    op = RegWiseOp.EQ


class RegNeq(Neq, RegLogicExp):
    op = RegWiseOp.NEQ


class RegLt(ConstPredicate, RegLogicExp):
    op = RegWiseOp.LT

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] < args[1]


class RegGt(ConstPredicate, RegLogicExp):
    op = RegWiseOp.GT

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] > args[1]


class RegLeq(ConstPredicate, RegLogicExp):
    op = RegWiseOp.LEQ

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] <= args[1]


class RegGeq(ConstPredicate, RegLogicExp):
    op = RegWiseOp.GEQ

    @staticmethod
    def _const_eval(args: List[Constant]) -> Constant:
        return args[0] >= args[1]


# utility to define register comparison methods
def gen_const_regpredicate(
    op: Ops,
) -> Callable[[Union[RegLogicExp, BitRegister], Constant], ConstPredicate]:
    def const_predicate(
        register: Union[RegLogicExp, BitRegister], value: Constant
    ) -> ConstPredicate:
        return cast(Type[ConstPredicate], LogicExp.factory(op))(register, value)

    return const_predicate


reg_eq = gen_const_regpredicate(RegWiseOp.EQ)
reg_eq.__doc__ = """Function to express a BitRegister equality predicate, i.e.
    for a register ``r``, ``(r == 5)`` is expressed as ``reg_eq(r, 5)``"""
reg_neq = gen_const_regpredicate(RegWiseOp.NEQ)
reg_neq.__doc__ = """Function to express a BitRegister inequality predicate, i.e.
    for a register ``r``, ``(r != 5)`` is expressed as ``reg_neq(r, 5)``"""
reg_lt = gen_const_regpredicate(RegWiseOp.LT)
reg_lt.__doc__ = """Function to express a BitRegister less than predicate, i.e.
    for a register ``r``, ``(r < 5)`` is expressed as ``reg_lt(r, 5)``"""
reg_gt = gen_const_regpredicate(RegWiseOp.GT)
reg_gt.__doc__ = """Function to express a BitRegister greater than predicate, i.e.
    for a register ``r``, ``(r > 5)`` is expressed as ``reg_gt(r, 5)``"""
reg_leq = gen_const_regpredicate(RegWiseOp.LEQ)
reg_leq.__doc__ = """Function to express a BitRegister less than or equal to predicate,
    i.e. for a register ``r``, ``(r <= 5)`` is expressed as ``reg_leq(r, 5)``"""
reg_geq = gen_const_regpredicate(RegWiseOp.GEQ)
reg_geq.__doc__ = """Function to express a BitRegister greater than or equal to
    predicate, i.e. for a register ``r``, ``(r >= 5)`` is expressed as
    ``reg_geq(r, 5)``"""


def if_bit(bit: Union[Bit, BitLogicExp]) -> ConstPredicate:
    """Equivalent of ``if bit:``."""
    return BitEq(bit, 1)


def if_not_bit(bit: Union[Bit, BitLogicExp]) -> ConstPredicate:
    """Equivalent of ``if not bit:``."""
    return BitEq(bit, 0)
