# Copyright 2019-2022 Cambridge Quantum Computing
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

"""The circuit module provides an API to interact with the
 tket :py:class:`Circuit` data structure.
  This module is provided in binary form during the PyPI installation."""
from typing import Any, Tuple, Type, Union, cast, Callable

from pytket._tket.circuit import *  # type: ignore
from pytket._tket.circuit import Bit, BitRegister

# prefixes for assertion bits
from pytket._tket.circuit import _DEBUG_ZERO_REG_PREFIX, _DEBUG_ONE_REG_PREFIX  # type: ignore
from pytket._tket.pauli import Pauli  # type: ignore

from .logic_exp import (
    BitLogicExp,
    BitWiseOp,
    RegLogicExp,
    RegWiseOp,
    ArgType,
    LogicExp,
    BinaryOp,
    Ops,
    if_bit,
    if_not_bit,
    reg_eq,
    reg_geq,
    reg_gt,
    reg_leq,
    reg_lt,
    reg_neq,
)

# Add ability to compare Bit equality with arbitrary class
Bit.oldeq = Bit.__eq__


def overload_biteq(self: Bit, other: Any) -> bool:
    if not isinstance(other, Bit):
        return False
    return cast(bool, self.oldeq(other))


setattr(Bit, "__eq__", overload_biteq)

# overload operators for Bit, BitRegister and expressions over these
# such that the operation returns a LogicExp describing the operation
cls_enum_pairs: Tuple[Tuple[Type, Union[Type[BitWiseOp], Type[RegWiseOp]]], ...] = (
    (Bit, BitWiseOp),
    (BitLogicExp, BitWiseOp),
    (BitRegister, RegWiseOp),
    (RegLogicExp, RegWiseOp),
)


def gen_binary_method(op: Ops) -> Callable[[ArgType, ArgType], BinaryOp]:
    def logic_operation(self: ArgType, other: ArgType) -> BinaryOp:
        return cast(Type[BinaryOp], LogicExp.factory(op))(self, other)

    return logic_operation


for clas, enum in cls_enum_pairs:
    setattr(clas, "__and__", gen_binary_method(enum.AND))
    setattr(clas, "__rand__", gen_binary_method(enum.AND))
    setattr(clas, "__or__", gen_binary_method(enum.OR))
    setattr(clas, "__ror__", gen_binary_method(enum.OR))
    setattr(clas, "__xor__", gen_binary_method(enum.XOR))
    setattr(clas, "__rxor__", gen_binary_method(enum.XOR))
