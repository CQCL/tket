# Copyright Quantinuum
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

from typing import Callable, Union

from pytket._tket.unit_id import *
from pytket._tket.unit_id import (
    _TEMP_BIT_NAME,
    _TEMP_BIT_REG_BASE,
    Bit,
    BitRegister,
    Qubit,
    QubitRegister,
)
from pytket.circuit.logic_exp import (
    BitLogicExp,
    BitWiseOp,
    Constant,
    LogicExp,
    RegLogicExp,
    RegWiseOp,
    create_bit_logic_exp,
    create_reg_logic_exp,
)


def _bitregister_next(self: BitRegister) -> Bit:
    if self._current < self.size:
        result = self[self._current]
        self._current += 1
        return result
    raise StopIteration


def _qubitregister_next(self: QubitRegister) -> Qubit:
    if self._current < self.size:
        result = self[self._current]
        self._current += 1
        return result
    raise StopIteration


setattr(BitRegister, "__next__", _bitregister_next)
BitRegister.__next__.__name__ = "__next__"

setattr(QubitRegister, "__next__", _qubitregister_next)
QubitRegister.__next__.__name__ = "__next__"

# overload operators for Bit, BitRegister and expressions over these
# such that the operation returns a LogicExp describing the operation

BitArgType = Union[LogicExp, Bit, Constant]
RegArgType = Union[LogicExp, BitRegister, Constant]


def gen_binary_method_bit(
    op: BitWiseOp, name: str
) -> Callable[[BitArgType, BitArgType], BitLogicExp]:
    def logic_operation(self: BitArgType, other: BitArgType) -> BitLogicExp:
        return create_bit_logic_exp(op, [self, other])

    logic_operation.__name__ = name
    return logic_operation


def gen_binary_method_reg(
    op: RegWiseOp, name: str
) -> Callable[[RegArgType, RegArgType], RegLogicExp]:
    def logic_operation(self: RegArgType, other: RegArgType) -> RegLogicExp:
        return create_reg_logic_exp(op, [self, other])

    logic_operation.__name__ = name
    return logic_operation


setattr(Bit, "__and__", gen_binary_method_bit(BitWiseOp.AND, "__and__"))
setattr(Bit, "__rand__", gen_binary_method_bit(BitWiseOp.AND, "__rand__"))
setattr(Bit, "__or__", gen_binary_method_bit(BitWiseOp.OR, "__or__"))
setattr(Bit, "__ror__", gen_binary_method_bit(BitWiseOp.OR, "__ror__"))
setattr(Bit, "__xor__", gen_binary_method_bit(BitWiseOp.XOR, "__xor__"))
setattr(Bit, "__rxor__", gen_binary_method_bit(BitWiseOp.XOR, "__rxor__"))
setattr(BitRegister, "__and__", gen_binary_method_reg(RegWiseOp.AND, "__and__"))
setattr(BitRegister, "__rand__", gen_binary_method_reg(RegWiseOp.AND, "__rand__"))
setattr(BitRegister, "__or__", gen_binary_method_reg(RegWiseOp.OR, "__or__"))
setattr(BitRegister, "__ror__", gen_binary_method_reg(RegWiseOp.OR, "__ror__"))
setattr(BitRegister, "__xor__", gen_binary_method_reg(RegWiseOp.XOR, "__xor__"))
setattr(BitRegister, "__rxor__", gen_binary_method_reg(RegWiseOp.XOR, "__rxor__"))
setattr(BitRegister, "__add__", gen_binary_method_reg(RegWiseOp.ADD, "__add__"))
setattr(BitRegister, "__sub__", gen_binary_method_reg(RegWiseOp.SUB, "__sub__"))
setattr(BitRegister, "__mul__", gen_binary_method_reg(RegWiseOp.MUL, "__mul__"))
setattr(
    BitRegister, "__floordiv__", gen_binary_method_reg(RegWiseOp.DIV, "__floordiv__")
)
setattr(BitRegister, "__pow__", gen_binary_method_reg(RegWiseOp.POW, "__pow__"))
setattr(BitRegister, "__lshift__", gen_binary_method_reg(RegWiseOp.LSH, "__lshift__"))
setattr(BitRegister, "__rshift__", gen_binary_method_reg(RegWiseOp.RSH, "__rshift__"))
