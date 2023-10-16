# Copyright 2019-2023 Cambridge Quantum Computing
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
from typing import (
    TYPE_CHECKING,
    Any,
    Tuple,
    Type,
    Union,
    Callable,
    Optional,
    Sequence,
)

from pytket._tket.circuit import *
from pytket._tket.circuit import Circuit

from pytket._tket.unit_id import *
from pytket._tket.unit_id import Bit, BitRegister

# prefixes for assertion bits
from pytket._tket.unit_id import _DEBUG_ZERO_REG_PREFIX, _DEBUG_ONE_REG_PREFIX
from pytket._tket.pauli import Pauli

from pytket import wasm

from .logic_exp import (
    BitLogicExp,
    BitWiseOp,
    RegLogicExp,
    RegWiseOp,
    Constant,
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


def overload_add_wasm(
    self: Circuit,
    funcname: str,
    filehandler: wasm.WasmFileHandler,
    list_i: Sequence[int],
    list_o: Sequence[int],
    args: Union[Sequence[int], Sequence[Bit]],
    args_wasm: Optional[Sequence[int]] = None,
    **kwargs: Any,
) -> Circuit:
    """Add a classical function call from a wasm file to the circuit.
    \n\n:param funcname: name of the function that is called
    \n:param filehandler: wasm file handler to identify the wasm file
    \n:param list_i: list of the number of bits in the input variables
    \n:param list_o: list of the number of bits in the output variables
    \n:param args: vector of circuit bits the wasm op should be added to
    \n:param args_wasm: vector of wasmstates the wasm op should be added to
    \n:param kwargs: additional arguments passed to `add_gate_method` .
     Allowed parameters are `opgroup`,  `condition` , `condition_bits`,
     `condition_value`
    \n:return: the new :py:class:`Circuit`"""

    if args_wasm is None:
        args_wasm = [0]

    for x in list_i:
        if x > filehandler._int_size:
            raise ValueError(
                f"only functions with i{filehandler._int_size} type are allowed"
            )

    for x in list_o:
        if x > filehandler._int_size:
            raise ValueError(
                f"only functions with i{filehandler._int_size} type are allowed"
            )

    if filehandler.check_function(funcname, len(list_i), len(list_o)):
        if (len(args_wasm)) > 0:
            self._add_w_register(max(args_wasm) + 1)
        return self._add_wasm(
            funcname, str(filehandler), list_i, list_o, args, args_wasm, **kwargs
        )

    raise ValueError(f"{funcname} not found, check {repr(filehandler)}")


setattr(Circuit, "add_wasm", overload_add_wasm)


def overload_add_wasm_to_reg(
    self: Circuit,
    funcname: str,
    filehandler: wasm.WasmFileHandler,
    list_i: Sequence[BitRegister],
    list_o: Sequence[BitRegister],
    args_wasm: Optional[Sequence[int]] = None,
    **kwargs: Any,
) -> Circuit:
    """Add a classical function call from a wasm file to the circuit.
    \n\n:param funcname: name of the function that is called
    \n:param filehandler: wasm file handler to identify the wasm file
    \n:param list_i: list of the classical registers assigned to
     the input variables of the function call
    \n:param list_o: list of the classical registers assigned to
     the output variables of the function call
    \n:param args_wasm: vector of wasmstates the wasm op should be added to
    \n:param kwargs: additional arguments passed to `add_gate_method` .
     Allowed parameters are `opgroup`,  `condition` , `condition_bits`,
     `condition_value`
    \n:return: the new :py:class:`Circuit`"""

    if args_wasm is None:
        args_wasm = [0]

    if filehandler.check_function(funcname, len(list_i), len(list_o)):
        if (len(args_wasm)) > 0:
            self._add_w_register(max(args_wasm) + 1)
        return self._add_wasm(
            funcname, str(filehandler), list_i, list_o, args_wasm, **kwargs
        )

    raise ValueError(f"{funcname} not found, check {repr(filehandler)}")


setattr(Circuit, "add_wasm_to_reg", overload_add_wasm_to_reg)


# overload operators for Bit, BitRegister and expressions over these
# such that the operation returns a LogicExp describing the operation
cls_enum_pairs_bit: Tuple[Tuple[Type, Union[Type[BitWiseOp], Type[RegWiseOp]]], ...] = (
    (Bit, BitWiseOp),
)

cls_enum_pairs_reg: Tuple[Tuple[Type, Union[Type[BitWiseOp], Type[RegWiseOp]]], ...] = (
    (BitRegister, RegWiseOp),
)

BitArgType = Union[BitLogicExp, Bit, Constant]
RegArgType = Union[RegLogicExp, BitRegister, Constant]


def gen_binary_method_bit(op: Ops) -> Callable[[BitArgType, BitArgType], BitLogicExp]:
    def logic_operation(self: BitArgType, other: BitArgType) -> BitLogicExp:
        return BitLogicExp(op, [self, other])

    return logic_operation


def gen_binary_method_reg(op: Ops) -> Callable[[RegArgType, RegArgType], RegLogicExp]:
    def logic_operation(self: RegArgType, other: RegArgType) -> RegLogicExp:
        return RegLogicExp(op, [self, other])

    return logic_operation


for clas, enum in cls_enum_pairs_bit:
    setattr(clas, "__and__", gen_binary_method_bit(enum.AND))
    setattr(clas, "__rand__", gen_binary_method_bit(enum.AND))
    setattr(clas, "__or__", gen_binary_method_bit(enum.OR))
    setattr(clas, "__ror__", gen_binary_method_bit(enum.OR))
    setattr(clas, "__xor__", gen_binary_method_bit(enum.XOR))
    setattr(clas, "__rxor__", gen_binary_method_bit(enum.XOR))

for clas, enum in cls_enum_pairs_reg:
    setattr(clas, "__and__", gen_binary_method_reg(enum.AND))
    setattr(clas, "__rand__", gen_binary_method_reg(enum.AND))
    setattr(clas, "__or__", gen_binary_method_reg(enum.OR))
    setattr(clas, "__ror__", gen_binary_method_reg(enum.OR))
    setattr(clas, "__xor__", gen_binary_method_reg(enum.XOR))
    setattr(clas, "__rxor__", gen_binary_method_reg(enum.XOR))
    setattr(clas, "__add__", gen_binary_method_reg(RegWiseOp.ADD))
    setattr(clas, "__sub__", gen_binary_method_reg(RegWiseOp.SUB))
    setattr(clas, "__mul__", gen_binary_method_reg(RegWiseOp.MUL))
    setattr(clas, "__floordiv__", gen_binary_method_reg(RegWiseOp.DIV))
    setattr(clas, "__pow__", gen_binary_method_reg(RegWiseOp.POW))
    setattr(clas, "__lshift__", gen_binary_method_reg(RegWiseOp.LSH))
    setattr(clas, "__rshift__", gen_binary_method_reg(RegWiseOp.RSH))
