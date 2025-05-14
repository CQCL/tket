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

"""The circuit module provides an API to interact with the
tket :py:class:`~.Circuit` data structure.
This module is provided in binary form during the PyPI installation."""

from collections.abc import Callable, Sequence
from typing import (
    Any,
    Optional,
    Union,
)

from pytket import wasm
from pytket._tket.circuit import *
from pytket._tket.circuit import Circuit
from pytket._tket.pauli import Pauli
from pytket._tket.unit_id import *

# prefixes for assertion bits
from pytket._tket.unit_id import (
    _DEBUG_ONE_REG_PREFIX,
    _DEBUG_ZERO_REG_PREFIX,
    Bit,
    BitRegister,
)

from .clexpr import wired_clexpr_from_logic_exp
from .logic_exp import (
    BinaryOp,
    LogicExp,
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


def add_wasm(  # noqa: PLR0913
    self: Circuit,
    funcname: str,
    filehandler: wasm.WasmModuleHandler,
    list_i: Sequence[int],
    list_o: Sequence[int],
    args: Union[Sequence[int], Sequence[Bit]],  # noqa: UP007
    args_wasm: Optional[Sequence[int]] = None,  # noqa: UP007
    **kwargs: Any,
) -> Circuit:
    """Add a classical function call from a wasm file to the circuit.
    \n\n:param funcname: name of the function that is called
    \n:param filehandler: wasm file or module handler to identify the wasm module
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
        if x > filehandler._int_size:  # noqa: SLF001
            raise ValueError(
                f"only functions with i{filehandler._int_size} type are allowed"  # noqa: SLF001
            )

    for x in list_o:
        if x > filehandler._int_size:  # noqa: SLF001
            raise ValueError(
                f"only functions with i{filehandler._int_size} type are allowed"  # noqa: SLF001
            )

    if filehandler.check_function(funcname, len(list_i), len(list_o)):
        if (len(args_wasm)) > 0:
            self._add_w_register(max(args_wasm) + 1)
        return self._add_wasm(
            funcname, str(filehandler), list_i, list_o, args, args_wasm, **kwargs
        )

    raise ValueError(f"{funcname} not found, check {filehandler!r}")


setattr(Circuit, "add_wasm", add_wasm)  # noqa: B010


def add_wasm_to_reg(  # noqa: PLR0913
    self: Circuit,
    funcname: str,
    filehandler: wasm.WasmModuleHandler,
    list_i: Sequence[BitRegister],
    list_o: Sequence[BitRegister],
    args_wasm: Optional[Sequence[int]] = None,  # noqa: UP007
    **kwargs: Any,
) -> Circuit:
    """Add a classical function call from a wasm file to the circuit.
    \n\n:param funcname: name of the function that is called
    \n:param filehandler: wasm file or module handler to identify the wasm module
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

    if filehandler.checked:
        for reg in list_i:
            if reg.size > 32:  # noqa: PLR2004
                raise ValueError(
                    """wasm is only supporting 32 bit size registers,
please use only registers of at most 32 bits"""
                )

        for reg in list_o:
            if reg.size > 32:  # noqa: PLR2004
                raise ValueError(
                    """wasm is only supporting 32 bit size registers,
please use only registers of at most 32 bits"""
                )

    # If the filehandler has not been checked we allow it to
    # be added without checking the function arity.
    if not filehandler.checked or filehandler.check_function(
        funcname, len(list_i), len(list_o)
    ):
        if (len(args_wasm)) > 0:
            self._add_w_register(max(args_wasm) + 1)
        return self._add_wasm(
            funcname, str(filehandler), list_i, list_o, args_wasm, **kwargs
        )

    raise ValueError(f"{funcname} not found, check {filehandler!r}")


setattr(Circuit, "add_wasm_to_reg", add_wasm_to_reg)  # noqa: B010


def add_clexpr_from_logicexp(
    circ: Circuit, exp: LogicExp, output_bits: list[Bit], **kwargs: Any
) -> Circuit:
    """Append a :py:class:`~.ClExprOp` defined in terms of a logical expression.
    \n\nExample:
    \n>>> c = Circuit()\n>>> x_reg = c.add_c_register('x', 3)\n>>> y_reg = c.add_c_register('y', 3)\n>>> z_reg = c.add_c_register('z', 3)\n>>> c.add_clexpr_from_logicexp(x_reg | y_reg, z_reg.to_list())\n[ClExpr x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]; ]
    \n:param exp: logical expression
    \n:param output_bits: list of bits in output
    \n:return: the updated circuit"""
    wexpr, args = wired_clexpr_from_logic_exp(exp, output_bits)
    circ.add_clexpr(wexpr, args, **kwargs)
    return circ


setattr(Circuit, "add_clexpr_from_logicexp", add_clexpr_from_logicexp)  # noqa: B010
