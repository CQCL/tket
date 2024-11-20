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

from dataclasses import dataclass
from typing import Any

from pytket.circuit import (
    Bit,
    BitRegister,
    Circuit,
    ClBitVar,
    ClExpr,
    ClExprOp,
    ClOp,
    ClRegVar,
    OpType,
    WiredClExpr,
)
from pytket.circuit.logic_exp import BitWiseOp, LogicExp, Ops, RegWiseOp

_reg_output_clops = set(
    [
        ClOp.RegAnd,
        ClOp.RegOr,
        ClOp.RegXor,
        ClOp.RegNot,
        ClOp.RegZero,
        ClOp.RegOne,
        ClOp.RegAdd,
        ClOp.RegSub,
        ClOp.RegMul,
        ClOp.RegDiv,
        ClOp.RegPow,
        ClOp.RegLsh,
        ClOp.RegRsh,
        ClOp.RegNeg,
    ]
)


def has_reg_output(op: ClOp) -> bool:
    return op in _reg_output_clops


def clop_from_ops(op: Ops) -> ClOp:
    match op:
        case BitWiseOp.AND:
            return ClOp.BitAnd
        case BitWiseOp.OR:
            return ClOp.BitOr
        case BitWiseOp.XOR:
            return ClOp.BitXor
        case BitWiseOp.EQ:
            return ClOp.BitEq
        case BitWiseOp.NEQ:
            return ClOp.BitNeq
        case BitWiseOp.NOT:
            return ClOp.BitNot
        case BitWiseOp.ZERO:
            return ClOp.BitZero
        case BitWiseOp.ONE:
            return ClOp.BitOne
        case RegWiseOp.AND:
            return ClOp.RegAnd
        case RegWiseOp.OR:
            return ClOp.RegOr
        case RegWiseOp.XOR:
            return ClOp.RegXor
        case RegWiseOp.EQ:
            return ClOp.RegEq
        case RegWiseOp.NEQ:
            return ClOp.RegNeq
        case RegWiseOp.LT:
            return ClOp.RegLt
        case RegWiseOp.GT:
            return ClOp.RegGt
        case RegWiseOp.LEQ:
            return ClOp.RegLeq
        case RegWiseOp.GEQ:
            return ClOp.RegGeq
        case RegWiseOp.ADD:
            return ClOp.RegAdd
        case RegWiseOp.SUB:
            return ClOp.RegSub
        case RegWiseOp.MUL:
            return ClOp.RegMul
        case RegWiseOp.DIV:
            return ClOp.RegDiv
        case RegWiseOp.POW:
            return ClOp.RegPow
        case RegWiseOp.LSH:
            return ClOp.RegLsh
        case RegWiseOp.RSH:
            return ClOp.RegRsh
        case RegWiseOp.NOT:
            return ClOp.RegNot
        case RegWiseOp.NEG:
            return ClOp.RegNeg


@dataclass
class ExpressionConverter:
    bit_indices: dict[Bit, int]
    reg_indices: dict[BitRegister, int]

    def convert(self, exp: LogicExp) -> ClExpr:
        op: ClOp = clop_from_ops(exp.op)
        args: list[int | ClBitVar | ClRegVar | ClExpr] = []
        for arg in exp.args:
            if isinstance(arg, LogicExp):
                args.append(self.convert(arg))
            elif isinstance(arg, Bit):
                args.append(ClBitVar(self.bit_indices[arg]))
            elif isinstance(arg, BitRegister):
                args.append(ClRegVar(self.reg_indices[arg]))
            else:
                assert isinstance(arg, int)
                args.append(arg)
        return ClExpr(op, args)


def wired_clexpr_from_logic_exp(
    exp: LogicExp, output_bits: list[Bit]
) -> tuple[WiredClExpr, list[Bit]]:
    """Convert a :py:class:`LogicExp` to a :py:class:`WiredClExpr`

    :param exp: the LogicExp
    :param output_bits: list of output bits of the LogicExp
    :return: the WiredClExpr and its full list of arguments
    """
    # 1. Construct lists of input bits and registers (where the positions of the items
    # in each list will be the indices of the corresponding variables in the ClExpr):
    all_vars = exp.all_inputs_ordered()
    input_bits: list[Bit] = [var for var in all_vars if isinstance(var, Bit)]
    input_regs: list[BitRegister] = [
        var for var in all_vars if isinstance(var, BitRegister)
    ]
    # 2. Order the arguments: first the input bits, then all the bits in the input
    # registers then any remaining output bits:
    args = []
    args.extend(input_bits)
    for r in input_regs:
        args.extend(r.to_list())
    args.extend(b for b in output_bits if b not in args)
    # 3. Construct the WiredClExpr and return it with the argument list:
    return (
        WiredClExpr(
            ExpressionConverter(
                {b: i for i, b in enumerate(input_bits)},
                {r: i for i, r in enumerate(input_regs)},
            ).convert(exp),
            {i: args.index(b) for i, b in enumerate(input_bits)},
            {i: [args.index(b) for b in r.to_list()] for i, r in enumerate(input_regs)},
            [args.index(b) for b in output_bits],
        ),
        args,
    )


def check_register_alignments(circ: Circuit) -> bool:
    """Check whether all `ClExprOp` operations in the circuit are register-aligned.

    This means that all register variables and outputs occurring in `ClExprOp` comprise
    whole registers with the bits in the correct order.

    :param circ: circuit to check
    :return: True iff all `ClExprOp` operations are register-aligned
    """
    cregs: set[tuple[Bit, ...]] = set(
        tuple(creg.to_list()) for creg in circ.c_registers
    )
    for cmd in circ:
        op = cmd.op
        if op.type == OpType.ClExpr:
            assert isinstance(op, ClExprOp)
            wexpr: WiredClExpr = op.expr
            args = cmd.args
            if any(
                tuple(args[i] for i in poslist) not in cregs
                for poslist in wexpr.reg_posn.values()
            ) or (
                has_reg_output(wexpr.expr.op)
                and tuple(args[i] for i in wexpr.output_posn) not in cregs
            ):
                return False
    return True


def _add_clexpr_to_circuit_from_logicexp(
    circ: Circuit, exp: LogicExp, output_bits: list[Bit], **kwargs: Any
) -> None:
    wexpr, args = wired_clexpr_from_logic_exp(exp, output_bits)
    circ.add_clexpr(wexpr, args, **kwargs)
