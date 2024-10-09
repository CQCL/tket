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

from pytket.circuit import Circuit, ClOp, ClBitVar, ClRegVar, ClExpr, WiredClExpr


def test_op() -> None:
    reg_add = ClOp.RegAdd
    assert str(reg_add) == "ClOp.RegAdd"
    reg_sub = ClOp.RegSub
    assert reg_add != reg_sub


def test_vars() -> None:
    bvar3 = ClBitVar(3)
    assert bvar3.i == 3
    assert str(bvar3) == "b3"
    bvar4 = ClBitVar(4)
    rvar3 = ClRegVar(3)
    assert rvar3.i == 3
    assert str(rvar3) == "r3"
    rvar3a = ClRegVar(3)
    assert bvar3 != bvar4
    assert bvar3 != rvar3
    assert rvar3 == rvar3a


def test_expr() -> None:
    b0 = ClBitVar(0)
    r0 = ClRegVar(0)
    three = 3
    expr0 = ClExpr(op=ClOp.RegAdd, args=[r0, three])
    expr = ClExpr(op=ClOp.BitXor, args=[expr0, b0])
    assert str(expr) == "xor(add(r0, 3), b0)"
    assert expr.op == ClOp.BitXor
    args = expr.args
    assert len(args) == 2
    assert args[0] == expr0
    assert args[1] == b0


def test_wexpr() -> None:
    expr = ClExpr(
        op=ClOp.RegDiv,
        args=[ClRegVar(0), ClExpr(op=ClOp.RegAdd, args=[2, ClBitVar(0)])],
    )
    wexpr = WiredClExpr(
        expr=expr, bit_posn={0: 1}, reg_posn={0: [2, 0]}, output_posn=[2, 0]
    )
    assert str(wexpr) == "div(r0, add(2, b0)) [b0:1, r0:(2,0) --> (2,0)]"
    assert wexpr.expr == expr
    assert wexpr.bit_posn == {0: 1}
    assert wexpr.reg_posn == {0: [2, 0]}
    assert wexpr.output_posn == [2, 0]


def test_adding_to_circuit() -> None:
    expr = ClExpr(op=ClOp.BitXor, args=[ClBitVar(0), ClBitVar(1)])
    wexpr = WiredClExpr(expr=expr, bit_posn={0: 0, 1: 1}, output_posn=[2])
    c = Circuit(0, 3)
    c.add_clexpr(wexpr, c.bits)
    cmds = c.get_commands()
    assert len(cmds) == 1
    assert cmds[0].op.expr == wexpr
    d = c.to_dict()
    c1 = Circuit.from_dict(d)
    assert c == c1
    d1 = c1.to_dict()
    assert d == d1
