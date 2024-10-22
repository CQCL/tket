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

from pytket.circuit import (
    Bit,
    CircBox,
    Circuit,
    ClBitVar,
    ClExpr,
    ClExprOp,
    ClOp,
    ClRegVar,
    OpType,
    WiredClExpr,
)
from pytket.qasm import circuit_to_qasm_str, circuit_from_qasm_str


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
    wexpr_dict = wexpr.to_dict()
    assert wexpr_dict == {
        "bit_posn": [[0, 1]],
        "expr": {
            "args": [
                {
                    "input": {
                        "term": {"type": "reg", "var": {"index": 0}},
                        "type": "var",
                    },
                    "type": "term",
                },
                {
                    "input": {
                        "args": [
                            {"input": {"term": 2, "type": "int"}, "type": "term"},
                            {
                                "input": {
                                    "term": {"type": "bit", "var": {"index": 0}},
                                    "type": "var",
                                },
                                "type": "term",
                            },
                        ],
                        "op": "RegAdd",
                    },
                    "type": "expr",
                },
            ],
            "op": "RegDiv",
        },
        "output_posn": [2, 0],
        "reg_posn": [[0, [2, 0]]],
    }
    wexpr1 = WiredClExpr.from_dict(wexpr_dict)
    assert wexpr == wexpr1


def test_adding_to_circuit() -> None:
    expr = ClExpr(op=ClOp.BitXor, args=[ClBitVar(0), ClBitVar(1)])
    wexpr = WiredClExpr(expr=expr, bit_posn={0: 0, 1: 1}, output_posn=[2])
    c = Circuit(0, 3)
    c.add_clexpr(wexpr, c.bits)
    cmds = c.get_commands()
    assert len(cmds) == 1
    op = cmds[0].op
    assert isinstance(op, ClExprOp)
    assert op.expr == wexpr
    d = c.to_dict()
    c1 = Circuit.from_dict(d)
    assert c == c1
    d1 = c1.to_dict()
    assert d == d1
    c2 = c.copy()
    assert c2 == c


def test_qasm_conversion() -> None:
    c = Circuit()
    c.add_c_register("a", 3)
    c.add_c_register("b", 3)
    c.add_c_register("c", 3)
    c.add_c_register("d", 3)
    expr = ClExpr(
        op=ClOp.RegSub,
        args=[
            ClExpr(
                op=ClOp.RegDiv,
                args=[ClExpr(ClOp.RegAdd, args=[ClRegVar(0), ClRegVar(1)]), 2],
            ),
            ClRegVar(2),
        ],
    )
    wexpr = WiredClExpr(
        expr=expr,
        reg_posn={0: [0, 1, 2], 1: [3, 4, 5], 2: [6, 7, 8]},
        output_posn=[9, 10, 11],
    )
    c.add_clexpr(wexpr, c.bits)
    qasm = circuit_to_qasm_str(c, header="hqslib1")
    assert (
        qasm
        == """OPENQASM 2.0;
include "hqslib1.inc";

creg a[3];
creg b[3];
creg c[3];
creg d[3];
d = (((a + b) / 2) - c);
"""
    )
    c1 = circuit_from_qasm_str(qasm, use_clexpr=True)
    assert c == c1


def make_circ() -> Circuit:
    c = Circuit()
    c.add_bit(Bit("x", 0))
    c.add_bit(Bit("x", 1))
    c.add_bit(Bit("y", 0))
    c.add_clexpr(
        WiredClExpr(
            expr=ClExpr(op=ClOp.BitXor, args=[ClBitVar(0), ClBitVar(1)]),
            bit_posn={0: 0, 1: 1},
            output_posn=[2],
        ),
        [Bit("x", 0), Bit("x", 1), Bit("y", 0)],
    )
    return c


def test_copy_and_flatten() -> None:
    # See https://github.com/CQCL/tket/issues/1544
    c0 = make_circ()
    c1 = make_circ()
    assert c0 == c1
    c2 = c1.copy()
    c2.flatten_registers()
    assert c0 == c1
    assert c2.get_commands()[0].op == c0.get_commands()[0].op
    qasm = circuit_to_qasm_str(c2, header="hqslib1")
    assert (
        qasm
        == """OPENQASM 2.0;
include "hqslib1.inc";

creg c[3];
c[2] = (c[0] ^ c[1]);
"""
    )


def test_circbox() -> None:
    # See https://github.com/CQCL/tket/issues/1544
    c0 = make_circ()
    cbox = CircBox(c0)
    c1 = Circuit(0, 3)
    c1.add_circbox(cbox, [0, 1, 2])
    c2 = c1.copy()
    c2.flatten_registers()
    assert c1 == c2
