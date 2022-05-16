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

import json
from jsonschema import validate  # type: ignore
from pathlib import Path
import pickle

from pytket.circuit import (  # type: ignore
    Circuit,
    Op,
    OpType,
    Command,
    fresh_symbol,
    CircBox,
    Unitary1qBox,
    Unitary2qBox,
    Unitary3qBox,
    ExpBox,
    PauliExpBox,
    QControlBox,
    PhasePolyBox,
    CustomGateDef,
    CustomGate,
    Qubit,
    Bit,
    BitRegister,
    QubitRegister,
)
from pytket.circuit.display import render_circuit_as_html

from pytket.pauli import Pauli  # type: ignore
from pytket.passes import PauliSimp, CliffordSimp, SynthesiseTket, DecomposeBoxes, RemoveRedundancies  # type: ignore
from pytket.predicates import CompilationUnit  # type: ignore
from pytket.transform import Transform, PauliSynthStrat  # type: ignore

import numpy as np
import sympy  # type: ignore
from sympy import Symbol, pi, sympify, functions  # type: ignore
from math import sqrt

import pytest  # type: ignore

from hypothesis import given, settings
import strategies as st  # type: ignore


curr_file_path = Path(__file__).resolve().parent


with open(curr_file_path.parent.parent / "schemas/circuit_v1.json", "r") as f:
    schema = json.load(f)


def test_op_free_symbols() -> None:
    c = Circuit(2)
    c.add_barrier([0, 1])
    op = c.get_commands()[0].op
    assert op.free_symbols() == set()
    alpha = Symbol("alpha")  # type: ignore
    c.Rx(alpha, 0)
    op = c.get_commands()[1].op
    assert op.free_symbols() == {alpha}


def test_circuit_transpose() -> None:
    c = Circuit(2)
    u = np.asarray([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0]])
    ubox = Unitary2qBox(u)
    c.add_unitary2qbox(ubox, 1, 0)
    c.CX(0, 1)
    c_t = c.transpose()
    commands = c_t.get_commands()
    assert str(commands[0]) == "CX q[0], q[1];"
    assert str(commands[1]) == "Unitary2qBox q[1], q[0];"
    assert commands[1].op.get_matrix().all() == u.transpose().all()


def test_circuit_dagger() -> None:
    c = Circuit(2)
    u = np.asarray([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1j, 0], [1, 0, 0, 0]])
    ubox = Unitary2qBox(u)
    c.add_unitary2qbox(ubox, 1, 0)
    c.add_gate(OpType.CnRy, 0.3, [0, 1])
    c_d = c.dagger()
    commands = c_d.get_commands()
    assert str(commands[0]) == "CnRy(3.7) q[0], q[1];"
    assert commands[0].qubits == [Qubit(0), Qubit(1)]
    assert str(commands[1]) == "Unitary2qBox q[1], q[0];"
    assert commands[1].op.get_matrix().all() == u.conj().transpose().all()


# TKET-1365 bug
def test_cnx_dagger() -> None:
    c = Circuit(2)
    c.add_gate(OpType.CnX, [0, 1])
    for d in (c.dagger(), c.transpose()):
        d_dict = d.to_dict()

        args = d_dict["commands"][0]["args"]
        # before bugfix args would be None, due to empty signature
        assert args == [["q", [0]], ["q", [1]]]


def test_circuit_name() -> None:
    name = "test"
    a = Circuit(name)
    b = Circuit(2, name=name)
    c = Circuit(2, 2, name)

    assert a.name == b.name == c.name

    a.name = "new_test"
    assert a.name == "new_test"
    assert a.name != b.name

    c = Circuit(4, 4, name)
    c.X(0)
    c.H(1)
    c.S(1)
    c.CX(2, 0)
    c.CRz(0.5, 0, 3)
    c.Measure(3, 3)
    c.Measure(1, 1)

    SynthesiseTket().apply(c)
    assert c.name == name


def test_circuit_gen() -> None:
    c = Circuit(4, 4)
    c.X(0)
    c.H(1)
    c.S(1)
    c.CX(2, 0)
    c.CRz(0.5, 0, 3)
    c.CRx(0.5, 0, 3)
    c.CRy(0.5, 0, 3)
    c.CV(0, 3)
    c.CVdg(1, 3)
    c.CSX(0, 2)
    c.CSXdg(1, 2)
    c.SX(3)
    c.SXdg(0)
    c.Measure(3, 3)
    c.Measure(1, 1)

    assert c.n_qubits == 4
    assert c._n_vertices() == 31
    assert c.n_gates == 15

    commands = c.get_commands()
    assert len(commands) == 15
    assert str(commands[0]) == "X q[0];"
    assert str(commands[2]) == "CX q[2], q[0];"
    assert str(commands[4]) == "CRz(0.5) q[0], q[3];"
    assert str(commands[5]) == "CRx(0.5) q[0], q[3];"
    assert str(commands[6]) == "CRy(0.5) q[0], q[3];"
    assert str(commands[7]) == "CV q[0], q[3];"
    assert str(commands[8]) == "CSX q[0], q[2];"
    assert str(commands[9]) == "CVdg q[1], q[3];"
    assert str(commands[10]) == "SXdg q[0];"
    assert str(commands[11]) == "CSXdg q[1], q[2];"
    assert str(commands[12]) == "SX q[3];"
    assert str(commands[13]) == "Measure q[1] --> c[1];"
    assert str(commands[14]) == "Measure q[3] --> c[3];"

    assert commands[14].qubits == [Qubit(3)]
    assert commands[14].bits == [Bit(3)]
    assert c.depth_by_type({OpType.CX, OpType.CRz}) == 2
    assert commands[0] == Command(Op.create(OpType.X), [Qubit(0)])


def test_circuit_gen_ids() -> None:
    c = Circuit()
    a = [Qubit("a", i) for i in range(4)]
    b = [Bit("b", i) for i in range(4)]
    for q in a:
        c.add_qubit(q)
    for cb in b:
        c.add_bit(cb)
    c.add_gate(OpType.X, [a[0]])
    c.H(a[1])
    c.add_gate(OpType.Rz, 0.5, [a[1]])
    c.CX(a[2], a[0])
    c.CRz(0.5, a[0], a[3])
    c.Measure(a[3], b[3])
    c.Measure(a[1], b[1])

    assert c.n_qubits == 4
    assert c._n_vertices() == 23
    assert c.n_gates == 7

    commands = c.get_commands()
    assert len(commands) == 7
    assert str(commands[0]) == "X a[0];"
    assert str(commands[2]) == "CX a[2], a[0];"
    assert str(commands[3]) == "Rz(0.5) a[1];"
    assert str(commands[6]) == "Measure a[3] --> b[3];"


def test_symbolic_ops() -> None:
    c = Circuit(2)
    alpha = Symbol("alpha")  # type: ignore
    c.Rx(alpha, 0)
    beta = fresh_symbol("alpha")
    c.CRz(beta * 2, 1, 0)
    s_map = {alpha: 0.5, beta: 3.2}
    assert c.is_symbolic()
    assert c.free_symbols() == set([alpha, beta])
    c.symbol_substitution(s_map)
    assert not c.is_symbolic()

    commands = c.get_commands()
    assert beta.__str__() == "alpha_1"
    assert np.allclose(commands[0].op.params, [0.5], atol=1e-10)
    assert np.allclose(commands[1].op.params, [2.4], atol=1e-10)


def test_subst_4() -> None:
    # https://github.com/CQCL/tket/issues/219
    m = fresh_symbol("m")
    c = Circuit(1)
    a = m / 4
    c.add_gate(OpType.Rx, a, [0])
    c.symbol_substitution({m: 4})
    angle = c.get_commands()[0].op.params[0]
    print(angle)
    assert np.isclose(angle, 1.0)


def test_sympy_conversion() -> None:
    def get_type_tree(expr: sympy.Expr) -> str:
        # Format e.g. "<class 'sympy.core.numbers.Pi'>" to "Pi"
        tree_str = str(type(expr)).rsplit(".", 1)[-1].split("'")[0]
        if len(expr.args) != 0:
            tree_str += " ("
            tree_str += ", ".join([get_type_tree(a) for a in expr.args])  # type: ignore
            tree_str += ")"
        return tree_str

    test_dict = {
        "Abs(x**2)": "Abs (Pow (Symbol, Integer))",
        "atan2(0.436, pi*I + b)": "atan2 (Float, Add (Symbol, Mul (Pi, ImaginaryUnit)))",
        "cos(a)": "cos (Symbol)",
        "a*oo": "Mul (Infinity, Symbol)",
        "a*-oo": "Mul (NegativeInfinity, Symbol)",
        "a*zoo": "Mul (ComplexInfinity, Symbol)",
        "a*E": "Mul (Exp1, Symbol)",
        "atan2(a, nan)": "atan2 (Symbol, NaN)",
    }
    for expr_string, type_tree in test_dict.items():
        c = Circuit(1)
        c.Rz(sympify(expr_string), 0)  # type: ignore
        com = c.get_commands()[0]
        assert get_type_tree(com.op.params[0]) == type_tree


def test_4x4_matrix_to_circ() -> None:
    u = np.asarray([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
    ubox = Unitary2qBox(u)
    c = Circuit(2)
    c.add_unitary2qbox(ubox, 0, 1)
    Transform.DecomposeBoxes().apply(c)
    Transform.OptimiseCliffords().apply(c)
    mat = c.get_unitary()
    assert np.allclose(mat, u)
    assert c.n_gates_of_type(OpType.CX) == 1


def test_8x8_matrix_to_circ() -> None:
    u = np.asarray(
        [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, -1],
        ]
    )
    ubox = Unitary3qBox(u)
    c = Circuit(3)
    c.add_unitary3qbox(ubox, 0, 1, 2)
    Transform.DecomposeBoxes().apply(c)
    Transform.OptimiseCliffords().apply(c)
    mat = c.get_unitary()
    assert np.allclose(mat, u)
    assert c.n_gates_of_type(OpType.CX) <= 10


def test_exp_to_circ() -> None:
    PI = float(pi.evalf())
    u = np.asarray([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]) * -PI / 4
    ebox = ExpBox(u, 1.0)
    c = Circuit(2)
    c.add_expbox(ebox, 0, 1)
    Transform.DecomposeBoxes().apply(c)
    mat = c.get_unitary()
    correct = np.asarray(
        [[1 - 1j, 0, 0, 0], [0, 1 + 1j, 0, 0], [0, 0, 1 + 1j, 0], [0, 0, 0, 1 - 1j]]
    ) / np.sqrt(2)
    assert np.allclose(mat, correct)


def test_implicit_swaps() -> None:
    c = Circuit(2)
    c.CX(0, 1)
    c.CX(1, 0)
    perm = c.implicit_qubit_permutation()
    assert all(a == b for (a, b) in perm.items())
    CliffordSimp().apply(c)
    perm1 = c.implicit_qubit_permutation()
    assert all(a != b for (a, b) in perm1.items())


def test_boxes() -> None:
    c = Circuit(4, 4)
    c.X(0)
    c.H(1)
    c.add_gate(OpType.ESWAP, 0.4, [0, 1])
    c.add_gate(OpType.FSim, [0.7, 0.1], [2, 3])
    c.S(1)
    c.CX(2, 0)
    c.CRz(0.5, 0, 3)
    cbox = CircBox(c)
    assert cbox.type == OpType.CircBox
    c1 = cbox.get_circuit()
    assert len(c1.get_commands()) == 7
    d = Circuit(4, 4)
    d.add_circbox(cbox, [0, 2, 1, 3, 3, 2, 1, 0])

    m = np.asarray([[1 / 2, sqrt(3) / 2], [sqrt(3) / 2, -1 / 2]])
    mbox = Unitary1qBox(m)
    m1 = mbox.get_matrix()
    assert np.allclose(m, m1)
    assert mbox.type == OpType.Unitary1qBox
    d.add_unitary1qbox(mbox, 3)

    u2q = np.asarray([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
    u2qbox = Unitary2qBox(u2q)
    u2q_ = u2qbox.get_matrix()
    assert np.allclose(u2q, u2q_)
    assert u2qbox.type == OpType.Unitary2qBox
    d.add_unitary2qbox(u2qbox, 1, 2)

    u3q = np.asarray(
        [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1],
        ]
    )
    u3qbox = Unitary3qBox(u3q)
    u3q_ = u3qbox.get_matrix()
    assert np.allclose(u3q, u3q_)
    assert u3qbox.type == OpType.Unitary3qBox
    d.add_unitary3qbox(u3qbox, 1, 2, 0)

    A = np.asarray(
        [[1, 2, 3, 4 + 1j], [2, 0, 1j, -1], [3, -1j, 2, 1j], [4 - 1j, -1, -1j, 1]]
    )
    ebox = ExpBox(A, 0.5)
    assert ebox.type == OpType.ExpBox
    d.add_expbox(ebox, 3, 2)

    paulis = [Pauli.X, Pauli.Z, Pauli.X]
    pbox = PauliExpBox(paulis, Symbol("alpha"))  # type: ignore
    assert pbox.type == OpType.PauliExpBox
    d.add_pauliexpbox(pbox, [3, 2, 1])

    qcbox = QControlBox(Op.create(OpType.S), 2)
    assert qcbox.type == OpType.QControlBox
    assert qcbox.get_op().type == OpType.S
    assert qcbox.get_n_controls() == 2
    d.add_qcontrolbox(qcbox, [1, 2, 3])

    assert d.n_gates == 7

    pauli_exps = [cmd.op for cmd in d if cmd.op.type == OpType.PauliExpBox]
    assert len(pauli_exps) == 1
    assert pauli_exps[0].get_paulis() == paulis
    assert pauli_exps[0].get_phase() == Symbol("alpha")  # type: ignore

    boxes = (cbox, mbox, u2qbox, u3qbox, ebox, pbox, qcbox)
    assert all(box == box for box in boxes)
    assert all(isinstance(box, Op) for box in boxes)


def test_u1q_stability() -> None:
    # https://github.com/CQCL/tket/issues/222
    u = np.array(
        [
            [
                -1.0000000000000000e00 + 0.0000000000000000e00j,
                -4.7624091282918654e-10 + 2.0295010872500105e-16j,
            ],
            [
                4.5447577055178555e-10 - 1.4232772405184710e-10j,
                -9.5429791447115209e-01 + 2.9885697320961047e-01j,
            ],
        ]
    )
    ubox = Unitary1qBox(u)
    op = ubox.get_circuit().get_commands()[0].op
    assert op.type == OpType.TK1
    a, b, c = op.params
    assert np.isfinite(a)
    assert np.isfinite(b)
    assert np.isfinite(c)


def test_custom_gates() -> None:
    a = Symbol("a")  # type: ignore
    b = Symbol("b")  # type: ignore
    setup = Circuit(3)
    setup.CX(0, 1)
    setup.Rz(a + 0.5, 2)
    setup.CRz(b, 0, 2)
    gatedef = CustomGateDef.define("g", setup, [a, b])
    c = Circuit(4)
    c.add_custom_gate(gatedef, [0.2, 1.3], [0, 3, 1])
    coms = c.get_commands()
    assert len(coms) == 1
    cmd0 = coms[0]
    assert str(cmd0) == "g(0.2,1.3) q[0], q[3], q[1];"
    gate = CustomGate(gatedef, [0.2, 1.3])
    op0 = cmd0.op
    assert gate.type == op0.type
    assert gate.params == op0.params
    Transform.DecomposeBoxes().apply(c)
    coms = c.get_commands()
    assert str(coms[0]) == "CX q[0], q[3];"
    assert str(coms[1]) == "Rz(0.7) q[1];"
    assert str(coms[2]) == "CRz(1.3) q[0], q[1];"


def test_errors() -> None:
    # TKET-289
    c = Circuit(1)
    a = Symbol("a")  # type: ignore
    c.Rz(a, 0)
    c.Rz(0.5, 0)
    c.Rz(0, 0)
    with pytest.raises(TypeError):
        c.Rz(0, "a")
    assert c.get_commands()[0].free_symbols() == set([a])


def test_str() -> None:
    # TKET-437
    c = Circuit(2).CRz(0.5, 0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CRz(0.5)"
    # TKET-871
    c = Circuit(2).CRx(0.5, 0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CRx(0.5)"
    c = Circuit(2).CRy(0.5, 0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CRy(0.5)"
    # TKET-957
    c = Circuit(2).CV(0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CV"
    c = Circuit(2).CVdg(0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CVdg"
    c = Circuit(2).CSX(0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CSX"
    c = Circuit(2).CSXdg(0, 1)
    op = c.get_commands()[0].op
    assert op.__str__() == "CSXdg"
    c = Circuit(1).SX(0)
    op = c.get_commands()[0].op
    assert op.__str__() == "SX"
    c = Circuit(1).SXdg(0)
    op = c.get_commands()[0].op
    assert op.__str__() == "SXdg"


def test_qubit_to_bit_map() -> None:
    c = Circuit()
    a = [Qubit("a", i) for i in range(4)]
    b = [Bit("b", i) for i in range(4)]
    for q in a:
        c.add_qubit(q)
    for cb in b:
        c.add_bit(cb)
    c.add_gate(OpType.X, [a[0]])
    c.H(a[1])
    c.add_gate(OpType.Rz, 0.5, [a[1]])
    c.CX(a[2], a[0])
    c.CRz(0.5, a[0], a[3])
    c.Measure(a[3], b[3])
    c.Measure(a[1], b[1])

    qbmap = c.qubit_to_bit_map
    assert len(qbmap) == 2
    assert qbmap[a[3]] == b[3]
    assert qbmap[a[1]] == b[1]


def test_barrier_errors() -> None:
    c = Circuit(2)
    c.add_barrier([0, 1])
    with pytest.raises(RuntimeError) as err:
        c.add_gate(OpType.Barrier, [], [0, 1])
    assert "Please use `add_barrier`" in str(err.value)


def test_phase() -> None:
    c = Circuit(2).H(0).CX(0, 1)
    c1 = c.copy()
    c1.add_phase(0.125)
    assert c1.phase - c.phase == 0.125
    c1.add_phase(1.875)
    assert c == c1


def test_phase_return_circ() -> None:
    c = Circuit(2).H(0).CX(0, 1)
    c1 = c.copy()
    c1.add_phase(0.125).H(0)
    assert c1.phase - c.phase == 0.125
    c1.add_phase(1.875)
    c.H(0)
    assert c == c1


@given(st.circuits())
def test_circuit_to_serializable_json_roundtrip(circuit: Circuit) -> None:
    serializable_form = circuit.to_dict()
    assert json.loads(json.dumps(serializable_form)) == serializable_form


@given(st.circuits())
def test_circuit_pickle_roundtrip(circuit: Circuit) -> None:
    assert pickle.loads(pickle.dumps(circuit)) == circuit


@given(st.circuits())
def test_circuit_from_to_serializable(circuit: Circuit) -> None:
    serializable_form = circuit.to_dict()
    validate(instance=serializable_form, schema=schema)
    assert circuit == Circuit.from_dict(serializable_form)


@given(st.circuits())
@settings(deadline=None)
def test_circuit_display(circuit: Circuit) -> None:
    html_str_circ = render_circuit_as_html(circuit, jupyter=False)
    html_str_dict = render_circuit_as_html(circuit.to_dict(), jupyter=False)
    assert isinstance(html_str_circ, str)
    assert isinstance(html_str_dict, str)


def test_circuit_display_with_barrier() -> None:
    # TKET-1434
    c = Circuit(2).X(0).X(1).add_barrier([0, 1]).X(0).X(1)
    html_str_circ = render_circuit_as_html(c, jupyter=False)
    html_str_dict = render_circuit_as_html(c.to_dict(), jupyter=False)
    assert isinstance(html_str_circ, str)
    assert isinstance(html_str_dict, str)


def test_ops_of_type() -> None:
    c = Circuit(2).H(0).Rz(0.5, 1).CX(0, 1).Rz(0.5, 1).Rz(-0.5, 0).CRy(0.3, 0, 1)
    ops_H = c.ops_of_type(OpType.H)
    ops_CX = c.ops_of_type(OpType.CX)
    ops_Rz = c.ops_of_type(OpType.Rz)
    ops_Rx = c.ops_of_type(OpType.Rx)
    ops_CRy = c.ops_of_type(OpType.CRy)
    assert len(ops_H) == 1
    assert ops_H[0].type == OpType.H
    assert len(ops_CX) == 1
    assert ops_CX[0].type == OpType.CX
    assert len(ops_Rz) == 3
    assert all(op.type == OpType.Rz for op in ops_Rz)
    assert len(ops_Rx) == 0
    assert len(ops_CRy) == 1
    assert ops_CRy[0].type == OpType.CRy


def test_commands_of_type() -> None:
    c = Circuit(2).H(0).Rz(0.5, 1).CX(0, 1).Rz(0.5, 1).CZ(0, 1).Rz(-0.5, 0)
    cmds_H = c.commands_of_type(OpType.H)
    cmds_CX = c.commands_of_type(OpType.CX)
    cmds_Rz = c.commands_of_type(OpType.Rz)
    cmds_Rx = c.commands_of_type(OpType.Rx)
    assert len(cmds_H) == 1
    assert cmds_H[0].op.type == OpType.H
    assert cmds_H[0].qubits == [c.qubits[0]]
    assert len(cmds_CX) == 1
    assert cmds_CX[0].op.type == OpType.CX
    assert cmds_CX[0].qubits == c.qubits
    assert len(cmds_Rz) == 3
    assert all(cmd.op.type == OpType.Rz for cmd in cmds_Rz)
    assert [cmd.qubits[0] for cmd in cmds_Rz] == [c.qubits[i] for i in [1, 1, 0]]
    assert len(cmds_Rx) == 0


def test_cempty_circuit() -> None:
    circ = Circuit(0)
    circt_dict = circ.to_dict()
    assert type(circt_dict) == type({})
    assert len(circt_dict) > 0
    assert Circuit(0) == Circuit(0)


def with_empty_qubit(op: Op) -> CircBox:
    n_qb = op.n_qubits
    return CircBox(Circuit(n_qb + 1).add_gate(op, list(range(1, n_qb + 1))))


def with_control_qubit(op: Op) -> QControlBox:
    return QControlBox(op, 1)


def test_opgroups() -> None:
    # Replace with op
    c = Circuit(1)
    c.Rz(0.2, 0, opgroup="Z rotations")
    c.H(0)
    c.Rz(0.3, 0, opgroup="Z rotations")
    newop = Op.create(OpType.T)
    assert c.substitute_named(newop, "Z rotations")

    # Replace with circuit
    newcirc = Circuit(1).H(0)
    assert c.substitute_named(newcirc, "Z rotations")

    # Add a control to a bunch of H and CX ops:
    c = Circuit(3)
    h_op = Op.create(OpType.H)
    cx_op = Op.create(OpType.CX)
    h_0_cbox = with_empty_qubit(h_op)
    h_q_qbox = with_control_qubit(h_op)
    cx_0_cbox = with_empty_qubit(cx_op)
    cx_q_qbox = with_control_qubit(cx_op)
    c.X(0).Y(1)
    c.add_circbox(h_0_cbox, [2, 0], opgroup="hgroup")
    c.add_circbox(cx_0_cbox, [2, 0, 1], opgroup="cxgroup")
    c.Y(0).X(1)
    c.add_circbox(h_0_cbox, [2, 1], opgroup="hgroup")
    c.add_circbox(cx_0_cbox, [2, 1, 0], opgroup="cxgroup")
    c.X(0).Y(1)
    assert c.substitute_named(h_q_qbox, "hgroup")
    assert c.substitute_named(cx_q_qbox, "cxgroup")

    # Add a pi/4 rotation to an Rz op:
    c = Circuit(2).H(0).H(1)
    c.Rz(0.5, 0, opgroup="rot")
    c.CX(0, 1)
    cbox = CircBox(Circuit(1).Rz(0.5, 0, opgroup="rot").Rz(0.25, 0, opgroup="rot"))
    assert c.substitute_named(cbox, "rot")
    assert DecomposeBoxes().apply(c)
    assert RemoveRedundancies().apply(c)
    assert c.n_gates_of_type(OpType.Rz) == 1
    # When all squashed vertices have the same name, the squashed vertex will have that
    # name too.
    # Swap the Rz for an Ry
    ry_op = Op.create(OpType.Ry, 0.75)
    assert c.substitute_named(ry_op, "rot")
    assert c.n_gates_of_type(OpType.Rz) == 0
    assert c.n_gates_of_type(OpType.Ry) == 1

    # Remove a redundant gate
    c = Circuit(3).H(0)
    assert len(c.opgroups) == 0
    c.CX(0, 1, opgroup="cx0")
    c.CX(1, 2, opgroup="cx1")
    c.CX(2, 0, opgroup="cx2")
    c.CX(0, 1, opgroup="cx3")
    assert c.opgroups == {"cx0", "cx1", "cx2", "cx3"}
    c.substitute_named(Circuit(2), "cx3")
    assert c.n_gates == 4
    assert c.n_gates_of_type(OpType.CX) == 3
    assert c.opgroups == {"cx0", "cx1", "cx2"}


def test_phase_polybox() -> None:
    c = Circuit(1, 1)
    n_qb = 1
    qubit_indices = {Qubit(0): 0}
    phase_polynomial = {(True,): 0.1}
    linear_transformation = np.array([[1]])
    p_box = PhasePolyBox(n_qb, qubit_indices, phase_polynomial, linear_transformation)

    b = p_box.get_circuit()

    p_box_ii = PhasePolyBox(b)

    c.add_phasepolybox(p_box, [0])
    c.add_phasepolybox(p_box, [Qubit(0)])

    c.add_phasepolybox(p_box_ii, [0])

    assert p_box.n_qubits == n_qb
    assert p_box_ii.n_qubits == n_qb
    assert p_box.qubit_indices == qubit_indices
    assert p_box_ii.qubit_indices == qubit_indices
    assert p_box.phase_polynomial == phase_polynomial
    assert p_box_ii.phase_polynomial == phase_polynomial
    assert np.array_equal(p_box.linear_transformation, linear_transformation)
    assert np.array_equal(p_box_ii.linear_transformation, linear_transformation)
    assert DecomposeBoxes().apply(c)


def test_phase_polybox_II() -> None:
    c = Circuit(1, 1)
    n_qb = 1
    qubit_indices = {Qubit(0): 0}
    phase_polynomial = {(True,): 0.1, (True,): 0.3}
    linear_transformation = np.array([[1]])
    p_box = PhasePolyBox(n_qb, qubit_indices, phase_polynomial, linear_transformation)

    b = p_box.get_circuit()

    p_box_ii = PhasePolyBox(b)

    c.add_phasepolybox(p_box, [0])
    c.add_phasepolybox(p_box, [Qubit(0)])

    c.add_phasepolybox(p_box_ii, [0])

    assert p_box.n_qubits == n_qb
    assert p_box_ii.n_qubits == n_qb
    assert p_box.qubit_indices == qubit_indices
    assert p_box_ii.qubit_indices == qubit_indices
    assert p_box.phase_polynomial == phase_polynomial
    assert p_box_ii.phase_polynomial == phase_polynomial
    assert np.array_equal(p_box.linear_transformation, linear_transformation)
    assert np.array_equal(p_box_ii.linear_transformation, linear_transformation)
    assert DecomposeBoxes().apply(c)


def test_phase_polybox_big() -> None:
    c = Circuit(3, 3)
    n_qb = 3
    qubit_indices = {Qubit(0): 0, Qubit(1): 1, Qubit(2): 2}
    phase_polynomial = {
        (True, False, True): 0.333,
        (False, False, True): 0.05,
        (False, True, False): 1.05,
    }
    linear_transformation = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    p_box = PhasePolyBox(n_qb, qubit_indices, phase_polynomial, linear_transformation)

    b = p_box.get_circuit()

    p_box_ii = PhasePolyBox(b)

    c.add_phasepolybox(p_box, [0, 1, 2])
    c.add_phasepolybox(p_box_ii, [0, 1, 2])
    assert p_box.n_qubits == n_qb
    assert p_box_ii.n_qubits == n_qb
    assert p_box.qubit_indices == qubit_indices
    assert p_box_ii.qubit_indices == qubit_indices
    assert p_box.phase_polynomial == phase_polynomial
    assert p_box_ii.phase_polynomial == phase_polynomial
    assert np.array_equal(p_box.linear_transformation, linear_transformation)
    assert np.array_equal(p_box_ii.linear_transformation, linear_transformation)
    assert DecomposeBoxes().apply(c)


def test_depth() -> None:
    c = Circuit(3)
    c.H(0).H(1).CX(1, 2).CZ(0, 1).H(1).CZ(1, 2)
    c.add_barrier([0, 1, 2])
    c.CX(1, 2).CZ(0, 1).H(1).CX(1, 2)
    assert c.depth() == 9
    assert c.depth_by_type(OpType.H) == 3
    assert c.depth_by_type(OpType.CX) == 3
    assert c.depth_by_type(OpType.CZ) == 3
    assert c.depth_by_type({OpType.CX, OpType.CZ}) == 6
    assert c.depth_by_type({OpType.CX, OpType.H}) == 6
    assert c.depth_by_type({OpType.CZ, OpType.H}) == 6
    assert c.depth_by_type(set()) == 0


def test_op_dagger_transpose() -> None:
    sx = Op.create(OpType.SX)
    sxdg = Op.create(OpType.SXdg)
    assert sx.dagger == sxdg
    assert sx.transpose == sx


def test_clifford_checking() -> None:
    c = Circuit(2, 1)
    c.H(0).CX(0, 1).T(1).Rz(0.5, 1).Rz(0.3, 1).Measure(1, 0)
    h = c.get_commands()[0].op
    assert h.is_clifford_type()
    cx = c.get_commands()[1].op
    assert cx.is_clifford_type()
    t = c.get_commands()[2].op
    assert t.is_clifford_type() == False
    rz1 = c.get_commands()[3].op
    assert rz1.is_clifford_type() == False
    rz2 = c.get_commands()[4].op
    assert rz2.is_clifford_type() == False
    m = c.get_commands()[5].op
    assert m.is_clifford_type() == False


def test_getting_registers() -> None:
    c = Circuit(2, 1)
    c_regs = c.c_registers
    assert len(c_regs) == 1
    assert c_regs[0] == BitRegister("c", 1)
    q_regs = c.q_registers
    assert len(q_regs) == 1
    assert q_regs[0] == QubitRegister("q", 2)
    q_err_msg = "Cannot find quantum register with name"
    c_err_msg = "Cannot find classical register with name"
    with pytest.raises(RuntimeError) as e:
        c.get_c_register("q")
    assert c_err_msg in str(e.value)
    with pytest.raises(RuntimeError) as e:
        c.get_q_register("c")
    assert q_err_msg in str(e.value)
    assert c.get_c_register("c").name == "c"
    assert c.get_c_register("c").size == 1
    assert c.get_q_register("q").name == "q"
    assert c.get_q_register("q").size == 2
    c.add_q_register("test_qr", 10)
    c.add_c_register("test_cr", 8)
    assert c.get_c_register("test_cr").name == "test_cr"
    assert c.get_c_register("test_cr").size == 8
    assert c.get_q_register("test_qr").name == "test_qr"
    assert c.get_q_register("test_qr").size == 10

    c_regs = c.c_registers
    c_regs.sort()
    assert len(c_regs) == 2
    assert c_regs[0] == BitRegister("c", 1)
    assert c_regs[1] == BitRegister("test_cr", 8)
    q_regs = c.q_registers
    q_regs.sort()
    assert len(q_regs) == 2
    assert q_regs[0] == QubitRegister("q", 2)
    assert q_regs[1] == QubitRegister("test_qr", 10)


def test_measuring_registers() -> None:
    c = Circuit()
    with pytest.raises(RuntimeError) as e:
        qreg = QubitRegister("qr", 2)
        c.measure_register(qreg, "cr")
    assert "The given QubitRegister is not in use" in str(e.value)
    qreg = c.add_q_register("qr", 2)
    c.measure_register(qreg, "cr")
    assert len(c.bits) == 2
    assert c.n_qubits == 2
    commands = c.get_commands()
    assert len(commands) == 2
    assert str(commands[0]) == "Measure qr[0] --> cr[0];"
    assert str(commands[1]) == "Measure qr[1] --> cr[1];"
    qreg2 = c.add_q_register("qr2", 2)
    with pytest.raises(RuntimeError) as e:
        c.measure_register(qreg2, "cr")
    assert 'A register with name "cr" already exists' in str(e.value)


def test_zzmax() -> None:
    c = Circuit(5)
    c.ZZMax(0, 1)
    assert c.depth() == 1


if __name__ == "__main__":
    test_circuit_gen()
    test_symbolic_ops()
    test_4x4_matrix_to_circ()
    test_exp_to_circ()
    test_boxes()
    test_errors()
    test_str()
    test_phase()
    test_clifford_checking()
    test_measuring_registers()
