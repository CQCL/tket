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

from enum import Enum, unique
from math import pi
from typing import NamedTuple

from lark import Lark, Transformer, Tree

from pytket.circuit import CircBox, Circuit, OpType

# The Lark grammar, transformer and type definitions below are adapted from the
# code in Eddie Schoute's `quippy` project
# (https://github.com/eddieschoute/quippy). The main enhancements are support
# for multi-qubit gates and correct handling of negative controls on qubit 0.


# Types
class _Wire(NamedTuple):
    i: int


class _ControlWire(NamedTuple):
    wire: _Wire
    negative: bool


class _Control(NamedTuple):
    controlled: list[_ControlWire]
    no_control: bool


@unique
class _TypeAssignment_Type(Enum):
    Qbit = 1
    Cbit = 2


class _TypeAssignment(NamedTuple):
    wire: _Wire
    type: _TypeAssignment_Type


class _Gate:
    pass


@unique
class _QGate_Op(Enum):
    Not = 1  # Pauli X
    H = 2  # Hadamard
    MultiNot = 3  # Multi-target not
    Y = 4  # Pauli Y
    Z = 5  # Pauli Z
    S = 6  # Clifford S
    T = 7  # Clifford T = sqrt(S)
    E = 8  # Clifford E = H S (omega)^3
    Omega = 9  # scalar = exp(i pi / 4)
    V = 10  # V = sqrt(X)
    Swap = 11
    W = 12  # W is self-inverse and diagonalizes Swap
    IX = 13  # i X


class _QGate(
    _Gate,
    NamedTuple(
        "QGate",
        [
            ("op", _QGate_Op),
            ("inverted", bool),
            ("wires", list[_Wire]),
            ("control", _Control),
        ],
    ),
):
    pass


@unique
class _QRot_Op(Enum):
    ExpZt = 1  # exp(−i Z t)  # noqa: RUF003
    R = 2  # R(2 / 2^t) (notation is confusing but see Monad.hs)


class _QRot(
    _Gate,
    NamedTuple(
        "QRot",
        [("op", _QRot_Op), ("inverted", bool), ("timestep", float), ("wire", _Wire)],
    ),
):
    pass


class _QInit(_Gate, NamedTuple("QInit", [("value", bool), ("wire", _Wire)])):
    pass


class _CInit(_Gate, NamedTuple("CInit", [("value", bool), ("wire", _Wire)])):
    pass


class _QTerm(_Gate, NamedTuple("QTerm", [("value", bool), ("wire", _Wire)])):
    pass


class _CTerm(_Gate, NamedTuple("CTerm", [("value", bool), ("wire", _Wire)])):
    pass


class _QMeas(_Gate, NamedTuple("QMeas", [("wire", _Wire)])):
    pass


class _QDiscard(_Gate, NamedTuple("QDiscard", [("wire", _Wire)])):
    pass


class _CDiscard(_Gate, NamedTuple("CDiscard", [("wire", _Wire)])):
    pass


class _SubroutineCall(
    _Gate,
    NamedTuple(
        "SubroutineCall",
        [
            ("repetitions", int),
            ("name", str),
            ("shape", str),
            ("inverted", bool),
            ("inputs", list[_Wire]),
            ("outputs", list[_Wire]),
            ("control", _Control),
        ],
    ),
):
    pass


class _Comment(
    _Gate,
    NamedTuple(
        "Comment",
        [
            ("comment", str),
            ("inverted", bool),
            ("wire_comments", list[tuple[_Wire, str]]),
        ],
    ),
):
    pass


class _Program(NamedTuple):
    inputs: list[_TypeAssignment]
    gates: list[_Gate]
    outputs: list[_TypeAssignment]


@unique
class _Subroutine_Control(Enum):
    yes = 1
    no = 2
    classically = 3


class _Subroutine(NamedTuple):
    name: str
    shape: str
    controllable: _Subroutine_Control
    circuit: _Program


class _Start(NamedTuple):
    circuit: _Program
    subroutines: list[_Subroutine]


# Transformer
class _QuipperTransformer(Transformer):
    def int(self, t: list) -> int:
        return int(t[0])

    def float(self, t: list) -> float:
        return float(t[0])

    def string(self, t: list) -> str:
        return str(t[0][1:-1])

    def wire(self, t: list) -> _Wire:
        return _Wire(t[0])

    wire_list = list

    def wire_string_list(self, t: list) -> list[tuple[_Wire, str]]:
        wires = (el for i, el in enumerate(t) if i % 2 == 0)
        labels = (el for i, el in enumerate(t) if i % 2 == 1)
        return list(zip(wires, labels, strict=False))

    def pos_control_wire(self, t: list) -> _ControlWire:
        return _ControlWire(t[0], False)

    control_wire_list = list

    def neg_control_wire(self, t: list) -> _ControlWire:
        return _ControlWire(t[0], True)

    def type_assignment(self, t: list) -> _TypeAssignment:
        ty = _TypeAssignment_Type.Qbit if t[1] == "Qbit" else _TypeAssignment_Type.Cbit
        return _TypeAssignment(t[0], ty)

    def arity(self, t: list) -> list[Tree]:
        return list(t)

    def qgate(self, t: list) -> _QGate:  # noqa: PLR0912
        ops = _QGate_Op
        n = t[0]
        if n in {"not", "x", "X"}:
            op = ops.Not
        elif n == "H":
            op = ops.H
        elif n == "multinot":
            op = ops.MultiNot
        elif n == "Y":
            op = ops.Y
        elif n == "Z":
            op = ops.Z
        elif n == "S":
            op = ops.S
        elif n == "E":
            op = ops.E
        elif n == "T":
            op = ops.T
        elif n == "V":
            op = ops.V
        elif n == "swap":
            op = ops.Swap
        elif n == "omega":
            op = ops.Omega
        elif n == "iX":
            op = ops.IX
        elif n == "W":
            op = ops.W
        else:
            raise RuntimeError(f"Unknown QGate operation: {n}")
        return _QGate(op=op, inverted=len(t[1].children) > 0, wires=t[2], control=t[3])

    def qrot1(self, t: list) -> _QRot:
        return _QRot(
            op=_QRot_Op.ExpZt, timestep=t[0], inverted=len(t[1].children) > 0, wire=t[2]
        )

    def qrot2(self, t: list) -> _QRot:
        return _QRot(
            op=_QRot_Op.R, timestep=t[0], inverted=len(t[1].children) > 0, wire=t[2]
        )

    def qinit(self, t: list) -> _QInit:
        return _QInit(value=(t[0] == "QInit1"), wire=t[1])

    def cinit(self, t: list) -> _CInit:
        return _CInit(value=(t[0] == "CInit1"), wire=t[1])

    def qterm(self, t: list) -> _QTerm:
        return _QTerm(value=(t[0] == "QTerm1"), wire=t[1])

    def cterm(self, t: list) -> _CTerm:
        return _CTerm(value=(t[0] == "CTerm1"), wire=t[1])

    def qmeas(self, t: list) -> _QMeas:
        return _QMeas(wire=t[0])

    def qdiscard(self, t: list) -> _QDiscard:
        return _QDiscard(wire=t[0])

    def cdiscard(self, t: list) -> _CDiscard:
        return _CDiscard(wire=t[0])

    def subroutine_call(self, t: list) -> _SubroutineCall:
        repetitions = 1
        if t[0] is not None:
            assert isinstance(t[0], int)
            repetitions = t[0]
        return _SubroutineCall(
            repetitions=repetitions,
            name=t[1],
            shape=t[2],
            inverted=len(t[3].children) > 0,
            inputs=t[4],
            outputs=t[5],
            control=t[6],
        )

    def comment(self, t: list) -> _Comment:
        wire_comments = t[2] if len(t) > 2 else []  # noqa: PLR2004
        return _Comment(
            comment=t[0], inverted=len(t[1].children) > 0, wire_comments=wire_comments
        )

    def control_app(self, t: list) -> _Control:
        if not t:
            return _Control(controlled=[], no_control=False)
        if len(t) == 2:  # noqa: PLR2004
            return _Control(controlled=t[0], no_control=True)
        if t[0] == "with nocontrol":
            return _Control(controlled=[], no_control=True)
        return _Control(controlled=t[0], no_control=False)

    def circuit(self, t: list) -> _Program:
        return _Program(inputs=t[0], gates=t[1:-1], outputs=t[-1])

    def subroutine(self, t: list) -> _Subroutine:
        if t[2] == "yes":
            controllable = _Subroutine_Control.yes
        elif t[2] == "no":
            controllable = _Subroutine_Control.no
        else:
            controllable = _Subroutine_Control.classically
        return _Subroutine(
            name=t[0], shape=t[1], controllable=controllable, circuit=t[3]
        )

    def start(self, t: list) -> _Start:
        circuit = t.pop(0)
        return _Start(circuit, list(t))


# Utility function
def _allowed(op: str, arity: int) -> bool:
    if op in ["Not", "IX", "H", "Y", "Z", "S", "T", "E", "Omega", "V"]:
        return arity == 1
    if op in ["Swap", "W"]:
        return arity == 2  # noqa: PLR2004
    # MultiNot
    return True


# Class for constructing a pytket Circuit from a parsed Quipper program
class _CircuitMaker:
    def __init__(self, subr: list[_Subroutine]) -> None:
        self.subrd = {s.name: s for s in subr}
        if len(self.subrd) != len(subr):
            raise TypeError("Repeated subroutine names")

    def make_circuit(self, circ: _Program) -> Circuit:  # noqa: PLR0912, PLR0915
        inps, outs, gates = circ.inputs, circ.outputs, circ.gates
        if inps != outs:
            raise TypeError("Inputs don't match outputs")
        qbits = [inp.wire.i for inp in inps if inp.type.name == "Qbit"]
        cbits = [inp.wire.i for inp in inps if inp.type.name == "Cbit"]
        n_qbits = len(qbits)
        n_cbits = len(cbits)
        # Construct mappings from wire labels to tket indices.
        tkqbits = {qbits[i]: i for i in range(n_qbits)}
        tkcbits = {cbits[i]: i for i in range(n_cbits)}
        # Construct circuit in tket.
        c = Circuit(n_qbits, n_cbits)
        for gate in gates:
            if isinstance(gate, _SubroutineCall):
                s = self.subrd[gate.name]
                if gate.shape != s.shape:
                    # Just a sanity check, otherwise we assume 'shape' OK.
                    raise TypeError("Mismatched shape in subroutine call")
                if gate.inputs != gate.outputs:
                    raise TypeError("Mismatched outputs in subroutine call")
                if gate.control.controlled:
                    raise NotImplementedError("Controls on subroutine")
                subcirc = self.make_circuit(s.circuit)  # recurse
                if gate.inverted:
                    subcirc = subcirc.dagger()
                cbox = CircBox(subcirc)
                for _ in range(gate.repetitions):
                    c.add_circbox(cbox, [wire.i for wire in gate.inputs])
            elif isinstance(gate, _QGate):
                # The `nocontrol' flag is irrelevant here.
                op = gate.op.name
                qctrls, neg_qctrls, cctrls = [], [], []
                for w in gate.control.controlled:
                    idx = w.children[0].wire.i  # type: ignore
                    is_neg = w.children[0].negative  # type: ignore
                    if idx in qbits:
                        tkqidx = tkqbits[idx]
                        qctrls.append(tkqidx)
                        if is_neg:
                            neg_qctrls.append(tkqidx)
                    elif idx in cbits:
                        tkcidx = tkcbits[idx]
                        cctrls.append(tkcidx)
                if cctrls:
                    raise NotImplementedError("Classical controls")
                inv = gate.inverted
                wires = [tkqbits[wire.i] for wire in gate.wires]  # all must be qubits
                n_wires = len(wires)
                if not _allowed(op, n_wires):
                    raise TypeError("'%s' gate with %d wires" % (op, n_wires))  # noqa: UP031
                n_ctrls = len(qctrls)
                # Negative control values must be handled using NOT gates
                # either side.
                for ctrl in neg_qctrls:
                    c.X(ctrl)
                if op == "Not":
                    if n_ctrls == 0:
                        c.X(wires[0])
                    elif n_ctrls == 1:
                        c.CX(qctrls[0], wires[0])
                    elif n_ctrls == 2:  # noqa: PLR2004
                        c.CCX(qctrls[0], qctrls[1], wires[0])
                    else:
                        c.add_gate(OpType.CnX, qctrls + wires)
                elif op == "IX":
                    if n_ctrls == 0:
                        c.X(wires[0])  # ignore phase
                    elif n_ctrls == 1:
                        c.S(qctrls[0])
                        c.CX(qctrls[0], wires[0])
                    else:
                        raise NotImplementedError("IX with more than 1 control")
                elif op == "MultiNot":
                    # Implement as a series of (controlled) X gates.
                    for wire in wires:
                        if n_ctrls == 0:
                            c.X(wire)
                        elif n_ctrls == 1:
                            c.CX(qctrls[0], wire)
                        elif n_ctrls == 2:  # noqa: PLR2004
                            c.CCX(qctrls[0], qctrls[1], wire)
                        else:
                            c.add_gate(OpType.CnX, [*qctrls, wire])
                elif op == "H":
                    if n_ctrls == 0:
                        c.H(wires[0])
                    elif n_ctrls == 1:
                        c.CH(qctrls[0], wires[0])
                    else:
                        raise NotImplementedError("H with more than 1 control")
                elif op == "Y":
                    if n_ctrls == 0:
                        c.Y(wires[0])
                    elif n_ctrls == 1:
                        c.CY(qctrls[0], wires[0])
                    else:
                        raise NotImplementedError("Y with more than 1 control")
                elif op == "Z":
                    if n_ctrls == 0:
                        c.Z(wires[0])
                    elif n_ctrls == 1:
                        c.CZ(qctrls[0], wires[0])
                    else:
                        raise NotImplementedError("Z with more than 1 control")
                elif op == "S":
                    if n_ctrls == 0:
                        if inv:
                            c.Sdg(wires[0])
                        else:
                            c.S(wires[0])
                    elif n_ctrls == 1:
                        # S = U1(1/2)
                        if inv:
                            c.add_gate(OpType.CU1, -0.5, [qctrls[0], wires[0]])
                        else:
                            c.add_gate(OpType.CU1, 0.5, [qctrls[0], wires[0]])
                    else:
                        raise NotImplementedError("S with more than 1 control")
                elif op == "T":
                    if n_ctrls == 0:
                        if inv:
                            c.Tdg(wires[0])
                        else:
                            c.T(wires[0])
                    elif n_ctrls == 1:
                        # T = U1(1/4)
                        if inv:
                            c.add_gate(OpType.CU1, -0.25, [qctrls[0], wires[0]])
                        else:
                            c.add_gate(OpType.CU1, 0.25, [qctrls[0], wires[0]])
                    else:
                        raise NotImplementedError("T with more than 1 control")
                elif op == "E":
                    # E = H S^3 omega^3 = H Sdg up to a scalar.
                    if n_ctrls == 0:
                        if inv:
                            c.H(wires[0])
                            c.S(wires[0])
                        else:
                            c.Sdg(wires[0])
                            c.H(wires[0])
                    else:
                        raise NotImplementedError("Controlled E")
                elif op == "Omega":
                    if n_ctrls == 0:
                        pass  # ignore phase
                    elif n_ctrls == 1:
                        c.Rz(0.25, qctrls[0])
                    else:
                        raise NotImplementedError("Omega with more than 1 control")
                elif op == "V":
                    if n_ctrls == 0:
                        if inv:
                            c.add_gate(OpType.Vdg, wires)
                        else:
                            c.add_gate(OpType.V, wires)
                    else:
                        raise NotImplementedError("Controlled V")
                elif op == "Swap":
                    if n_ctrls == 0:
                        c.SWAP(wires[0], wires[1])
                    elif n_ctrls == 1:
                        c.CSWAP(qctrls[0], wires[0], wires[1])
                    else:
                        raise NotImplementedError("SWAP with more than 2 controls")
                elif op == "W":
                    # W is self-inverse.
                    # Is there a simpler way to synthesize this?
                    if n_ctrls == 0:
                        c.Rz(1.25, wires[0])
                        c.Ry(1.0, wires[0])
                        c.T(wires[1])
                        c.Ry(1.0, wires[1])
                        c.CX(wires[1], wires[0])
                        c.Ry(-0.25, wires[1])
                        c.CX(wires[0], wires[1])
                        c.Ry(0.25, wires[1])
                        c.CX(wires[1], wires[0])
                        c.T(wires[0])
                        c.Ry(1.0, wires[0])
                        c.Rz(1.25, wires[1])
                        c.Ry(1.0, wires[1])
                    else:
                        raise NotImplementedError("Controlled W")
                else:
                    raise TypeError(f"Unknown op type: {op}")
                # Apply the NOT gates again for the negative controls.
                for ctrl in neg_qctrls:
                    c.X(ctrl)
            elif isinstance(gate, _QRot):
                op = gate.op.name
                inv = gate.inverted
                t = gate.timestep
                wire = tkqbits[gate.wire.i]  # must be quantum
                if op == "ExpZt":
                    if inv:
                        c.Rz(-2 * t / pi, wire)
                    else:
                        c.Rz(2 * t / pi, wire)
                elif op == "R":
                    if inv:
                        c.Rz(-2 / t, wire)
                    else:
                        c.Rz(2 / t, wire)
                else:
                    raise TypeError(f"Unknown op type: {op}")
            # QInit, QTerm, CInit, CTerm represent 'temporary' wires that
            # only occupy part of a circuit (and can be initialized with 0/1).
            # Not supported in pytket.
            elif isinstance(gate, _QInit):
                wire = tkqbits[gate.wire.i]  # must be quantum
                raise NotImplementedError("QInit")
            elif isinstance(gate, _CInit):
                wire = tkcbits[gate.wire.i]  # must be classical
                raise NotImplementedError("CInit")
            elif isinstance(gate, _QTerm):
                wire = tkqbits[gate.wire.i]  # must be quantum
                raise NotImplementedError("QTerm")
            elif isinstance(gate, _CTerm):
                wire = tkcbits[gate.wire.i]  # must be classical
                raise NotImplementedError("CTerm")
            elif isinstance(gate, _QMeas):
                # QMeas turns a qubit into a bit.
                # Not supported in pytket.
                wire = tkqbits[gate.wire.i]  # must be quantum
                raise NotImplementedError("QMeas")
            elif isinstance(gate, _QDiscard):
                # QDiscard discards a qubit.
                # Not supported in pytket.
                wire = tkqbits[gate.wire.i]  # must be quantum
                raise NotImplementedError("QDiscard")
            elif isinstance(gate, _CDiscard):
                # CDiscard discards a bit.
                # Not supported in pytket.
                wire = tkcbits[gate.wire.i]  # must be classical
                raise NotImplementedError("CDiscard")
            elif isinstance(gate, _Comment):
                pass
            else:
                raise TypeError(f"Unknown gate type: {type(gate)}")
        return c


def circuit_from_quipper(input_file: str) -> Circuit:
    # pylint: disable=C0301
    """Generate a pytket :py:class:`~.Circuit` given a program in Quipper ASCII format.
    Limitations (due to current limitations of `pytket`):

    * Subroutines must be defined over the full set of qubits.

    * Global phases are ignored.

    * Only limited support for controlled gates (depending on the gate type).

    * No support for 'QInit' and 'QTerm'. All qubits must run from the beginning
      to the end of the circuit.

    * No support for the legacy keywords ('QNot', 'QMultinot', 'QHad', 'QSwap',
      'QW'). These are now represented in Quipper as types of 'QGate'.

    * No support for 'QMeas' (which in Quipper converts a Q wire to a C wire).

    * No support for: 'GPhase' (global phase), 'QPrep', 'QUnprep', 'QDiscard',
      or 'DTerm'.

    * No support for classical operations ('CNot', 'CGate', 'CSwap', 'CInit',
      'CTerm', 'CDiscard').

    :param input_file: name of file to read
    """

    # Read Quipper program from file.
    with open(input_file) as f:
        quip = f.read()

    # Parse the circuit using the _QuipperTransformer.
    x = Lark(
        """
        start : circuit subroutine* _NEWLINE*
        circuit : "Inputs:" arity (gate _NEWLINE)* "Outputs:" arity
        subroutine: _NEWLINE "Subroutine:" string _NEWLINE "Shape:" string _NEWLINE "Controllable:" SUB_CONTROL _NEWLINE circuit
        SUB_CONTROL : "yes"
            | "no"
            | "classically"
        arity : type_assignment ("," type_assignment)* ","? _NEWLINE
        type_assignment : wire ":" TYPE
        TYPE: "Qbit"
            | "Cbit"
        control_app : controlled? NO_CONTROL?
        ?controlled : "with controls=[" control_wire_list "]"
        NO_CONTROL : "with nocontrol"
        ?gate : qgate
            | qrot1
            | qrot2
            | gphase
            | cnot
            | cgate
            | cswap
            | qprep
            | qunprep
            | qinit
            | cinit
            | qterm
            | cterm
            | qmeas
            | qdiscard
            | cdiscard
            | dterm
            | subroutine_call
            | comment
        !inversion  : "*"?
        qgate       : "QGate[" string "]" inversion "(" wire_list ")" control_app
        qrot1       : "QRot[\\"exp(-i%Z)\\"," float "]" inversion "(" wire ")"
        qrot2       : "QRot[\\"R(2pi/%)\\"," int "]" inversion "(" wire ")"
        gphase      : "Gphase() with t=" float control_app "with anchors=[" wire_list "]"
        cnot        : "CNot(" wire ")" control_app
        cgate       : "CGate[" string "]" inversion "(" wire_list ")" NO_CONTROL?
        cswap       : "CSwap(" wire "," wire ")" control_app
        qprep       : "QPrep(" wire ")" NO_CONTROL?
        qunprep     : "QUnprep(" wire ")" NO_CONTROL?
        qinit       : QINIT_STATE "(" wire ")" NO_CONTROL?
        QINIT_STATE : "QInit0" | "QInit1"
        cinit       :  CINIT_STATE "(" wire ")" NO_CONTROL?
        CINIT_STATE : "CInit0" | "CInit1"
        qterm       : QTERM_STATE "(" wire ")" NO_CONTROL?
        QTERM_STATE : "QTerm0" | "QTerm1"
        cterm       : CTERM_STATE "(" wire ")" NO_CONTROL?
        CTERM_STATE : "CTerm0" | "CTerm1"
        qmeas       : "QMeas(" wire ")"
        qdiscard    : "QDiscard(" wire ")"
        cdiscard    : "CDiscard(" wire ")"
        dterm       : DTERM_STATE "(" wire ")"
        DTERM_STATE : "DTerm0" | "DTerm1"
        subroutine_call : "Subroutine" ["(x" int ")"] "[" string ", shape" string "]" inversion "(" wire_list ") -> (" wire_list ")" control_app
        comment : "Comment[" string "]" inversion "(" wire_string_list? ")"
        wire_string_list: wire ":" string ("," wire ":" string)*
        wire : int
        wire_list : wire ("," wire)*
        pos_control_wire : "+" wire
        neg_control_wire : "-" wire
        control_wire : pos_control_wire | neg_control_wire
        control_wire_list : control_wire ("," control_wire)*
        %import common.WS_INLINE -> WS
        %ignore WS
        %import common.CR
        %import common.LF
        _NEWLINE: CR? LF
        %import common.ESCAPED_STRING
        string: ESCAPED_STRING
        %import common.SIGNED_FLOAT
        float : SIGNED_FLOAT
        %import common.INT
        int : INT
        """,
        parser="lalr",
        transformer=_QuipperTransformer(),
    ).parse(quip)

    # Load the subroutine list.
    maker = _CircuitMaker(x.subroutines)  # type: ignore

    # Make the tket circuit.
    return maker.make_circuit(x.circuit)  # type: ignore
