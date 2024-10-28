# Copyright 2019-2024 Quantinuum
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


"""Methods to allow conversion between Qiskit and pytket circuit classes"""

from typing import (
    Optional,
    Any,
    cast,
    Callable,
    List,
    Tuple,
)
from inspect import signature
from uuid import UUID

import numpy as np
from numpy.typing import NDArray
from symengine import sympify  # type: ignore
from symengine.lib import symengine_wrapper  # type: ignore

import sympy
import qiskit.circuit.library.standard_gates as qiskit_gates  # type: ignore
from qiskit.circuit.quantumcircuitdata import QuantumCircuitData  # type: ignore

from qiskit import (
    ClassicalRegister,
    QuantumCircuit,
    QuantumRegister,
)
from qiskit.circuit import (
    Barrier,
    Instruction,
    InstructionSet,
    Gate,
    ControlledGate,
    Measure,
    Parameter,
    ParameterExpression,
    Reset,
    Clbit,
)
from qiskit.circuit.library import (
    CRYGate,
    RYGate,
    PauliEvolutionGate,
    StatePreparation,
    UnitaryGate,
    Initialize,
)


from qiskit.transpiler import PassManager, CouplingMap  # type: ignore
from qiskit.transpiler.preset_passmanagers.builtin_plugins import SabreLayoutPassManager  # type: ignore
from qiskit.transpiler.passmanager_config import PassManagerConfig  # type: ignore

from pytket.circuit import (
    CircBox,
    Circuit,
    Node,
    Op,
    OpType,
    Unitary1qBox,
    Unitary2qBox,
    Unitary3qBox,
    UnitType,
    UnitID,
    Bit,
    Qubit,
    QControlBox,
    StatePreparationBox,
)


from pytket.unit_id import _TEMP_BIT_NAME
from pytket.pauli import Pauli, QubitPauliString
from pytket.architecture import Architecture
from pytket.utils import (
    QubitPauliOperator,
    gen_term_sequence_circuit,
    permute_rows_cols_in_unitary,
)
from pytket.passes import AutoRebase

_qiskit_gates_1q = {
    # Exact equivalents (same signature except for factor of pi in each parameter):
    qiskit_gates.HGate: OpType.H,
    qiskit_gates.IGate: OpType.noop,
    qiskit_gates.PhaseGate: OpType.U1,
    qiskit_gates.RGate: OpType.PhasedX,
    qiskit_gates.RXGate: OpType.Rx,
    qiskit_gates.RYGate: OpType.Ry,
    qiskit_gates.RZGate: OpType.Rz,
    qiskit_gates.SdgGate: OpType.Sdg,
    qiskit_gates.SGate: OpType.S,
    qiskit_gates.SXdgGate: OpType.SXdg,
    qiskit_gates.SXGate: OpType.SX,
    qiskit_gates.TdgGate: OpType.Tdg,
    qiskit_gates.TGate: OpType.T,
    qiskit_gates.U1Gate: OpType.U1,
    qiskit_gates.U2Gate: OpType.U2,
    qiskit_gates.U3Gate: OpType.U3,
    qiskit_gates.UGate: OpType.U3,
    qiskit_gates.XGate: OpType.X,
    qiskit_gates.YGate: OpType.Y,
    qiskit_gates.ZGate: OpType.Z,
}

_qiskit_gates_2q = {
    # Exact equivalents (same signature except for factor of pi in each parameter):
    qiskit_gates.CHGate: OpType.CH,
    qiskit_gates.CPhaseGate: OpType.CU1,
    qiskit_gates.CRXGate: OpType.CRx,
    qiskit_gates.CRYGate: OpType.CRy,
    qiskit_gates.CRZGate: OpType.CRz,
    qiskit_gates.CUGate: OpType.CU3,
    qiskit_gates.CU1Gate: OpType.CU1,
    qiskit_gates.CU3Gate: OpType.CU3,
    qiskit_gates.CXGate: OpType.CX,
    qiskit_gates.CSXGate: OpType.CSX,
    qiskit_gates.CYGate: OpType.CY,
    qiskit_gates.CZGate: OpType.CZ,
    qiskit_gates.ECRGate: OpType.ECR,
    qiskit_gates.iSwapGate: OpType.ISWAPMax,
    qiskit_gates.RXXGate: OpType.XXPhase,
    qiskit_gates.RYYGate: OpType.YYPhase,
    qiskit_gates.RZZGate: OpType.ZZPhase,
    qiskit_gates.SwapGate: OpType.SWAP,
}

_qiskit_gates_other = {
    # Exact equivalents (same signature except for factor of pi in each parameter):
    qiskit_gates.C3XGate: OpType.CnX,
    qiskit_gates.C4XGate: OpType.CnX,
    qiskit_gates.CCXGate: OpType.CCX,
    qiskit_gates.CCZGate: OpType.CnZ,
    qiskit_gates.CSwapGate: OpType.CSWAP,
    # Multi-controlled gates (qiskit expects a list of controls followed by the target):
    qiskit_gates.MCXGate: OpType.CnX,
    qiskit_gates.MCXGrayCode: OpType.CnX,
    qiskit_gates.MCXRecursive: OpType.CnX,
    qiskit_gates.MCXVChain: OpType.CnX,
    # Special types:
    Barrier: OpType.Barrier,
    Instruction: OpType.CircBox,
    Gate: OpType.CircBox,
    Measure: OpType.Measure,
    Reset: OpType.Reset,
    Initialize: OpType.StatePreparationBox,
    StatePreparation: OpType.StatePreparationBox,
}

_known_qiskit_gate = {**_qiskit_gates_1q, **_qiskit_gates_2q, **_qiskit_gates_other}

# Some qiskit gates are aliases (e.g. UGate and U3Gate).
# In such cases this reversal will select one or the other.
_known_qiskit_gate_rev = {v: k for k, v in _known_qiskit_gate.items()}

# Ensure U3 maps to UGate. (U3Gate deprecated in Qiskit but equivalent.)
_known_qiskit_gate_rev[OpType.U3] = qiskit_gates.UGate

# There is a bijective mapping, but requires some special parameter conversions
# tk1(a, b, c) = U(b, a-1/2, c+1/2) + phase(-(a+c)/2)
_known_qiskit_gate_rev[OpType.TK1] = qiskit_gates.UGate

# some gates are only equal up to global phase, support their conversion
# from tket -> qiskit
_known_gate_rev_phase = {
    optype: (qgate, 0.0) for optype, qgate in _known_qiskit_gate_rev.items()
}

_known_gate_rev_phase[OpType.V] = (qiskit_gates.SXGate, -0.25)
_known_gate_rev_phase[OpType.Vdg] = (qiskit_gates.SXdgGate, 0.25)

# use minor signature hacks to figure out the string names of qiskit Gate objects
_gate_str_2_optype: dict[str, OpType] = dict()
for gate, optype in _known_qiskit_gate.items():
    if gate in (
        UnitaryGate,
        Instruction,
        Gate,
        qiskit_gates.MCXGate,  # all of these have special (c*n)x names
        qiskit_gates.MCXGrayCode,
        qiskit_gates.MCXRecursive,
        qiskit_gates.MCXVChain,
    ):
        continue
    sig = signature(gate.__init__)
    # name is only a property of the instance, not the class
    # so initialize with the correct number of dummy variables
    n_params = len([p for p in sig.parameters.values() if p.default is p.empty]) - 1
    name = gate(*([1] * n_params)).name
    _gate_str_2_optype[name] = optype

_gate_str_2_optype_rev = {v: k for k, v in _gate_str_2_optype.items()}
# the aliasing of the name is ok in the reverse map
_gate_str_2_optype_rev[OpType.Unitary1qBox] = "unitary"


def _qpo_from_peg(peg: PauliEvolutionGate, qubits: list[Qubit]) -> QubitPauliOperator:
    op = peg.operator
    t = peg.params[0]
    qpodict = {}
    for p, c in zip(op.paulis, op.coeffs):
        if np.iscomplex(c):
            raise ValueError("Coefficient for Pauli {} is non-real.".format(p))
        coeff = param_to_tk(t) * c
        qpslist = []
        pstr = p.to_label()
        for a in pstr:
            if a == "X":
                qpslist.append(Pauli.X)
            elif a == "Y":
                qpslist.append(Pauli.Y)
            elif a == "Z":
                qpslist.append(Pauli.Z)
            else:
                assert a == "I"
                qpslist.append(Pauli.I)
        qpodict[QubitPauliString(qubits, qpslist)] = coeff
    return QubitPauliOperator(qpodict)


def _string_to_circuit(
    circuit_string: str,
    n_qubits: int,
    qiskit_prep: Initialize | StatePreparation,
) -> Circuit:
    """Helper function to generate circuits for Initialize
    and StatePreparation objects built with strings"""

    circ = Circuit(n_qubits)
    # Check if Instruction is Initialize or Statepreparation
    # If Initialize, add resets
    if isinstance(qiskit_prep, Initialize):
        for qubit in circ.qubits:
            circ.Reset(qubit)

    # We iterate through the string in reverse to add the
    # gates in the correct order (endian-ness).
    for qubit_index, character in enumerate(reversed(circuit_string)):
        match character:
            case "0":
                pass
            case "1":
                circ.X(qubit_index)
            case "+":
                circ.H(qubit_index)
            case "-":
                circ.X(qubit_index)
                circ.H(qubit_index)
            case "r":
                circ.H(qubit_index)
                circ.S(qubit_index)
            case "l":
                circ.H(qubit_index)
                circ.Sdg(qubit_index)
            case _:
                raise ValueError(
                    f"Cannot parse string for character {character}. "
                    + "The supported characters are {'0', '1', '+', '-', 'r', 'l'}."
                )

    return circ


def _get_pytket_ctrl_state(bitstring: str, n_bits: int) -> tuple[bool, ...]:
    "Converts a little endian string '001'=1 (LE) to (1, 0, 0)."
    assert set(bitstring).issubset({"0", "1"})
    padded_bitstring = bitstring.zfill(n_bits)
    pytket_ctrl_state = reversed([bool(int(b)) for b in padded_bitstring])
    return tuple(pytket_ctrl_state)


def _all_bits_set(integer: int, n_bits: int) -> bool:
    return integer.bit_count() == n_bits


def _get_controlled_tket_optype(c_gate: ControlledGate) -> OpType:
    """Get a pytket contolled OpType from a qiskit ControlledGate."""

    # If the control state is not "all |1>", use QControlBox
    if not _all_bits_set(c_gate.ctrl_state, c_gate.num_ctrl_qubits):
        return OpType.QControlBox

    elif c_gate.base_class in _known_qiskit_gate:
        # First we check if the gate is in _known_qiskit_gate
        # this avoids CZ being converted to CnZ
        return _known_qiskit_gate[c_gate.base_class]

    match c_gate.base_gate.base_class:
        case qiskit_gates.RYGate:
            return OpType.CnRy
        case qiskit_gates.YGate:
            return OpType.CnY
        case qiskit_gates.ZGate:
            return OpType.CnZ
        case _:
            if (
                c_gate.base_gate.base_class in _known_qiskit_gate
                or c_gate.base_gate.base_class is UnitaryGate
            ):
                return OpType.QControlBox
            else:
                raise NotImplementedError(
                    "Conversion of qiskit ControlledGate with base gate "
                    + f"base gate {c_gate.base_gate}"
                    + "not implemented."
                )


def _optype_from_qiskit_instruction(instruction: Instruction) -> OpType:
    """Get a pytket OpType from a qiskit Instruction."""
    if isinstance(instruction, ControlledGate):
        return _get_controlled_tket_optype(instruction)
    try:
        optype = _known_qiskit_gate[instruction.base_class]
        return optype
    except KeyError:
        raise NotImplementedError(
            f"Conversion of qiskit's {instruction.name} instruction is "
            + "currently unsupported by qiskit_to_tk. Consider "
            + "using QuantumCircuit.decompose() before attempting "
            + "conversion."
        )


UnitaryBox = Unitary1qBox | Unitary2qBox | Unitary3qBox


def _get_unitary_box(unitary: NDArray[np.complex128], num_qubits: int) -> UnitaryBox:
    match num_qubits:
        case 1:
            assert unitary.shape == (2, 2)
            return Unitary1qBox(unitary)
        case 2:
            assert unitary.shape == (4, 4)
            return Unitary2qBox(unitary)
        case 3:
            assert unitary.shape == (8, 8)
            return Unitary3qBox(unitary)
        case _:
            raise NotImplementedError(
                f"Conversion of {num_qubits}-qubit unitary gates not supported."
            )


def _get_qcontrol_box(c_gate: ControlledGate, params: list[float]) -> QControlBox:
    qiskit_ctrl_state: str = bin(c_gate.ctrl_state)[2:]
    pytket_ctrl_state: tuple[bool, ...] = _get_pytket_ctrl_state(
        bitstring=qiskit_ctrl_state, n_bits=c_gate.num_ctrl_qubits
    )
    if isinstance(c_gate.base_gate, UnitaryGate):
        unitary = c_gate.base_gate.params[0]
        # Here we reverse the order of the columns to correct for endianness.
        new_unitary: NDArray[np.complex128] = permute_rows_cols_in_unitary(
            matrix=unitary,
            permutation=tuple(reversed(range(c_gate.base_gate.num_qubits))),
        )
        base_op: Op = _get_unitary_box(new_unitary, c_gate.base_gate.num_qubits)
    else:
        base_tket_gate: OpType = _known_qiskit_gate[c_gate.base_gate.base_class]

        base_op: Op = Op.create(base_tket_gate, params)  # type: ignore

    return QControlBox(
        base_op, n_controls=c_gate.num_ctrl_qubits, control_state=pytket_ctrl_state
    )


def _add_state_preparation(
    tkc: Circuit, qubits: list[Qubit], prep: Initialize | StatePreparation
) -> None:
    """Handles different cases of Initialize and StatePreparation
    and appends the appropriate state preparation to a Circuit instance."""

    # Check how Initialize or StatePrep is constructed
    # With a string, an int or an array of amplitudes
    if len(prep.params) != 1:
        if isinstance(prep.params[0], str):
            # Parse string to get the right single qubit gates
            circuit_string: str = "".join(prep.params)
            circuit = _string_to_circuit(
                circuit_string, prep.num_qubits, qiskit_prep=prep
            )
            tkc.add_circuit(circuit, qubits)
        else:
            amplitude_array: NDArray[np.complex128] = np.array(prep.params)
            pytket_state_prep_box = StatePreparationBox(
                amplitude_array, with_initial_reset=(type(prep) is Initialize)
            )

            # Need to reverse qubits here (endian-ness)
            reversed_qubits = list(reversed(qubits))
            tkc.add_gate(pytket_state_prep_box, reversed_qubits)
    elif isinstance(prep.params[0], complex):
        # convert int to a binary string and apply X for |1>
        integer_parameter = int(prep.params[0].real)
        bit_string = bin(integer_parameter)[2:]
        circuit = _string_to_circuit(bit_string, prep.num_qubits, qiskit_prep=prep)
        tkc.add_circuit(circuit, qubits)
    else:
        raise TypeError(
            "Unrecognised type of Instruction.params "
            + "when trying to convert Initialize or StatePreparation instruction."
        )


def _get_pytket_condition_kwargs(
    instruction: Instruction,
    cregmap: dict[str, ClassicalRegister],
    circuit: QuantumCircuit,
) -> dict[str, Any]:
    if type(instruction.condition[0]) is ClassicalRegister:
        cond_reg = cregmap[instruction.condition[0]]
        condition_kwargs = {
            "condition_bits": [cond_reg[k] for k in range(len(cond_reg))],
            "condition_value": instruction.condition[1],
        }
        return condition_kwargs
    elif type(instruction.condition[0]) is Clbit:
        # .find_bit() returns type:
        #    tuple[index, list[tuple[ClassicalRegister, index]]]
        # We assume each bit belongs to exactly one register.
        index = circuit.find_bit(instruction.condition[0])[0]
        register = circuit.find_bit(instruction.condition[0])[1][0][0]
        cond_reg = cregmap[register]
        condition_kwargs = {
            "condition_bits": [cond_reg[index]],
            "condition_value": instruction.condition[1],
        }
        return condition_kwargs
    else:
        raise NotImplementedError("condition must contain classical bit or register")


class CircuitBuilder:
    def __init__(
        self,
        qregs: list[QuantumRegister],
        cregs: Optional[list[ClassicalRegister]] = None,
        name: Optional[str] = None,
        phase: Optional[sympy.Expr] = None,
    ):
        self.qregs = qregs
        self.cregs = [] if cregs is None else cregs
        self.qbmap = {}
        self.cbmap = {}
        if name is not None:
            self.tkc = Circuit(name=name)
        else:
            self.tkc = Circuit()
        if phase is not None:
            self.tkc.add_phase(phase)
        for reg in qregs:
            self.tkc.add_q_register(reg.name, len(reg))
            for i, qb in enumerate(reg):
                self.qbmap[qb] = Qubit(reg.name, i)
        self.cregmap = {}
        for reg in self.cregs:
            tk_reg = self.tkc.add_c_register(reg.name, len(reg))
            self.cregmap.update({reg: tk_reg})
            for i, cb in enumerate(reg):
                self.cbmap[cb] = Bit(reg.name, i)

    def circuit(self) -> Circuit:
        return self.tkc

    def add_qiskit_data(
        self, circuit: QuantumCircuit, data: Optional["QuantumCircuitData"] = None
    ) -> None:
        data = data or circuit.data
        for datum in data:
            instr, qargs, cargs = datum.operation, datum.qubits, datum.clbits

            qubits: list[Qubit] = [self.qbmap[qbit] for qbit in qargs]
            bits: list[Bit] = [self.cbmap[bit] for bit in cargs]

            condition_kwargs = {}
            if instr.condition is not None:
                condition_kwargs = _get_pytket_condition_kwargs(
                    instruction=instr,
                    cregmap=self.cregmap,
                    circuit=circuit,
                )

            optype = None
            if type(instr) not in (PauliEvolutionGate, UnitaryGate):
                # Handling of PauliEvolutionGate and UnitaryGate below
                optype = _optype_from_qiskit_instruction(instruction=instr)

            if optype == OpType.QControlBox:
                params = [param_to_tk(p) for p in instr.base_gate.params]
                q_ctrl_box = _get_qcontrol_box(c_gate=instr, params=params)
                self.tkc.add_qcontrolbox(q_ctrl_box, qubits)

            elif isinstance(instr, (Initialize, StatePreparation)):
                # Append OpType found by stateprep helpers
                _add_state_preparation(self.tkc, qubits, instr)

            elif type(instr) is PauliEvolutionGate:
                qpo = _qpo_from_peg(instr, qubits)
                empty_circ = Circuit(len(qargs))
                circ = gen_term_sequence_circuit(qpo, empty_circ)
                ccbox = CircBox(circ)
                self.tkc.add_circbox(ccbox, qubits)

            elif type(instr) is UnitaryGate:
                unitary = cast(NDArray[np.complex128], instr.params[0])
                if len(qubits) == 0:
                    # If the UnitaryGate acts on no qubits, we add a phase.
                    self.tkc.add_phase(np.angle(unitary[0][0]) / np.pi)
                else:
                    unitary_box = _get_unitary_box(
                        unitary=unitary, num_qubits=instr.num_qubits
                    )
                    self.tkc.add_gate(
                        unitary_box,
                        list(reversed(qubits)),
                        **condition_kwargs,
                    )

            elif optype == OpType.Barrier:
                self.tkc.add_barrier(qubits)
            elif optype == OpType.CircBox:
                qregs = (
                    [QuantumRegister(instr.num_qubits, "q")]
                    if instr.num_qubits > 0
                    else []
                )
                cregs = (
                    [ClassicalRegister(instr.num_clbits, "c")]
                    if instr.num_clbits > 0
                    else []
                )
                builder = CircuitBuilder(qregs, cregs)
                builder.add_qiskit_data(circuit, instr.definition)
                subc = builder.circuit()
                subc.name = instr.name
                self.tkc.add_circbox(CircBox(subc), qubits + bits, **condition_kwargs)  # type: ignore

            elif optype == OpType.CU3 and type(instr) is qiskit_gates.CUGate:
                if instr.params[-1] == 0:
                    self.tkc.add_gate(
                        optype,
                        [param_to_tk(p) for p in instr.params[:-1]],
                        qubits,
                        **condition_kwargs,
                    )
                else:
                    raise NotImplementedError("CUGate with nonzero phase")
            else:
                params = [param_to_tk(p) for p in instr.params]
                self.tkc.add_gate(optype, params, qubits + bits, **condition_kwargs)  # type: ignore


def qiskit_to_tk(qcirc: QuantumCircuit, preserve_param_uuid: bool = False) -> Circuit:
    """
    Converts a qiskit :py:class:`qiskit.QuantumCircuit` to a pytket :py:class:`Circuit`.

    :param qcirc: A circuit to be converted
    :param preserve_param_uuid: Whether to preserve symbolic Parameter uuids
        by appending them to the tket Circuit symbol names as "_UUID:<uuid>".
        This can be useful if you want to reassign Parameters after conversion
        to tket and back, as it is necessary for Parameter object equality
        to be preserved.
    :return: The converted circuit
    """
    circ_name = qcirc.name
    # Parameter uses a hidden _uuid for equality check
    # we optionally preserve this in parameter name for later use
    if preserve_param_uuid:
        updates = {p: Parameter(f"{p.name}_UUID:{p._uuid}") for p in qcirc.parameters}
        qcirc = cast(QuantumCircuit, qcirc.assign_parameters(updates))

    builder = CircuitBuilder(
        qregs=qcirc.qregs,
        cregs=qcirc.cregs,
        name=circ_name,
        phase=param_to_tk(qcirc.global_phase),
    )
    builder.add_qiskit_data(qcirc)
    return builder.circuit()


def _get_qiskit_control_state(bool_list: list[bool]) -> str:
    return "".join(str(int(b)) for b in bool_list)[::-1]


def param_to_tk(p: float | ParameterExpression) -> sympy.Expr:
    if isinstance(p, ParameterExpression):
        symexpr = p._symbol_expr
        try:
            return symexpr._sympy_() / sympy.pi
        except AttributeError:
            return symexpr / sympy.pi
    else:
        return p / sympy.pi


def param_to_qiskit(
    p: sympy.Expr, symb_map: dict[Parameter, sympy.Symbol]
) -> float | ParameterExpression:
    ppi = p * sympy.pi
    if len(ppi.free_symbols) == 0:
        return float(ppi.evalf())
    else:
        return ParameterExpression(symb_map, sympify(ppi))


def _get_params(
    op: Op, symb_map: dict[Parameter, sympy.Symbol]
) -> list[float | ParameterExpression]:
    return [param_to_qiskit(p, symb_map) for p in op.params]


def append_tk_command_to_qiskit(
    op: "Op",
    args: list["UnitID"],
    qcirc: QuantumCircuit,
    qregmap: dict[str, QuantumRegister],
    cregmap: dict[str, ClassicalRegister],
    symb_map: dict[Parameter, sympy.Symbol],
    range_preds: dict[Bit, tuple[list["UnitID"], int]],
) -> InstructionSet:
    optype = op.type
    if optype == OpType.Measure:
        qubit = args[0]
        bit = args[1]
        qb = qregmap[qubit.reg_name][qubit.index[0]]
        b = cregmap[bit.reg_name][bit.index[0]]
        return qcirc.measure(qb, b)

    if optype == OpType.Reset:
        qb = qregmap[args[0].reg_name][args[0].index[0]]
        return qcirc.reset(qb)

    if optype in [OpType.CircBox, OpType.ExpBox, OpType.PauliExpBox, OpType.CustomGate]:
        subcircuit = op.get_circuit()  # type: ignore
        subqc = tk_to_qiskit(subcircuit)
        qargs = []
        cargs = []
        for a in args:
            if a.type == UnitType.qubit:
                qargs.append(qregmap[a.reg_name][a.index[0]])
            else:
                cargs.append(cregmap[a.reg_name][a.index[0]])
        if optype == OpType.CustomGate:
            instruc = subqc.to_gate()
            instruc.name = op.get_name()
        else:
            instruc = subqc.to_instruction()
        return qcirc.append(instruc, qargs, cargs)
    if optype in [OpType.Unitary1qBox, OpType.Unitary2qBox, OpType.Unitary3qBox]:
        qargs = [qregmap[q.reg_name][q.index[0]] for q in args]
        u = op.get_matrix()  # type: ignore
        g = UnitaryGate(u, label="unitary")
        # Note reversal of qubits, to account for endianness (pytket unitaries are
        # ILO-BE == DLO-LE; qiskit unitaries are ILO-LE == DLO-BE).
        return qcirc.append(g, qargs=list(reversed(qargs)))
    if optype == OpType.StatePreparationBox:
        qargs = [qregmap[q.reg_name][q.index[0]] for q in args]
        statevector_array = op.get_statevector()  # type: ignore
        # check if the StatePreparationBox contains resets
        if op.with_initial_reset():  # type: ignore
            initializer = Initialize(statevector_array)
            return qcirc.append(initializer, qargs=list(reversed(qargs)))
        else:
            qiskit_state_prep_box = StatePreparation(statevector_array)
            return qcirc.append(qiskit_state_prep_box, qargs=list(reversed(qargs)))

    if optype == OpType.QControlBox:
        assert isinstance(op, QControlBox)
        qargs = [qregmap[q.reg_name][q.index[0]] for q in args]
        pytket_control_state: list[bool] = op.get_control_state_bits()
        qiskit_control_state: str = _get_qiskit_control_state(pytket_control_state)
        try:
            gatetype, phase = _known_gate_rev_phase[op.get_op().type]
        except KeyError:
            raise NotImplementedError(
                "Conversion of QControlBox with base gate"
                + f"{op.get_op()} not supported by tk_to_qiskit."
            )
        params = _get_params(op.get_op(), symb_map)
        operation = gatetype(*params)
        return qcirc.append(
            operation.control(
                num_ctrl_qubits=op.get_n_controls(), ctrl_state=qiskit_control_state
            ),
            qargs=qargs,
        )

    if optype == OpType.Barrier:
        if any(q.type == UnitType.bit for q in args):
            raise NotImplementedError(
                "Qiskit Barriers are not defined for classical bits."
            )
        qargs = [qregmap[q.reg_name][q.index[0]] for q in args]
        g = Barrier(len(args))
        return qcirc.append(g, qargs=qargs)
    if optype == OpType.RangePredicate:
        if op.lower != op.upper:  # type: ignore
            raise NotImplementedError
        range_preds[args[-1]] = (args[:-1], op.lower)  # type: ignore
        # attach predicate to bit,
        # subsequent conditional will handle it
        return Instruction("", 0, 0, [])
    if optype == OpType.Conditional:
        if op.op.type == OpType.Phase:  # type: ignore
            # conditional phase not supported
            return InstructionSet()
        if args[0] in range_preds:
            assert op.value == 1  # type: ignore
            condition_bits, value = range_preds[args[0]]  # type: ignore
            del range_preds[args[0]]  # type: ignore
            args = condition_bits + args[1:]
            width = len(condition_bits)
        else:
            width = op.width  # type: ignore
            value = op.value  # type: ignore
        regname = args[0].reg_name
        for i, a in enumerate(args[:width]):
            if a.reg_name != regname:
                raise NotImplementedError("Conditions can only use a single register")
        instruction = append_tk_command_to_qiskit(
            op.op,
            args[width:],
            qcirc,
            qregmap,
            cregmap,
            symb_map,
            range_preds,  # type: ignore
        )
        if len(cregmap[regname]) == width:
            for i, a in enumerate(args[:width]):
                if a.index != [i]:
                    raise NotImplementedError(
                        """Conditions must be an entire register in\
 order or only one bit of one register"""
                    )

            instruction.c_if(cregmap[regname], value)
        elif width == 1:
            instruction.c_if(cregmap[regname][args[0].index[0]], value)
        else:
            raise NotImplementedError(
                """Conditions must be an entire register in\
order or only one bit of one register"""
            )

        return instruction
    # normal gates
    qargs = [qregmap[q.reg_name][q.index[0]] for q in args]
    if optype == OpType.CnX:
        return qcirc.mcx(qargs[:-1], qargs[-1])
    if optype == OpType.CnY:
        return qcirc.append(qiskit_gates.YGate().control(len(qargs) - 1), qargs)
    if optype == OpType.CnZ:
        new_gate = qiskit_gates.ZGate().control(len(qargs) - 1)
        new_gate.name = "mcz"
        return qcirc.append(new_gate, qargs)
    if optype == OpType.CnRy:
        # might as well do a bit more checking
        assert len(op.params) == 1
        alpha = param_to_qiskit(op.params[0], symb_map)
        assert len(qargs) >= 2
        if len(qargs) == 2:
            # presumably more efficient; single control only
            new_gate = CRYGate(alpha)
        else:
            new_gate = RYGate(alpha).control(len(qargs) - 1)
        qcirc.append(new_gate, qargs)
        return qcirc

    if optype == OpType.CU3:
        params = _get_params(op, symb_map) + [0]
        return qcirc.append(qiskit_gates.CUGate(*params), qargs=qargs)

    if optype == OpType.TK1:
        params = _get_params(op, symb_map)
        half = ParameterExpression(symb_map, sympify(sympy.pi / 2))
        qcirc.global_phase += -params[0] / 2 - params[2] / 2
        return qcirc.append(
            qiskit_gates.UGate(params[1], params[0] - half, params[2] + half),
            qargs=qargs,
        )

    if optype == OpType.Phase:
        params = _get_params(op, symb_map)
        assert len(params) == 1
        qcirc.global_phase += params[0]
        return InstructionSet()

    # others are direct translations
    try:
        gatetype, phase = _known_gate_rev_phase[optype]
    except KeyError as error:
        raise NotImplementedError(
            "Cannot convert tket Op to Qiskit gate: " + op.get_name()
        ) from error
    params = _get_params(op, symb_map)
    g = gatetype(*params)
    if type(phase) is float:
        qcirc.global_phase += phase * np.pi
    else:
        qcirc.global_phase += sympify(phase * sympy.pi)
    return qcirc.append(g, qargs=qargs)


# The set of tket gates that can be converted directly to qiskit gates
_supported_tket_gates = set(_known_gate_rev_phase.keys())

_additional_multi_controlled_gates = {OpType.CnY, OpType.CnZ, OpType.CnRy}

# tket gates which are protected from being decomposed in the rebase
_protected_tket_gates = (
    _supported_tket_gates
    | _additional_multi_controlled_gates
    | {
        OpType.Unitary1qBox,
        OpType.Unitary2qBox,
        OpType.Unitary3qBox,
        OpType.QControlBox,
    }
    | {OpType.CustomGate}
)

# This is a rebase to the set of tket gates which have an exact substitution in qiskit
supported_gate_rebase = AutoRebase(_protected_tket_gates)


def tk_to_qiskit(
    tkcirc: Circuit, replace_implicit_swaps: bool = False
) -> QuantumCircuit:
    """
    Converts a pytket :py:class:`Circuit` to a qiskit :py:class:`qiskit.QuantumCircuit`.

    In many cases there will be a qiskit gate to exactly replace each tket gate.
    If no exact replacement can be found for a part of the circuit then an equivalent
    circuit will be returned using the tket gates which are supported in qiskit.

    :param tkcirc: A :py:class:`Circuit` to be converted
    :param replace_implicit_swaps: Implement implicit permutation by adding SWAPs
        to the end of the circuit.
    :return: The converted circuit
    """
    tkc = tkcirc.copy()  # Make a local copy of tkcirc
    if replace_implicit_swaps:
        tkc.replace_implicit_wire_swaps()
    qcirc = QuantumCircuit(name=tkc.name)
    qreg_sizes: dict[str, int] = {}
    for qb in tkc.qubits:
        if len(qb.index) != 1:
            raise NotImplementedError("Qiskit registers must use a single index")
        if (qb.reg_name not in qreg_sizes) or (qb.index[0] >= qreg_sizes[qb.reg_name]):
            qreg_sizes.update({qb.reg_name: qb.index[0] + 1})
    c_regs = tkcirc.c_registers
    if set(bit for reg in c_regs for bit in reg) != set(tkcirc.bits):
        raise NotImplementedError("Bit registers must be singly indexed from zero")
    qregmap = {}
    for reg_name, size in qreg_sizes.items():
        qis_reg = QuantumRegister(size, reg_name)
        qregmap.update({reg_name: qis_reg})
        qcirc.add_register(qis_reg)
    cregmap = {}
    for c_reg in c_regs:
        if c_reg.name != _TEMP_BIT_NAME:
            qis_reg = ClassicalRegister(c_reg.size, c_reg.name)
            cregmap.update({c_reg.name: qis_reg})
            qcirc.add_register(qis_reg)
    symb_map = {Parameter(str(s)): s for s in tkc.free_symbols()}
    range_preds: dict[Bit, tuple[list["UnitID"], int]] = dict()

    # Apply a rebase to the set of pytket gates which have replacements in qiskit
    supported_gate_rebase.apply(tkc)

    for command in tkc:
        append_tk_command_to_qiskit(
            command.op, command.args, qcirc, qregmap, cregmap, symb_map, range_preds
        )
    qcirc.global_phase += param_to_qiskit(tkc.phase, symb_map)

    # if UUID stored in name, set parameter uuids accordingly (see qiskit_to_tk)
    updates = dict()
    for p in qcirc.parameters:
        name_spl = p.name.split("_UUID:", 2)
        if len(name_spl) == 2:
            p_name, uuid_str = name_spl
            uuid = UUID(uuid_str)
            # See Parameter.__init__() in qiskit/circuit/parameter.py.
            new_p = Parameter(p_name)
            new_p._uuid = uuid
            new_p._parameter_keys = frozenset(
                ((symengine_wrapper.Symbol(p_name), uuid),)
            )
            new_p._hash = hash((new_p._parameter_keys, new_p._symbol_expr))
            updates[p] = new_p
    qcirc.assign_parameters(updates, inplace=True)

    return qcirc


def _architecture_to_couplingmap(architecture: Architecture) -> CouplingMap:
    """
    Converts a pytket Architecture object to a Qiskit CouplingMap object.

    :param architecture: Architecture to be converted
    """
    # we can make some assumptions from how the Architecture object is
    # originally constructed from the Qiskit CouplingMap:
    # 1) All nodes are single indexed
    # 2) All nodes are default register
    # 3) Node with index "i" corresponds to integer "i" in the original coupling map
    # We confirm assumption 1) and 2) while producing the coupling map
    coupling_map: List[Tuple[int, int]] = []
    for edge in architecture.coupling:
        assert len(edge[0].index) == 1
        assert len(edge[1].index) == 1
        assert edge[0].reg_name == "node"
        assert edge[1].reg_name == "node"
        coupling_map.append((edge[0].index[0], edge[1].index[0]))
    return CouplingMap(coupling_map)


def _gen_lightsabre_transformation(
    architecture: Architecture,
    optimization_level: int = 2,
    seed=0,
) -> Callable[Circuit, Circuit]:
    """
    Generates a function that can be passed to CustomPass for running
    LightSABRE routing.

    :param coupling_map: Architecture LightSABRE routes circuits to match
    :param optimization_level: Corresponds to qiskit optmization levels
    :param seed: LightSABRE routing is stochastic, with this parameter setting the seed
    """
    config: PassManagerConfig = PassManagerConfig(
        coupling_map=_architecture_to_couplingmap(architecture),
        routing_method="sabre",
        seed_transpiler=seed,
    )

    def lightsabre(circuit: Circuit) -> Circuit:
        sabre_pass: PassManager = SabreLayoutPassManager().pass_manager(
            config, optimization_level=optimization_level
        )
        c: Circuit = qiskit_to_tk(sabre_pass.run(tk_to_qiskit(circuit)))
        c.remove_blank_wires()
        c.rename_units({q: Node(q.index[0]) for q in c.qubits})
        # RebaseTket().apply(c)
        # Transform.DecomposeCXDirected(architecture).apply(c)
        return c

    return lightsabre
