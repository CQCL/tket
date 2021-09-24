# Copyright 2019-2021 Cambridge Quantum Computing
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

# TODO: Output custom gates
# TODO: Figure out nice way to make these class methods of Circuit
import io
import os
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    TextIO,
    Tuple,
    Type,
    TypeVar,
    Union,
)
from itertools import groupby
from sympy import sympify, pi  # type: ignore
from pytket import Circuit, OpType, Qubit, Bit
from pytket.circuit import (  # type: ignore
    CustomGateDef,
    UnitID,
    BitRegister,
    QubitRegister,
    Op,
)

NOPARAM_COMMANDS = {
    "CX": OpType.CX,  # built-in gate equivalent to "cx"
    "cx": OpType.CX,
    "x": OpType.X,
    "y": OpType.Y,
    "z": OpType.Z,
    "h": OpType.H,
    "s": OpType.S,
    "sdg": OpType.Sdg,
    "t": OpType.T,
    "tdg": OpType.Tdg,
    "sx": OpType.SX,
    "sxdg": OpType.SXdg,
    "cz": OpType.CZ,
    "cy": OpType.CY,
    "ch": OpType.CH,
    "csx": OpType.CSX,
    "ccx": OpType.CCX,
    "ZZ": OpType.ZZMax,
    "measure": OpType.Measure,
    "reset": OpType.Reset,
    "id": OpType.noop,
    "barrier": OpType.Barrier,
    "swap": OpType.SWAP,
    "cswap": OpType.CSWAP,
    "ecr": OpType.ECR,
}

PARAM_COMMANDS = {
    "U": OpType.U3,  # built-in gate equivalent to "u3"
    "u3": OpType.U3,
    "u2": OpType.U2,
    "u1": OpType.U1,
    "rx": OpType.Rx,
    "ry": OpType.Ry,
    "rz": OpType.Rz,
    "Rz": OpType.Rz,
    "U1q": OpType.PhasedX,
    "crz": OpType.CRz,
    "crx": OpType.CRx,
    "cry": OpType.CRy,
    "cu1": OpType.CU1,
    "cu3": OpType.CU3,
}

included_gates = {
    "qelib1": set(
        (
            "CX",
            "cx",
            "x",
            "y",
            "z",
            "h",
            "s",
            "sdg",
            "t",
            "tdg",
            "sx",
            "sxdg",
            "cz",
            "cy",
            "ch",
            "csx",
            "ccx",
            "measure",
            "reset",
            "id",
            "barrier",
            "U",
            "u3",
            "u2",
            "u1",
            "rx",
            "ry",
            "rz",
            "crz",
            "crx",
            "cry",
            "cu1",
            "cu3",
            "swap",
            "cswap",
            "ecr",
        )
    ),
    "oqc": set(
        (
            "sx",
            "rz",
            "ecr",
            "barrier",
            "measure",
        )
    ),
}
included_gates["hqslib1"] = included_gates["qelib1"].copy()
included_gates["hqslib1"].update(("U1q", "rz", "ZZ"))
included_gates["hqslib1"].difference_update(
    ("crx", "cry", "sx", "sxdg", "csx", "swap", "cswap")
)
_tk_to_qasm_noparams = dict(((item[1], item[0]) for item in NOPARAM_COMMANDS.items()))
_tk_to_qasm_noparams[OpType.CX] = "cx"  # prefer "cx" to "CX"
_tk_to_qasm_params = dict(((item[1], item[0]) for item in PARAM_COMMANDS.items()))
_tk_to_qasm_params[OpType.U3] = "u3"  # prefer "u3" to "U"
_tk_to_qasm_params[OpType.Rz] = "rz"  # prefer "rz" to "Rz"

_classical_gatestr_map = {"AND": "&", "OR": "|", "XOR": "^"}


class QASMUnsupportedError(Exception):
    pass


class QASMParseError(Exception):
    pass


class QASMParser(object):
    """Class for parsing OpenQASM files into CQC tket Circuits."""

    def __init__(self) -> None:
        self.circuit = Circuit()
        self.gate_dict: Dict[str, CustomGateDef] = dict()
        self.reg_map: Dict[str, UnitID] = dict()
        self.include = ""

    def parse_qasm(self, qasm: str) -> Circuit:
        lines = qasm.splitlines()
        rows = []

        # first, get rid of comments and whitespace lines
        for l in lines:
            i = l.find("//")
            if i != -1:
                s = l[0:i].strip()
            else:
                s = l.strip()
            if s:
                rows.append(s)

        # now, throw away OPENQASM descriptor etc.
        if not (
            rows[0].startswith("OPENQASM 2.0")
            and rows[1].startswith('include "')
            and rows[1].endswith('.inc";')
        ):
            raise QASMParseError("File must declare OPENQASM version and its includes.")
        self.include = rows[1][len('include "') : -len('".inc;')]
        if self.include not in ("qelib1", "hqslib1"):
            raise QASMParseError("Header {}.inc not recognised".format(self.include))
        data = "\n".join(rows[2:])

        # now, separate out the custom gates to deal with elsewhere
        while True:
            i = data.find("gate ")
            if i == -1:
                break
            j = data.find("}", i)
            if j == -1:
                raise QASMParseError("Custom gate definition is invalid.")
            self.parse_custom_gate(data[i : j + 1])  # TODO: deal with custom gate
            data = data[:i] + data[j + 1 :]

        # now, parse the regular instructions
        instructions: List[str] = [s.strip() for s in data.split(";") if s.strip()]
        for instr in instructions:
            self.parse_instruction(instr, self.circuit, self.reg_map)
        return self.circuit

    def parse_custom_gate(self, data: str) -> None:
        signature, rest = data.split("{", 1)
        _, signature = signature.split(" ", 1)  # ignore "gate"
        if signature.find("(") != -1:
            gatename, other = signature.split("(")
            symbol_list, arg_list = other.split(")")
        else:
            gatename, arg_list = signature.split(" ", 1)
            symbol_list = ""
        gatename = gatename.strip()
        symbols = [sympify(s.strip()) for s in symbol_list.split(",")]
        args = [a.strip() for a in arg_list.split(",")]
        rename_map = {}
        qb_map = {}
        circ = Circuit()
        for i, a in enumerate(args):
            circ.add_qubit(Qubit(a))
            rename_map.update({Qubit(a): Qubit(i)})
            qb_map[a] = [Qubit(a)]
        command_block, _ = rest.split("}", 1)
        commands = [c.strip() for c in command_block.split(";") if c.strip()]
        for com in commands:
            self.parse_instruction(com, circ, qb_map)
        circ.rename_units(rename_map)
        symbol_map = {sym: sym * pi for sym in symbols}
        circ.symbol_substitution(symbol_map)  # qasm arguments are given in radians
        self.gate_dict[gatename] = CustomGateDef.define(gatename, circ, symbols)

    def parse_instruction(
        self, instruction: str, circuit: Circuit, reg_map: Dict[str, List[UnitID]]
    ) -> None:
        gate_kwargs: Dict[str, Any] = {}
        if instruction.find("if") == 0:
            ###parse condition
            if_phrase, rest = instruction.split("(", 1)
            if if_phrase.strip() != "if":
                raise QASMParseError(
                    'Error in parsing: cannot match "{}" against "if"'.format(if_phrase)
                )
            condition, rest = rest.split(")", 1)
            creg, eq_value = condition.split("==", 1)
            gate_kwargs.update({"condition_bits": reg_map[creg.strip()]})
            value = int(eq_value.strip())
            gate_kwargs.update({"condition_value": value})
            instruction = rest.strip()
        if instruction.find("->") != -1:
            ###handle measure gates
            ###currently assumes that there is just 1 qb being read to 1 bit
            name_and_qbs, bits = instruction.split("->", 1)
            if name_and_qbs.find("measure") == -1:
                raise QASMParseError(
                    "Error in parsing: cannot accept a non-Measure gate writing to "
                    "classical register"
                )
            name_and_qbs = name_and_qbs.replace("measure", "")
            name_and_qbs = name_and_qbs.replace(" ", "")

            name_and_qbs.strip()
            qubits_list: List[Bit]
            if "[" in name_and_qbs:
                qregname, qbindex = name_and_qbs.split("[")
                qbindex, _ = qbindex.split("]")
                qubits_list = [Qubit(qregname, int(qbindex))]
            else:
                qubits_list = reg_map[name_and_qbs]

            bits = bits.replace(" ", "")
            bits_list: List[Bit]
            if "[" in bits:
                bitreg, bitindex = bits.split("[")
                bitindex, _ = bitindex.split("]")
                bits_list = [Bit(bitreg, int(bitindex))]
            else:
                bits_list = reg_map[bits]

            for q, b in zip(qubits_list, bits_list):
                circuit.Measure(q, b, **gate_kwargs)
            return

        index = _find_respecting_brackets(instruction, " ")
        name = instruction[:index]
        rest = instruction[index + 1 :]
        args = [s.strip() for s in rest.split(",") if s.strip()]

        # deal with qubit register declarations
        if name == "qreg" or name == "creg":
            regname, rest = args[0].split("[", 1)
            regname.strip()
            size = int(rest[:-1])
            if name == "qreg":
                dict_map = circuit.add_q_register(regname, size)
            else:
                dict_map = circuit.add_c_register(regname, size)
            reg_map[regname] = [dict_map[i] for i in range(size)]
            return

        # get qubits to append operation to
        qubits = []
        for a in args:
            if "[" in a:
                regname, rest = a.split("[", 1)
                val = int(rest[:-1])
                qubits.append([Qubit(regname, val)])
            else:
                qubits.append(reg_map[a])

        # if the gate is parameterised, get these parameters
        if name.find("(") != -1:
            name, params = name.split("(", 1)
            params = params[:-1]  # cut off final close bracket
            angle_start = 0
            angle_end = _find_respecting_brackets(params, ",")
            angles = []
            while angle_end != -1:
                angles.append(params[angle_start:angle_end].strip())
                angle_start = angle_end + 1
                angle_end = _find_respecting_brackets(params, ",", angle_start)
            angles.append(params[angle_start:].strip())
            halfturn_angles = []
            for ang in angles:
                try:
                    halfturns = sympify(ang) / pi
                    halfturn_angles.append(halfturns)
                except:
                    raise QASMParseError("Cannot parse angle: {}".format(ang))
            if name in PARAM_COMMANDS:
                if (
                    self.include != "hqslib1"
                    and name in included_gates["hqslib1"]
                    and name not in included_gates["qelib1"]
                ):
                    raise QASMParseError(
                        "Gate of type {} is not defined in header {}.inc".format(
                            name, self.include
                        )
                    )
                for qbs in zip(*qubits):
                    circuit.add_gate(
                        PARAM_COMMANDS[name], halfturn_angles, list(qbs), **gate_kwargs
                    )
            elif name in self.gate_dict:
                for qbs in zip(*qubits):
                    circuit.add_custom_gate(
                        self.gate_dict[name], halfturn_angles, list(qbs), **gate_kwargs
                    )
            else:
                raise QASMParseError("Cannot parse gate of type: {}".format(name))

        else:
            if name == "barrier":
                circuit.add_barrier([q for qbs in qubits for q in qbs])
            elif name in NOPARAM_COMMANDS:
                if (
                    self.include != "hqslib1"
                    and name in included_gates["hqslib1"]
                    and name not in included_gates["qelib1"]
                ):
                    raise QASMParseError(
                        "Gate of type {} is not defined in header {}.inc".format(
                            name, self.include
                        )
                    )

                for qbs in zip(*qubits):
                    circuit.add_gate(
                        NOPARAM_COMMANDS[name], [], list(qbs), **gate_kwargs
                    )
            elif name in self.gate_dict:
                for qbs in zip(*qubits):
                    circuit.add_custom_gate(
                        self.gate_dict[name], [], list(qbs), **gate_kwargs
                    )
            else:
                raise QASMParseError("Cannot parse gate of type: {}".format(name))


def circuit_from_qasm(input_file: Union[str, "os.PathLike[Any]"]) -> Circuit:
    """A method to generate a tket Circuit from a qasm file"""
    ext = os.path.splitext(input_file)[-1]
    if ext != ".qasm":
        raise TypeError("Can only convert .qasm files")
    with open(input_file, "r") as f:
        circ = circuit_from_qasm_io(f)
    return circ


def circuit_from_qasm_str(qasm_str: str) -> Circuit:
    """A method to generate a tket Circuit from a qasm str"""
    p = QASMParser()
    return p.parse_qasm(qasm_str)


def circuit_from_qasm_io(stream_in: TextIO) -> Circuit:
    """A method to generate a tket Circuit from a qasm text stream"""
    return circuit_from_qasm_str(stream_in.read())


def circuit_to_qasm(circ: Circuit, output_file: str, header: str = "qelib1") -> None:
    """A method to generate a qasm file from a tket Circuit"""
    with open(output_file, "w") as out:
        circuit_to_qasm_io(circ, out, header=header)


def circuit_to_qasm_str(circ: Circuit, header: str = "qelib1") -> str:
    """A method to generate a qasm str from a tket Circuit"""
    buffer = io.StringIO()
    circuit_to_qasm_io(circ, buffer, header=header)
    return buffer.getvalue()


TypeReg = TypeVar("TypeReg", BitRegister, QubitRegister)


def _retrieve_registers(
    units: List[UnitID], reg_type: Type[TypeReg]
) -> Dict[str, TypeReg]:
    if any(len(unit.index) != 1 for unit in units):
        raise NotImplementedError("OPENQASM registers must use a single index")

    maxunits = map(
        lambda x: max(x[1]), groupby(units, key=lambda un: un.reg_name)  # type:ignore
    )
    return {
        maxunit.reg_name: reg_type(maxunit.reg_name, maxunit.index[0] + 1)
        for maxunit in maxunits
    }


def _parse_range(minval: int, maxval: int) -> Tuple[str, int]:
    REGMAX = (1 << 32) - 1
    if minval == maxval:
        return ("==", minval)
    elif minval == 0:
        return ("<=", maxval)
    elif maxval == REGMAX:
        return (">=", minval)
    else:
        raise NotImplementedError("Range can only be bounded on one side.")


def _get_optype_and_params(op: Op) -> Tuple[OpType, Optional[List[float]]]:
    optype = op.type
    params = op.params if optype in _tk_to_qasm_params else None
    if optype == OpType.TK1:
        # convert to U3
        optype = OpType.U3
        params = [op.params[1], op.params[0] - 0.5, op.params[2] + 0.5]
    return (optype, params)


def circuit_to_qasm_io(
    circ: Circuit, stream_out: TextIO, header: str = "qelib1"
) -> None:
    """A method to generate a qasm text stream from a tket Circuit"""
    if (
        any(
            circ.n_gates_of_type(typ)
            for typ in (
                OpType.RangePredicate,
                OpType.MultiBit,
                OpType.ExplicitPredicate,
                OpType.ExplicitModifier,
                OpType.SetBits,
            )
        )
        and header != "hqslib1"
    ):
        raise QASMUnsupportedError(
            "Complex classical gates only supported with hqslib1."
        )
    if header == "oqc":
        stream_out.write('OPENQASM 2.0;\ninclude "qelib1.inc";\n\n')
    else:
        stream_out.write('OPENQASM 2.0;\ninclude "{}.inc";\n\n'.format(header))
    qregs = _retrieve_registers(circ.qubits, QubitRegister)
    cregs = _retrieve_registers(circ.bits, BitRegister)

    for reg in qregs.values():
        stream_out.write(f"qreg {reg.name}[{reg.size}];\n")
    for reg in cregs.values():
        stream_out.write(f"creg {reg.name}[{reg.size}];\n")

    range_preds = dict()
    for command in circ:
        op = command.op
        args = command.args
        optype, params = _get_optype_and_params(op)
        if optype == OpType.RangePredicate:
            range_preds[args[-1]] = command
            # attach predicate to bit,
            # subsequent conditional will handle it
            continue
        if optype == OpType.ConditionalGate:
            bits = args[: op.width]
            control_bit = bits[0]
            if control_bit in range_preds:
                # write range predicate in condition
                range_com = range_preds[control_bit]
                range_op = range_com.op
                comparator, value = _parse_range(range_op.lower, range_op.upper)
                if op.value == 0 and comparator == "==":
                    comparator = "!="
                if header != "hqslib1" and comparator != "==":
                    raise QASMUnsupportedError(
                        "OpenQASM conditions must be on a register's fixed value."
                    )
                bits = range_com.args[:-1]
                variable = range_com.args[0].reg_name
                del range_preds[control_bit]
            else:
                comparator = "=="
                value = op.value
                variable = control_bit if header == "hqslib1" else control_bit.reg_name
            if header == "hqslib1":
                if op.width != 1 or op.value not in (0, 1):
                    raise QASMUnsupportedError(
                        "HQS OpenQASM conditions must act on one bit."
                    )
            else:
                if op.width != cregs[variable].size:
                    raise QASMUnsupportedError(
                        "OpenQASM conditions must be an entire classical register"
                    )
                if sorted(bits) != list(cregs[variable]):
                    raise QASMUnsupportedError(
                        "OpenQASM conditions must be a single classical register"
                    )

            stream_out.write(f"if({variable}{comparator}{value}) ")
            args = args[op.width :]
            op = op.op
            optype, params = _get_optype_and_params(op)
        if optype == OpType.SetBits:
            creg_name = args[0].reg_name
            bits, vals = zip(*sorted(zip(args, op.values)))

            # check if whole register can be set at once
            if bits == tuple(cregs[creg_name]):
                value = int("".join(map(str, map(int, vals[::-1]))), 2)
                stream_out.write(f"{creg_name} = {value};\n")
            else:
                for bit, value in zip(bits, vals):
                    stream_out.write(f"{bit} = {int(value)};\n")
            continue
        if optype == OpType.MultiBit:
            op = op.basic_op
            optype = op.type
            registers_involved = [arg.reg_name for arg in args[:2]]
            if args[2].reg_name not in registers_involved:
                registers_involved.append(args[2].reg_name)
            args = [cregs[name] for name in registers_involved]
        if optype in (
            OpType.ExplicitPredicate,
            OpType.ExplicitModifier,
        ):
            # &, ^ and | gates
            opstr = str(op)
            if opstr not in _classical_gatestr_map:
                raise QASMUnsupportedError(f"Classical gate {opstr} not supported.")
            stream_out.write(
                f"{args[-1]} = {args[0]} {_classical_gatestr_map[opstr]} {args[1]};\n"
            )
            continue
        if optype in _tk_to_qasm_noparams:
            opstr = _tk_to_qasm_noparams[optype]
        elif optype in _tk_to_qasm_params:
            opstr = _tk_to_qasm_params[optype]
        else:
            raise QASMUnsupportedError(
                "Cannot print command of type: {}".format(op.get_name())
            )
        if opstr not in included_gates[header]:
            raise QASMUnsupportedError(
                "Gate of type {} is not defined in header {}.inc".format(opstr, header)
            )
        stream_out.write(opstr)
        if params is not None:
            stream_out.write("(")
            for i in range(len(params)):
                reduced = True
                try:
                    p = float(params[i])
                except TypeError:
                    reduced = False
                    p = params[i]
                if i < len(params) - 1:
                    if reduced:
                        stream_out.write("{}*pi,".format(p))
                    else:
                        stream_out.write("({})*pi,".format(p))

                else:
                    if reduced:
                        stream_out.write("{}*pi)".format(p))
                    else:
                        stream_out.write("({})*pi)".format(p))
        stream_out.write(" ")
        if optype == OpType.Measure:  # assume written to only 1 bit
            stream_out.write(
                "{q} -> {c};\n".format(q=args[0].__repr__(), c=args[1].__repr__())
            )
        else:
            for i in range(len(args)):
                stream_out.write(args[i].__repr__())
                if i < len(args) - 1:
                    stream_out.write(",")
                else:
                    stream_out.write(";\n")


def _find_respecting_brackets(full_string: str, phrase: str, start: int = 0) -> int:
    """Assuming `full_string` is well-bracketed (at no point is there a close
    without an unmatched open), returns the index of the first location of
    `phrase` in `full_string` that is not inside any pair of brackets
    """
    length = len(full_string)
    non_neg_fix: Callable[[int], int] = lambda x: (length if x == -1 else x)
    next_phrase = full_string.find(phrase, start)
    next_open = non_neg_fix(full_string.find("(", start))
    next_close = non_neg_fix(full_string.find(")", start))
    depth = 0
    while next_phrase != -1:
        if next_phrase < next_open and next_phrase < next_close:
            if depth == 0:
                return next_phrase  # found a match
            else:
                next_phrase = full_string.find(
                    phrase, next_phrase + 1
                )  # bad match, try next
        elif next_open < next_close:
            depth += 1
            next_open = non_neg_fix(full_string.find("(", next_open + 1))
        else:
            # there must be a close first, as
            # length > next_phrase >= next_open >= next_close
            depth -= 1
            next_close = non_neg_fix(full_string.find(")", next_close + 1))
    return -1
