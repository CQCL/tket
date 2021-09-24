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

"""Functions and template filters used for displaying circuits."""

import math
import sys
from typing import List, Tuple, Optional, Dict, Union, Any, Callable

from pytket._tket.circuit import (  # type: ignore
    Circuit,
    Command,
    Op,
    QControlBox,
    ExpBox,
    BitRegister,
    QubitRegister,
    Bit,
    Qubit,
    CustomGateDef,
)
from pytket.circuit.logic_exp import LogicExp

if sys.version_info >= (3, 8):
    from typing import TypedDict  # pylint: disable=no-name-in-module
else:
    from typing_extensions import TypedDict

# Type synonyms
Choice = Any
Chart = Dict[str, Any]
Arg = Tuple[int, Union[int, bool]]
OutcomeArray = Union[BitRegister, QubitRegister]
Register = Union[Bit, Qubit]


class ParsedOperation(TypedDict):
    op: Dict[str, Union[str, Optional[Op]]]
    params: str
    args: List[Arg]
    n_args: int


class ParsedCircuit(TypedDict):
    bits: List[Register]
    qubits: List[Register]
    layers: List[List[ParsedOperation]]


def format_complex_number(real: float, imag: float) -> str:
    """Get a string to display a given complex number nicely"""
    str_imag = (
        "" if imag == 0 else ("" if abs(imag) == 1 else str(round(abs(imag), 3))) + "i"
    )
    str_real = "" if real == 0 else str(round(real, 3))
    sign = "0" if real == 0 and imag == 0 else " + " if real * imag != 0 else ""
    # if both are defined, move imag's sign to sign.
    if imag < 0:
        sign = " - "

    return "{}{}{}".format(str_real, sign, str_imag)


def print_bitstring(width: int, bitset: List[int]) -> str:
    """Format and print a bitset of a certain width with padding."""
    return "".join(map(lambda x: bin(x)[2:].zfill(8)[:width], bitset))


def get_states_from_outcomearray(outcome_array: OutcomeArray) -> List[str]:
    """Get a list of strings representing the states associated
    with an outcome array."""
    return list(map(lambda x: print_bitstring(outcome_array.size, x), outcome_array))


def get_first_state_from_outcomearray(outcome_array: OutcomeArray) -> str:
    """Get a string representing the first state associated with an outcome array."""
    return print_bitstring(outcome_array.size, outcome_array[0])


def format_register(reg: Register) -> str:
    """Given a register return a viewer friendly version"""
    return str(reg.reg_name) + str(reg.index)


def format_mapping(
    mapping: Dict[Any, Any],
    coerce_from: Optional[str] = None,
    coerce_to: Optional[str] = None,
) -> List[str]:
    """Format a mapping for display"""
    fn_coerce: Dict[str, Callable[[Any], Any]] = {
        "register": format_register,
        "bool": lambda xs: [int(x) for x in xs],
    }
    return [
        "{} → {}".format(
            el_from if coerce_from is None else fn_coerce[coerce_from](el_from),
            el_to if coerce_to is None else fn_coerce[coerce_to](el_to),
        )
        for (el_from, el_to) in mapping.items()
    ]


def format_raw_matrix(
    unitary: List[List[complex]],
) -> Dict[str, Any]:
    """Extract the matrix from a box if applicable and format it for display."""
    n_qubits = round(math.log(len(unitary), 2))
    basis = [bin(n)[2:].zfill(n_qubits) for n in range(2 ** n_qubits)]
    return {
        "chart": zip(
            basis,
            [
                zip(basis, [format_complex_number(num.real, num.imag) for num in row])
                for row in unitary
            ],
        )
    }


def format_bool_matrix(
    matrix: List[List[bool]],
) -> Dict[str, Any]:
    """Format a boolean matrix for display."""
    n_qubits = round(math.log(len(matrix), 2))
    basis = [bin(n)[2:].zfill(n_qubits) for n in range(2 ** n_qubits)]
    return {
        "chart": zip(
            basis,
            [zip(basis, [int(val) for val in row]) for row in matrix],
        )
    }


def format_logic_exp(exp: Union[LogicExp, Bit, BitRegister, int]) -> str:
    if isinstance(exp, (int, Bit, BitRegister)):
        # Variable or constant
        return str(exp)
    elif isinstance(exp, LogicExp):
        # Recursive logic exp
        if len(exp.args) == 1:
            return "{}{}".format(
                exp.op.value,
                format_logic_exp(exp.args[0]),
            )
        else:
            return "({} {} {})".format(
                format_logic_exp(exp.args[0]),
                exp.op.value,
                format_logic_exp(exp.args[1]),
            )
    return ""


def get_gate_colour(op_type: str) -> str:
    """Get ZX-colour of the gate"""
    if op_type in {"H"}:
        return "h"
    if op_type in {
        "X",
        "Rx",
        "V",
        "Vdg",
        "S",
        "Sdg",
        "SX",
        "SXdg",
        "XXPhase",
        "PhasedX",
        "XXPhase3",
    }:
        return "x"
    if op_type in {
        "Y",
        "Ry",
        "YYPhase",
    }:
        return "y"
    if op_type in {
        "Z",
        "S",
        "Sdg",
        "T",
        "Tdg",
        "Rz",
        "PhasedISWAP",
        "ZZPhase",
        "ZZMax",
        "Measure",
        "Reset",
    }:
        return "z"
    return ""


def has_gate_info(op_type: str) -> bool:
    """Deal with gates that need non-circuit-shaped info"""
    return op_type in {
        "Unitary1qBox",
        "Unitary2qBox",
        "Unitary3qBox",
        "ExpBox",
        "PauliExpBox",
        "PhasePolyBox",
        "ClassicalExpBox",
        "Custom",
        "CircBox",
        "ConditionalGate",
        "QControlBox",
    }


def get_op_params(op: Op) -> Optional[List[str]]:
    if op.is_gate():
        return [str(param) for param in op.params]
    else:
        return None


# Deal with gates that define sub-circuits
def has_sub_circuit(op_type: str) -> bool:
    """Check whether this gate admits a sub circuit"""
    return op_type in {
        "CircBox",
        "Custom",
        "ExpBox",
        "QControlBox",
    }


def get_sub_circuit(op: Op) -> Optional[Circuit]:
    """Get the sub circuit defined for this gate"""
    try:
        return op.get_circuit()
    except Exception:
        return None


def is_classical_gate(op_type: str) -> bool:
    return op_type in {
        "ExplicitPredicate",
        "ExplicitModifier",
        "SetBits",
        "RangePredicate",
        "ClassicalTransform",
        "CopyBits",
        "MultiBit",
    }


# Deal with conditional gates
def is_control_gate(op_type: str) -> bool:
    """Get whether a gate type includes a control (q)bit."""
    return op_type in {
        "CX",
        "CY",
        "CZ",
        "CH",
        "CRx",
        "CRy",
        "CRz",
        "CU1",
        "CU3",
        "CV",
        "CVdg",
        "CSx",
        "CSXdg",
        "CSWAP",
        "CnRy",
        "CnX",
        "CCX",
        "Control",
        "QControlBox",
        "ConditionalGate",
    }


def get_target_args(args: List[Arg]) -> List[Arg]:
    """Given args for a controlled operation, get the args that are being acted on"""
    return [arg for arg in args if arg[0] < 2 and arg[1] > -1]


def get_skipped_args(args: List[Arg]) -> List[Arg]:
    """Given args for a controlled operation, get the args that are being skipped"""
    return [arg for arg in args if arg[1] < 0]


def get_op_name(op_type: str, operation: Optional[Op]) -> Tuple[str, Optional[Op]]:
    """Get the display name of the circuit operation."""

    def convert(control_type: str, op: Optional[Op]) -> Tuple[str, Optional[Op]]:
        if control_type in {
            "CX",
            "CY",
            "CZ",
            "CH",
            "CRx",
            "CRy",
            "CRz",
            "CU1",
            "CU3",
            "CV",
            "CVdg",
            "CSx",
            "CSXdg",
            "CSWAP",
        }:
            return control_type[1:], None
        elif control_type in {"CCX", "CnRy", "CnX"}:
            return control_type[2:], None

        if op is not None:
            if control_type == "QControlBox":
                sub_op = op.get_op()
                return sub_op.type.name, sub_op
            elif control_type == "ConditionalGate":
                return op.op.type.name, op.op
        return "Unknown", None

    while is_control_gate(op_type):
        op_type, operation = convert(op_type, operation)

    return op_type, operation


def get_op_display_name(
    op_type: str, op: Optional[Op], params: Optional[List[str]]
) -> str:
    """Get the display name of the circuit operation."""
    if op_type == "CircBox" and op is not None:
        name = "Circuit"
    elif op_type == "Custom" and op is not None and isinstance(op.gate, CustomGateDef):
        name = op.gate.name
    elif is_classical_gate(op_type) and op is not None:
        name = op.get_name()
    else:
        name = op_type

    return name + (str(params) if params is not None and len(params) > 0 else "")


# Controlled operations get special treatment:
# List obtained from https://cqcl.github.io/pytket/build/html/optype.html
def get_controlled_ops(op_type: str, command: Command) -> int:
    """Return the number of control bits involved in a given controlled operation"""
    if op_type in {
        "CX",
        "CY",
        "CZ",
        "CH",
        "CRx",
        "CRy",
        "CRz",
        "CU1",
        "CU3",
        "CV",
        "CVdg",
        "CSx",
        "CSXdg",
        "CSWAP",
    }:
        return 1
    elif op_type in {"CnRy", "CnX"}:
        return len(command.args) - 1
    elif op_type == "CCX":
        return 2
    elif op_type == "QControlBox":
        return int(command.op.get_n_controls())
    elif op_type == "Control" and isinstance(command.op.box, QControlBox):
        return int(command.op.box.n_controls)
    elif op_type == "ConditionalGate":
        return int(command.op.width)

    return 0


def get_expbox_details(expbox: ExpBox) -> Dict[str, Any]:
    temp_circ = Circuit(2)
    temp_circ.add_expbox(expbox, 0, 1)
    temp_dict = temp_circ.to_dict()
    info = temp_dict["commands"][0]["op"]["box"]
    details = {}
    if "matrix" in info:
        matrix = [
            [complex(entry[0], entry[1]) for entry in row] for row in info["matrix"]
        ]
        details["matrix"] = matrix
    if "phase" in info:
        details["phase"] = info["phase"]
    return details


def get_box_matrix(box: Op) -> Optional[Dict[str, Any]]:
    try:
        matrix = box.get_matrix()
    except Exception:
        return None

    return format_raw_matrix(matrix)


def format_op_params(params: Optional[list]) -> str:
    """format operation params nicely for display"""
    if params is not None:
        return "(" + ", ".join(map(str, params)) + ")" if len(params) > 0 else ""
    return ""


def make_operation(
    args: List[Arg],
    params: Optional[list] = None,
    op_type: Optional[str] = "ID",
    raw_op: Optional[Op] = None,
    id_num: Optional[str] = None,
) -> ParsedOperation:
    """Construct an operation from a command"""
    if params is None:
        params = []
    return ParsedOperation(
        op={"type": op_type, "raw": raw_op, "id": id_num},
        params=format_op_params(params),
        args=args,
        n_args=len(args),
    )


# Prepare a circuit for the templates
def parse_circuit(raw_circuit: Circuit) -> ParsedCircuit:
    """We pre-process the raw circuit for display:
    - Break the circuit down into a grid, displayed as a sequence
        of stacked blocks.
    - Each block is either an identity wire on the relevant (qu)bit,
        or [part of] a command.
    - The circuit is split into a list of (vertical) layers, each
        of which is a list of operations.
    - So each layer contains a command, possibly surrounded by
        identities on either side.
    - Each layer involves each qubit exactly once, and defines
        which blocks are to be rendered for that column.
    """
    parsed_circuit = ParsedCircuit(
        bits=raw_circuit.bits,
        qubits=raw_circuit.qubits,
        layers=[],
    )
    counter = 0  # Track how many non ID gates we have added

    # Note: arg format is (type, pos) where:
    # - type: bit = 0, qubit = 1, control bit = 2, control qubit = 3
    # - link: -1 if not included, otherwise the position in the operation.

    # Start with a layer of identities.
    id_layer = [
        make_operation(
            [(1, 1) for _ in raw_circuit.qubits] + [(0, 1) for _ in raw_circuit.bits]
        )
    ]
    parsed_circuit["layers"].append(id_layer)

    # Create a layer for each command.
    # We try to collapse layers together where we can.
    pos = 0  # keep track of which register we have filled the layer up to.
    layer = []
    registers = parsed_circuit["qubits"] + parsed_circuit["bits"]
    for command in raw_circuit.get_commands():
        # Go through each register and allocate it to the argument list for the
        # appropriate operation:
        # - Id until we get to the first (qu)bit featuring in the command
        # - and Id after we have covered the last (qu)bit in the command.
        # Identify the boundaries:
        sorted_args = sorted(command.args, key=registers.index)
        q_index = {
            "min": registers.index(sorted_args[0]),
            "max": registers.index(sorted_args[-1]),
        }

        args_top = []
        args_command = []

        # If we can't fit this command onto the current layer,
        # fill the previous with ids and prepare the next layer:
        if q_index["min"] < pos:
            # pad the rest our with ids and add the layer in.
            layer.append(
                make_operation(
                    [
                        (1 if i < len(parsed_circuit["qubits"]) else 0, 1)
                        for i in range(pos, len(registers))
                    ]
                )
            )
            parsed_circuit["layers"].append(layer)
            # reset the layer and pos
            layer = []
            pos = 0

        for i in range(pos, q_index["max"] + 1):
            register, n_is_quantum = (
                registers[i],
                1 if i < len(parsed_circuit["qubits"]) else 0,
            )

            # if this is a controlled operation, alter the wire type accordingly
            if is_control_gate(command.op.type.name):
                n_control = get_controlled_ops(command.op.type.name, command)
                if registers[i] in command.args[:n_control]:
                    n_is_quantum += 2

            if i < q_index["min"]:  # top id section
                args_top += [(n_is_quantum, 1)]
            elif i < q_index["max"] + 1:  # part of/overlapping the command
                args_command += [
                    (
                        n_is_quantum,
                        -1
                        if register not in command.args
                        else command.args.index(register),
                    )
                ]

        # Add the (up to) 3 operations to our layer:
        if len(args_top) > 0:
            layer.append(make_operation(args_top))
        if len(args_command) > 0:
            params = get_op_params(command.op)

            layer.append(
                make_operation(
                    args_command,
                    op_type=command.op.type.name,
                    params=params,
                    raw_op=command.op,
                    id_num=str(counter),
                )
            )
            counter += 1

        pos = q_index["max"] + 1

    # If we have a leftover layer, add it in
    if len(layer) > 0:
        # pad the rest our with ids and add the layer in.
        layer.append(
            make_operation(
                [
                    (1 if i < len(parsed_circuit["qubits"]) else 0, 1)
                    for i in range(pos, len(registers))
                ]
            )
        )
        parsed_circuit["layers"].append(layer)

    # Add a layer of id at the end of the circuit
    parsed_circuit["layers"].append(id_layer)
    return parsed_circuit
