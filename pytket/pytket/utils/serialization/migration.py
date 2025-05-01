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

from copy import deepcopy
from typing import Any

from pytket.circuit.logic_exp import (  # noqa: F401 pylint: disable=W0611
    BitWiseOp,
    LogicExp,
    RegWiseOp,
)


def _convert_op(exp_op: str) -> str:
    assert exp_op.startswith(("RegWiseOp.", "BitWiseOp."))
    return LogicExp.factory(eval(exp_op)).__name__


def _convert_regwise_terms(args: list, regs: list[tuple[str, int]]) -> list:
    terms = []
    for arg in args:
        assert isinstance(arg, dict)
        regname = arg.get("name")
        if regname is None:
            # It's a nested expression.
            terms.append(
                {
                    "type": "expr",
                    "input": {
                        "op": _convert_op(arg["op"]),
                        "args": _convert_regwise_terms(arg["args"], regs),
                    },
                }
            )
        else:
            # It's a register.
            i = regs.index((regname, arg["size"]))
            terms.append(
                {
                    "type": "term",
                    "input": {
                        "type": "var",
                        "term": {"type": "reg", "var": {"index": i}},
                    },
                }
            )
    return terms


def _convert_bitwise_terms(args: list, bits: list[tuple[str, tuple[int]]]) -> list:
    terms = []
    for arg in args:
        if isinstance(arg, dict):
            # It's a nested expression.
            terms.append(
                {
                    "type": "expr",
                    "input": {
                        "op": _convert_op(arg["op"]),
                        "args": _convert_bitwise_terms(arg["args"], bits),
                    },
                }
            )
        else:
            # It's a bit.
            assert isinstance(arg, list)
            assert len(arg) == 2
            name, index = arg
            i = bits.index((name, tuple(index)))
            terms.append(
                {
                    "type": "term",
                    "input": {
                        "type": "var",
                        "term": {"type": "bit", "var": {"index": i}},
                    },
                }
            )
    return terms


def _find_regs_in_expr(exp: dict[str, Any], regs: list[tuple[str, int]]) -> None:
    op = exp["op"]
    args = exp["args"]
    assert op.startswith("RegWiseOp.")
    for arg in args:
        assert isinstance(arg, dict)
        regname = arg.get("name")
        if regname is None:
            # It's a nested expression.
            _find_regs_in_expr(arg, regs)
        else:
            # It's a register.
            reg = (regname, arg["size"])
            if reg not in regs:
                regs.append(reg)


def _find_bits_in_expr(exp: dict[str, Any], bits: list[tuple[str, tuple[int]]]) -> None:
    op = exp["op"]
    args = exp["args"]
    assert op.startswith("BitWiseOp.")
    for arg in args:
        if isinstance(arg, dict):
            # It's a nested expression.
            _find_bits_in_expr(arg, bits)
        else:
            # It's a bit.
            assert isinstance(arg, list)
            bit = (arg[0], tuple(arg[1]))
            if bit not in bits:
                bits.append(bit)


def _clexpr_from_regwise_classicalexpbox(
    box: dict[str, Any], command_args: list
) -> dict[str, Any]:
    exp = box["exp"]
    regs: list[tuple[str, int]] = []
    _find_regs_in_expr(exp, regs)  # construct ordered list of reg vars
    reg_posn = [
        [i, [command_args.index([name, [j]]) for j in range(size)]]
        for i, (name, size) in enumerate(regs)
    ]
    terms = _convert_regwise_terms(exp["args"], regs)
    return {
        "expr": {"op": _convert_op(exp["op"]), "args": terms},
        "bit_posn": [],
        "reg_posn": reg_posn,
        "output_posn": list(range(box["n_i"], len(command_args))),
    }


def _clexpr_from_bitwise_classicalexpbox(
    box: dict[str, Any], command_args: list
) -> dict[str, Any]:
    exp = box["exp"]
    bits: list[tuple[str, tuple[int]]] = []
    _find_bits_in_expr(exp, bits)  # construct ordered list of bit vars
    bit_posn = [
        [i, command_args.index([name, list(index)])]
        for i, (name, index) in enumerate(bits)
    ]
    terms = _convert_bitwise_terms(exp["args"], bits)
    return {
        "expr": {"op": _convert_op(exp["op"]), "args": terms},
        "bit_posn": bit_posn,
        "reg_posn": [],
        "output_posn": list(range(box["n_i"], len(command_args))),
    }


def _clexpr_from_classicalexpbox(
    box: dict[str, Any], command_args: list
) -> dict[str, Any]:
    exp_op = box["exp"]["op"]
    if exp_op.startswith("RegWiseOp."):
        return _clexpr_from_regwise_classicalexpbox(box, command_args)
    assert exp_op.startswith("BitWiseOp.")
    return _clexpr_from_bitwise_classicalexpbox(box, command_args)


def _new_op(op: dict[str, Any], args: list) -> dict[str, Any]:
    match op["type"]:
        case "Conditional":
            new_data = deepcopy(op)
            new_data["conditional"]["op"] = _new_op(
                op["conditional"]["op"], args[op["conditional"]["width"] :]
            )
            return new_data
        case "QControlBox":
            new_data = deepcopy(op)
            new_data["box"]["op"] = _new_op(op["box"]["op"], args)
            return new_data
        case "CircBox":
            new_data = deepcopy(op)
            new_data["box"]["circuit"] = circuit_dict_from_pytket1_dict(
                op["box"]["circuit"]
            )
            return new_data
        case "CustomGate":
            new_data = deepcopy(op)
            new_data["box"]["gate"] = circuit_dict_from_pytket1_dict(op["box"]["gate"])
            return new_data
        case "ClassicalExpBox":
            return {
                "type": "ClExpr",
                "expr": _clexpr_from_classicalexpbox(op["box"], args),
            }
        case _:
            return op


def _new_command(command: dict[str, Any]) -> dict[str, Any]:
    new_data = deepcopy(command)
    new_data["op"] = _new_op(command["op"], command["args"])
    return new_data


def circuit_dict_from_pytket1_dict(circuit_data: dict[str, Any]) -> dict[str, Any]:
    """Update the serialization of a pytket 1 circuit to an equivalent pytket circuit.
    (This converts ClassicalExpBox to ClExprOp operations.)

    :param circuit_data: serialization of pytket 1 circuit
    :return: serialization of equivalent pytket circuit
    """
    # Replace all ClassicalExpBox ops. Need to recurse into CircBoxes etc.
    new_data = deepcopy(circuit_data)
    new_data["commands"] = [
        _new_command(command) for command in circuit_data["commands"]
    ]
    return new_data
