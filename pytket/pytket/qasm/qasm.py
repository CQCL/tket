# Copyright 2019-2023 Cambridge Quantum Computing
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
import io
import os
import re
import uuid

# TODO: Output custom gates
from importlib import import_module
from itertools import chain, groupby
from decimal import Decimal
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    Iterable,
    Iterator,
    List,
    NewType,
    Optional,
    Sequence,
    Set,
    TextIO,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)
from sympy import Symbol, pi  # type: ignore
from lark import Discard, Lark, Token, Transformer, Tree

from pytket._tket.circuit import (
    ClassicalExpBox,
    Command,
    Conditional,
    RangePredicateOp,
    SetBitsOp,
    CopyBitsOp,
    MultiBitOp,
    WASMOp,
    CustomGate,
    MetaOp,
)
from pytket._tket.unit_id import _TEMP_BIT_NAME
from pytket.circuit import (
    Bit,
    BitRegister,
    Circuit,
    Op,
    OpType,
    Qubit,
    QubitRegister,
    UnitID,
)
from pytket.circuit.decompose_classical import int_to_bools
from pytket.circuit.logic_exp import (
    BitLogicExp,
    BitWiseOp,
    PredicateExp,
    LogicExp,
    RegEq,
    RegLogicExp,
    RegNeg,
    RegWiseOp,
)
from pytket.qasm.grammar import grammar
from pytket.passes import auto_rebase_pass, RemoveRedundancies
from pytket.wasm import WasmFileHandler


class QASMParseError(Exception):
    """Error while parsing QASM input."""

    def __init__(
        self, msg: str, line: Optional[int] = None, fname: Optional[str] = None
    ):
        self.msg = msg
        self.line = line
        self.fname = fname

        ctx = "" if fname is None else f"\nFile:{fname}: "
        ctx += "" if line is None else f"\nLine:{line}. "

        super().__init__(f"{msg}{ctx}")


class QASMUnsupportedError(Exception):
    pass


Value = Union[int, float, str]
T = TypeVar("T")

_BITOPS = set(op.value for op in BitWiseOp)
_BITOPS.update(("+", "-"))  # both are parsed to XOR
_REGOPS = set(op.value for op in RegWiseOp)

Arg = Union[List, str]


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
    "c3x": OpType.CnX,
    "c4x": OpType.CnX,
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
    "p": OpType.U1,  # alias. https://github.com/Qiskit/qiskit-terra/pull/4765
    "u": OpType.U3,  # alias. https://github.com/Qiskit/qiskit-terra/pull/4765
    "U": OpType.U3,  # built-in gate equivalent to "u3"
    "u3": OpType.U3,
    "u2": OpType.U2,
    "u1": OpType.U1,
    "rx": OpType.Rx,
    "rxx": OpType.XXPhase,
    "ry": OpType.Ry,
    "rz": OpType.Rz,
    "RZZ": OpType.ZZPhase,
    "rzz": OpType.ZZPhase,
    "Rz": OpType.Rz,
    "U1q": OpType.PhasedX,
    "crz": OpType.CRz,
    "crx": OpType.CRx,
    "cry": OpType.CRy,
    "cu1": OpType.CU1,
    "cu3": OpType.CU3,
    "Rxxyyzz": OpType.TK2,
}

NOPARAM_EXTRA_COMMANDS = {
    "v": OpType.V,
    "vdg": OpType.Vdg,
    "cv": OpType.CV,
    "cvdg": OpType.CVdg,
    "csxdg": OpType.CSXdg,
    "bridge": OpType.BRIDGE,
    "iswapmax": OpType.ISWAPMax,
    "zzmax": OpType.ZZMax,
    "ecr": OpType.ECR,
}

PARAM_EXTRA_COMMANDS = {
    "tk2": OpType.TK2,
    "iswap": OpType.ISWAP,
    "phasediswap": OpType.PhasedISWAP,
    "yyphase": OpType.YYPhase,
    "xxphase3": OpType.XXPhase3,
    "eswap": OpType.ESWAP,
    "fsim": OpType.FSim,
}

_tk_to_qasm_noparams = dict(((item[1], item[0]) for item in NOPARAM_COMMANDS.items()))
_tk_to_qasm_noparams[OpType.CX] = "cx"  # prefer "cx" to "CX"
_tk_to_qasm_params = dict(((item[1], item[0]) for item in PARAM_COMMANDS.items()))
_tk_to_qasm_params[OpType.U3] = "u3"  # prefer "u3" to "U"
_tk_to_qasm_params[OpType.Rz] = "rz"  # prefer "rz" to "Rz"
_tk_to_qasm_extra_noparams = dict(
    ((item[1], item[0]) for item in NOPARAM_EXTRA_COMMANDS.items())
)
_tk_to_qasm_extra_params = dict(
    ((item[1], item[0]) for item in PARAM_EXTRA_COMMANDS.items())
)

_classical_gatestr_map = {"AND": "&", "OR": "|", "XOR": "^"}


_all_known_gates = (
    set(NOPARAM_COMMANDS.keys())
    .union(PARAM_COMMANDS.keys())
    .union(PARAM_EXTRA_COMMANDS.keys())
    .union(NOPARAM_EXTRA_COMMANDS.keys())
)
_all_string_maps = {
    key: val.name
    for key, val in chain(
        PARAM_COMMANDS.items(),
        NOPARAM_COMMANDS.items(),
        PARAM_EXTRA_COMMANDS.items(),
        NOPARAM_EXTRA_COMMANDS.items(),
    )
}

unit_regex = re.compile(r"([a-z][a-zA-Z0-9_]*)\[([\d]+)\]")


def _extract_reg(var: Token) -> Tuple[str, int]:
    match = unit_regex.match(var.value)
    assert match is not None
    return match.group(1), int(match.group(2))


def _load_include_module(
    header_name: str, flter: bool, decls_only: bool
) -> Dict[str, Dict]:
    try:
        if decls_only:
            include_def: Dict[str, Dict] = import_module(
                f"pytket.qasm.includes._{header_name}_decls"
            )._INCLUDE_DECLS
        else:
            include_def = import_module(
                f"pytket.qasm.includes._{header_name}_defs"
            )._INCLUDE_DEFS
    except ModuleNotFoundError as e:
        raise QASMParseError(
            f"Header {header_name} is not known and cannot be loaded."
        ) from e
    return {
        gate: include_def[gate]
        for gate in include_def
        if not flter or gate not in _all_known_gates
    }


def _bin_par_exp(op: "str") -> Callable[["CircuitTransformer", List[str]], str]:
    def f(self: "CircuitTransformer", vals: List[str]) -> str:
        return f"({vals[0]} {op} {vals[1]})"

    return f


def _un_par_exp(op: "str") -> Callable[["CircuitTransformer", List[str]], str]:
    def f(self: "CircuitTransformer", vals: List[str]) -> str:
        return f"({op}{vals[0]})"

    return f


def _un_call_exp(op: "str") -> Callable[["CircuitTransformer", List[str]], str]:
    def f(self: "CircuitTransformer", vals: List[str]) -> str:
        return f"{op}({vals[0]})"

    return f


def _hashable_uid(arg: List) -> Tuple[str, int]:
    return arg[0], arg[1][0]


Reg = NewType("Reg", str)
CommandDict = Dict[str, Any]


@dataclass
class ParsMap:
    pars: Iterable[str]

    def __iter__(self) -> Iterable[str]:
        return self.pars


class CircuitTransformer(Transformer):
    def __init__(self, return_gate_dict: bool = False) -> None:
        super().__init__()
        self.q_registers: Dict[str, int] = {}
        self.c_registers: Dict[str, int] = {}
        self.gate_dict: Dict[str, Dict] = {}
        self.wasm: Optional[WasmFileHandler] = None
        self.include = ""
        self.return_gate_dict = return_gate_dict

    def _fresh_temp_bit(self) -> List:
        if _TEMP_BIT_NAME in self.c_registers:
            idx = self.c_registers[_TEMP_BIT_NAME]
        else:
            idx = 0
        self.c_registers[_TEMP_BIT_NAME] = idx + 1

        return [_TEMP_BIT_NAME, [idx]]

    def _reset_context(self, reset_wasm: bool = True) -> None:
        self.q_registers = {}
        self.c_registers = {}
        self.gate_dict = {}
        self.include = ""
        if reset_wasm:
            self.wasm = None

    def _get_reg(self, name: str) -> Reg:
        return Reg(name)

    def _get_uid(self, iarg: Token) -> List:
        name, idx = _extract_reg(iarg)
        return [name, [idx]]

    def _get_arg(self, arg: Token) -> Arg:
        if arg.type == "IARG":
            return self._get_uid(arg)
        else:
            return self._get_reg(arg.value)

    def unroll_all_args(self, args: Iterable[Arg]) -> Iterator[List[Any]]:
        for arg in args:
            if isinstance(arg, str):
                size = (
                    self.q_registers[arg]
                    if arg in self.q_registers
                    else self.c_registers[arg]
                )
                yield [[arg, [idx]] for idx in range(size)]
            else:
                yield [arg]

    def margs(self, tree: Iterable[Token]) -> Iterator[Arg]:
        return map(self._get_arg, tree)

    def iargs(self, tree: Iterable[Token]) -> Iterator[List]:
        return map(self._get_uid, tree)

    def args(self, tree: Iterable[Token]) -> Iterator[List]:
        return ([tok.value, [0]] for tok in tree)

    def creg(self, tree: List[Token]) -> None:
        name, size = _extract_reg(tree[0])
        self.c_registers[Reg(name)] = size

    def qreg(self, tree: List[Token]) -> None:
        name, size = _extract_reg(tree[0])
        self.q_registers[Reg(name)] = size

    def meas(self, tree: List[Token]) -> Iterable[CommandDict]:
        for args in zip(*self.unroll_all_args(self.margs(tree))):
            yield {"args": list(args), "op": {"type": "Measure"}}

    def barr(self, tree: List[Arg]) -> Iterable[CommandDict]:
        args = [q for qs in self.unroll_all_args(tree[0]) for q in qs]
        yield {
            "args": args,
            "op": {"signature": ["Q"] * len(args), "type": "Barrier"},
        }

    def reset(self, tree: List[Token]) -> Iterable[CommandDict]:
        for qb in next(self.unroll_all_args(self.margs(tree))):
            yield {"args": [qb], "op": {"type": "Reset"}}

    def pars(self, vals: Iterable[str]) -> ParsMap:
        return ParsMap(map(str, vals))

    def mixedcall(self, tree: List) -> Iterator[CommandDict]:
        child_iter = iter(tree)

        optoken = next(child_iter)
        opstr = optoken.value

        next_tree = next(child_iter)
        try:
            args = next(child_iter)
            pars = cast(ParsMap, next_tree).pars
        except StopIteration:
            args = next_tree
            pars = []

        treat_as_barrier = [
            "sleep",
            "order2",
            "order3",
            "order4",
            "order5",
            "order6",
            "order7",
            "order8",
            "order9",
            "order10",
            "order11",
            "order12",
            "order13",
            "order14",
            "order15",
            "order16",
            "order17",
            "order18",
            "order19",
            "order20",
            "group2",
            "group3",
            "group4",
            "group5",
            "group6",
            "group7",
            "group8",
            "group9",
            "group10",
            "group11",
            "group12",
            "group13",
            "group14",
            "group15",
            "group16",
            "group17",
            "group18",
            "group19",
            "group20",
        ]
        # other opaque gates, which are not handled as barrier
        # ["RZZ", "Rxxyyzz", "Rxxyyzz_zphase", "cu", "cp", "rccx", "rc3x", "c3sqrtx"]

        args = list(args)

        if opstr in treat_as_barrier:
            params = [f"{par}" for par in pars]
        else:
            params = [f"({par})/pi" for par in pars]
        if opstr in self.gate_dict:
            op: Dict[str, Any] = {}
            if opstr in treat_as_barrier:
                op["type"] = "Barrier"
                param_sorted = ",".join(params)

                op["data"] = f"{opstr}({param_sorted})"

                op["signature"] = [arg[0] for arg in args]
            else:
                gdef = self.gate_dict[opstr]
                op["type"] = "CustomGate"
                box = {
                    "type": "CustomGate",
                    "id": str(uuid.uuid4()),
                    "gate": gdef,
                }
                box["params"] = params
                op["box"] = box
                params = []  # to stop duplication in to op
        else:
            try:
                optype = _all_string_maps[opstr]
            except KeyError as e:
                raise QASMParseError(
                    "Cannot parse gate of type: {}".format(opstr), optoken.line
                ) from e
            op = {"type": optype}
            if params:
                op["params"] = params
            # Operations needing special handling:
            if optype.startswith("Cn"):
                # n-controlled rotations have variable signature
                op["n_qb"] = len(args)
            elif optype == "Barrier":
                op["signature"] = ["Q"] * len(args)

        for arg in zip(*self.unroll_all_args(args)):
            yield {"args": list(arg), "op": op}

    def gatecall(self, tree: List) -> Iterable[CommandDict]:
        return self.mixedcall(tree)

    def exp_args(self, tree: Iterable[Token]) -> Iterable[Reg]:
        for arg in tree:
            if arg.type == "ARG":
                yield self._get_reg(arg.value)
            else:
                raise QASMParseError(
                    "Non register arguments not supported for extern call.", arg.line
                )

    def _logic_exp(self, tree: List, opstr: str) -> LogicExp:
        args, line = self._get_logic_args(tree)
        openum: Union[Type[BitWiseOp], Type[RegWiseOp]]
        if opstr in _BITOPS and opstr not in _REGOPS:
            openum = BitWiseOp
        elif opstr in _REGOPS and opstr not in _BITOPS:
            openum = RegWiseOp
        elif all(isinstance(arg, int) for arg in args):
            openum = RegWiseOp
        elif all(isinstance(arg, (Bit, BitLogicExp, int)) for arg in args):
            if all(arg in (0, 1) for arg in args if isinstance(arg, int)):
                openum = BitWiseOp
            else:
                raise QASMParseError(
                    "Bits can only be operated with (0, 1) literals."
                    f" Incomaptible arguments {args}",
                    line,
                )
        else:
            openum = RegWiseOp
        if openum is BitWiseOp and opstr in ("+", "-"):
            op: Union[BitWiseOp, RegWiseOp] = BitWiseOp.XOR
        else:
            op = openum(opstr)
        return LogicExp.factory(op)(*args)  # type: ignore

    def _get_logic_args(
        self, tree: Sequence[Union[Token, LogicExp]]
    ) -> Tuple[List[Union[LogicExp, Bit, BitRegister, int]], Optional[int]]:
        args: List[Union[LogicExp, Bit, BitRegister, int]] = []
        line = None
        for tok in tree:
            if isinstance(tok, LogicExp):
                args.append(tok)
            elif isinstance(tok, Token):
                line = tok.line
                if tok.type == "INT":
                    args.append(int(tok.value))
                elif tok.type == "IARG":
                    args.append(Bit(*_extract_reg(tok)))
                elif tok.type == "ARG":
                    args.append(BitRegister(tok.value, self.c_registers[tok.value]))
                else:
                    raise QASMParseError(f"Could not pass argument {tok}")
            else:
                raise QASMParseError(f"Could not pass argument {tok}")
        return args, line

    par_add = _bin_par_exp("+")
    par_sub = _bin_par_exp("-")
    par_mul = _bin_par_exp("*")
    par_div = _bin_par_exp("/")
    par_pow = _bin_par_exp("**")

    par_neg = _un_par_exp("-")

    sqrt = _un_call_exp("sqrt")
    sin = _un_call_exp("sin")
    cos = _un_call_exp("cos")
    tan = _un_call_exp("tan")
    ln = _un_call_exp("ln")

    b_and = lambda self, tree: self._logic_exp(tree, "&")
    b_not = lambda self, tree: self._logic_exp(tree, "~")
    b_or = lambda self, tree: self._logic_exp(tree, "|")
    xor = lambda self, tree: self._logic_exp(tree, "^")
    lshift = lambda self, tree: self._logic_exp(tree, "<<")
    rshift = lambda self, tree: self._logic_exp(tree, ">>")
    add = lambda self, tree: self._logic_exp(tree, "+")
    sub = lambda self, tree: self._logic_exp(tree, "-")
    mul = lambda self, tree: self._logic_exp(tree, "*")
    div = lambda self, tree: self._logic_exp(tree, "/")
    ipow = lambda self, tree: self._logic_exp(tree, "**")

    def neg(self, tree: List[Union[Token, LogicExp]]) -> RegNeg:
        return RegNeg(self._get_logic_args(tree)[0][0])

    def cond(self, tree: List[Token]) -> PredicateExp:
        if tree[1].type == "IARG":
            arg = Bit(*_extract_reg(tree[1]))
        else:
            arg = BitRegister(tree[1].value, self.c_registers[tree[1].value])  # type: ignore

        op_enum = BitWiseOp if isinstance(arg, Bit) else RegWiseOp
        comp = cast(
            Type[PredicateExp],
            LogicExp.factory(
                cast(
                    Union[BitWiseOp, RegWiseOp],
                    op_enum(str(tree[2])),
                )
            ),
        )
        return comp(arg, int(tree[3].value))

    def ifc(self, tree: Sequence) -> Iterable[CommandDict]:
        condition = cast(PredicateExp, tree[0])

        var, val = condition.args
        condition_bits = []

        if isinstance(var, Bit):
            assert condition.op in (BitWiseOp.EQ, BitWiseOp.NEQ)
            assert val in (0, 1)
            condition_bits = [var.to_list()]

        else:
            assert isinstance(var, BitRegister)
            reg_bits = next(self.unroll_all_args([var.name]))
            if isinstance(condition, RegEq):
                # special case for base qasm
                condition_bits = reg_bits
            else:
                pred_val = cast(int, val)
                minval = 0
                maxval = (1 << 32) - 1
                if condition.op == RegWiseOp.LT:
                    maxval = pred_val - 1
                elif condition.op == RegWiseOp.GT:
                    minval = pred_val + 1
                if condition.op in (RegWiseOp.LEQ, RegWiseOp.EQ, RegWiseOp.NEQ):
                    maxval = pred_val
                if condition.op in (RegWiseOp.GEQ, RegWiseOp.EQ, RegWiseOp.NEQ):
                    minval = pred_val

                condition_bit = self._fresh_temp_bit()
                yield {
                    "args": reg_bits + [condition_bit],
                    "op": {
                        "classical": {
                            "lower": minval,
                            "n_i": len(reg_bits),
                            "upper": maxval,
                        },
                        "type": "RangePredicate",
                    },
                }
                condition_bits = [condition_bit]
                val = int(condition.op != RegWiseOp.NEQ)

        for com in filter(lambda x: x is not None and x is not Discard, tree[1]):
            com["args"] = condition_bits + com["args"]
            com["op"] = {
                "conditional": {
                    "op": com["op"],
                    "value": val,
                    "width": len(condition_bits),
                },
                "type": "Conditional",
            }

            yield com

    def cop(self, tree: Sequence[Iterable[CommandDict]]) -> Iterable[CommandDict]:
        return tree[0]

    def _calc_exp_io(
        self, exp: LogicExp, out_args: List
    ) -> Tuple[List[List], Dict[str, Any]]:
        all_inps: Set[Tuple[str, int]] = set(
            map(
                lambda b: (b.reg_name, b.index[0]),
                chain.from_iterable(
                    (
                        [inp] if isinstance(inp, Bit) else iter(inp)
                        for inp in exp.all_inputs()
                    )
                ),
            )
        )
        outs = (_hashable_uid(arg) for arg in out_args)
        o = []
        io = []
        for out in outs:
            if out in all_inps:
                all_inps.remove(out)
                io.append(out)
            else:
                o.append(out)

        exp_args = list(
            map(lambda x: [x[0], [x[1]]], chain.from_iterable((all_inps, io, o)))
        )
        numbers_dict = {
            "n_i": len(all_inps),
            "n_io": len(io),
            "n_o": len(o),
        }
        return exp_args, numbers_dict

    def _cexpbox_dict(self, exp: LogicExp, args: List[List]) -> CommandDict:
        box = {
            "exp": exp.to_dict(),
            "id": str(uuid.uuid4()),
            "type": "ClassicalExpBox",
        }
        args, numbers = self._calc_exp_io(exp, args)
        box.update(numbers)
        return {
            "args": args,
            "op": {
                "box": box,
                "type": "ClassicalExpBox",
            },
        }

    def assign(self, tree: List) -> Iterable[CommandDict]:
        child_iter = iter(tree)
        out_args = list(next(child_iter))
        args_uids = list(self.unroll_all_args(out_args))

        exp_tree = next(child_iter)

        exp: Union[str, List, LogicExp, int] = ""
        line = None
        if isinstance(exp_tree, Token):
            if exp_tree.type == "INT":
                exp = int(exp_tree.value)
            elif exp_tree.type in ("ARG", "IARG"):
                exp = self._get_arg(exp_tree)
            line = exp_tree.line
        elif isinstance(exp_tree, Generator):
            # assume to be extern (wasm) call
            chained_uids = list(chain.from_iterable(args_uids))
            com = next(exp_tree)
            com["args"].pop()  # remove the wasmstate from the args
            com["args"] += chained_uids
            com["args"].append(["_w", [0]])
            com["op"]["wasm"]["n"] += len(chained_uids)
            com["op"]["wasm"]["width_o_parameter"] = [
                self.c_registers[reg] for reg in out_args
            ]

            yield com
            return
        else:
            exp = exp_tree

        assert len(out_args) == 1
        out_arg = out_args[0]
        args = args_uids[0]
        if isinstance(out_arg, List):
            if isinstance(exp, LogicExp):
                yield self._cexpbox_dict(exp, args)
            elif isinstance(exp, (int, bool)):
                assert exp in (0, 1, True, False)
                yield {
                    "args": args,
                    "op": {"classical": {"values": [bool(exp)]}, "type": "SetBits"},
                }
            elif isinstance(exp, List):
                yield {
                    "args": [exp] + args,
                    "op": {"classical": {"n_i": 1}, "type": "CopyBits"},
                }
            else:
                raise QASMParseError(f"Unexpected expression in assignment {exp}", line)
        else:
            reg = out_arg
            if isinstance(exp, RegLogicExp):
                yield self._cexpbox_dict(exp, args)
            elif isinstance(exp, BitLogicExp):
                yield self._cexpbox_dict(exp, args[:1])
            elif isinstance(exp, int):
                yield {
                    "args": args,
                    "op": {
                        "classical": {
                            "values": int_to_bools(exp, self.c_registers[reg])
                        },
                        "type": "SetBits",
                    },
                }

            elif isinstance(exp, str):
                width = min(self.c_registers[exp], len(args))
                yield {
                    "args": [[exp, [i]] for i in range(width)] + args[:width],
                    "op": {"classical": {"n_i": width}, "type": "CopyBits"},
                }

            else:
                raise QASMParseError(f"Unexpected expression in assignment {exp}", line)

    def extern(self, tree: List[Any]) -> Type[Discard]:
        # TODO parse extern defs
        return Discard

    def ccall(self, tree: List) -> Iterable[CommandDict]:
        return self.cce_call(tree)

    def cce_call(self, tree: List) -> Iterable[CommandDict]:
        nam = tree[0].value
        params = list(tree[1])
        if self.wasm is None:
            raise QASMParseError(
                "Cannot include extern calls without a wasm module specified.",
                tree[0].line,
            )
        n_i_vec = [self.c_registers[reg] for reg in params]

        wasm_args = list(chain.from_iterable(self.unroll_all_args(params)))

        wasm_args.append(["_w", [0]])

        yield {
            "args": wasm_args,
            "op": {
                "type": "WASM",
                "wasm": {
                    "func_name": nam,
                    "ww_n": 1,
                    "n": sum(n_i_vec),
                    "width_i_parameter": n_i_vec,
                    "width_o_parameter": [],  # this will be set in the assign function
                    "wasm_file_uid": str(self.wasm),
                },
            },
        }

    def transform(self, tree: Tree) -> Dict[str, Any]:
        self._reset_context()
        return cast(Dict[str, Any], super().transform(tree))

    def gdef(self, tree: List) -> None:
        child_iter = iter(tree)

        gate = next(child_iter).value

        next_tree = next(child_iter)
        symbols, args = [], []
        if isinstance(next_tree, ParsMap):
            symbols = list(next_tree.pars)
            args = list(next(child_iter))
        else:
            args = list(next_tree)

        symbol_map = {sym: sym * pi for sym in map(Symbol, symbols)}
        rename_map = {Qubit.from_list(qb): Qubit("q", i) for i, qb in enumerate(args)}

        new = CircuitTransformer()
        circ_dict = new.prog(child_iter)

        circ_dict["qubits"] = args
        gate_circ = Circuit.from_dict(circ_dict)

        # check to see whether gate definition was generated by pytket converter
        # if true, add op as pytket Op
        existing_op: bool = False
        if gate in NOPARAM_EXTRA_COMMANDS:
            qubit_args = [
                Qubit(gate + "q" + str(index), 0) for index in list(range(len(args)))
            ]
            comparison_circ = _get_gate_circuit(
                NOPARAM_EXTRA_COMMANDS[gate], qubit_args
            )
            if circuit_to_qasm_str(comparison_circ) == circuit_to_qasm_str(gate_circ):
                existing_op = True
        elif gate in PARAM_EXTRA_COMMANDS:
            qubit_args = [
                Qubit(gate + "q" + str(index), 0) for index in list(range(len(args)))
            ]
            comparison_circ = _get_gate_circuit(
                PARAM_EXTRA_COMMANDS[gate],
                qubit_args,
                [
                    Symbol("param" + str(index) + "/pi") for index in range(len(symbols))  # type: ignore
                ],
            )
            # checks that each command has same string
            existing_op = all(
                str(g) == str(c)
                for g, c in zip(
                    gate_circ.get_commands(), comparison_circ.get_commands()
                )
            )
        if not existing_op:
            gate_circ.symbol_substitution(symbol_map)
            gate_circ.rename_units(cast(Dict[UnitID, UnitID], rename_map))

            self.gate_dict[gate] = {
                "definition": gate_circ.to_dict(),
                "args": symbols,
                "name": gate,
            }

    opaq = gdef

    def oqasm(self, tree: List) -> Type[Discard]:
        return Discard

    def incl(self, tree: List[Token]) -> None:
        self.include = str(tree[0].value).split(".")[0]
        self.gate_dict.update(_load_include_module(self.include, True, False))

    def prog(self, tree: Iterable) -> Dict[str, Any]:
        outdict: Dict[str, Any] = {
            "commands": list(
                chain.from_iterable(
                    filter(lambda x: x is not None and x is not Discard, tree)
                )
            )
        }
        if self.return_gate_dict:
            return self.gate_dict
        outdict["qubits"] = [
            [reg, [i]] for reg, size in self.q_registers.items() for i in range(size)
        ]
        outdict["bits"] = [
            [reg, [i]] for reg, size in self.c_registers.items() for i in range(size)
        ]
        outdict["implicit_permutation"] = [[q, q] for q in outdict["qubits"]]
        outdict["phase"] = "0.0"
        self._reset_context()
        return outdict


parser = Lark(
    grammar,
    start="prog",
    debug=False,
    parser="lalr",
    cache=True,
    transformer=CircuitTransformer(),
)


def circuit_from_qasm(
    input_file: Union[str, "os.PathLike[Any]"], encoding: str = "utf-8"
) -> Circuit:
    """A method to generate a tket Circuit from a qasm file"""
    ext = os.path.splitext(input_file)[-1]
    if ext != ".qasm":
        raise TypeError("Can only convert .qasm files")
    with open(input_file, "r", encoding=encoding) as f:
        try:
            circ = circuit_from_qasm_io(f)
        except QASMParseError as e:
            raise QASMParseError(e.msg, e.line, str(input_file))
    return circ


def circuit_from_qasm_str(qasm_str: str) -> Circuit:
    """A method to generate a tket Circuit from a qasm str"""
    cast(CircuitTransformer, parser.options.transformer)._reset_context(
        reset_wasm=False
    )
    return Circuit.from_dict(parser.parse(qasm_str))  # type: ignore[arg-type]


def circuit_from_qasm_io(stream_in: TextIO) -> Circuit:
    """A method to generate a tket Circuit from a qasm text stream"""
    return circuit_from_qasm_str(stream_in.read())


def circuit_from_qasm_wasm(
    input_file: Union[str, "os.PathLike[Any]"],
    wasm_file: Union[str, "os.PathLike[Any]"],
    encoding: str = "utf-8",
) -> Circuit:
    """A method to generate a tket Circuit from a qasm str and external WASM module."""
    wasm_module = WasmFileHandler(str(wasm_file))
    cast(CircuitTransformer, parser.options.transformer).wasm = wasm_module
    return circuit_from_qasm(input_file, encoding=encoding)


def circuit_to_qasm(circ: Circuit, output_file: str, header: str = "qelib1") -> None:
    """A method to generate a qasm file from a tket Circuit.
    Note that this will not account for implicit qubit permutations in the Circuit."""
    with open(output_file, "w") as out:
        circuit_to_qasm_io(circ, out, header=header)


def _filtered_qasm_str(qasm: str) -> str:
    # remove any c registers starting with _TEMP_BIT_NAME
    # that are not being used somewhere else
    lines = qasm.split("\n")
    def_matcher = re.compile(r"creg ({}\_*\d*)\[\d+\]".format(_TEMP_BIT_NAME))
    arg_matcher = re.compile(r"({}\_*\d*)\[\d+\]".format(_TEMP_BIT_NAME))
    unused_regs = dict()
    for i, line in enumerate(lines):
        if reg := def_matcher.match(line):
            # Mark a reg temporarily as unused
            unused_regs[reg.group(1)] = i
        elif args := arg_matcher.findall(line):
            # If the line contains scratch bits that are used as arguments
            # mark these regs as used
            for arg in args:
                if arg in unused_regs:
                    unused_regs.pop(arg)
    # remove unused reg defs
    redundant_lines = sorted(unused_regs.values(), reverse=True)
    for line_index in redundant_lines:
        del lines[line_index]
    return "\n".join(lines)


def circuit_to_qasm_str(circ: Circuit, header: str = "qelib1") -> str:
    """A method to generate a qasm str from a tket Circuit.
    Note that this will not account for implicit qubit permutations in the Circuit."""
    buffer = io.StringIO()
    circuit_to_qasm_io(circ, buffer, header=header)
    return buffer.getvalue()


TypeReg = TypeVar("TypeReg", BitRegister, QubitRegister)


def _retrieve_registers(
    units: List[UnitID], reg_type: Type[TypeReg]
) -> Dict[str, TypeReg]:
    if any(len(unit.index) != 1 for unit in units):
        raise NotImplementedError("OPENQASM registers must use a single index")
    maxunits = map(lambda x: max(x[1]), groupby(units, key=lambda un: un.reg_name))
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
    params = (
        op.params
        if (optype in _tk_to_qasm_params) or (optype in _tk_to_qasm_extra_params)
        else None
    )
    if optype == OpType.TK1:
        # convert to U3
        optype = OpType.U3
        params = [op.params[1], op.params[0] - 0.5, op.params[2] + 0.5]
    elif optype == OpType.CustomGate:
        params = op.params
    return optype, params  # type: ignore


def _get_gate_circuit(
    optype: OpType, qubits: List[Qubit], symbols: Optional[List[Symbol]] = None
) -> Circuit:
    # create Circuit for constructing qasm from
    unitids = cast(List[UnitID], qubits)
    gate_circ = Circuit()
    for q in qubits:
        gate_circ.add_qubit(q)
    if symbols:
        exprs = [symbol.as_expr() for symbol in symbols]  # type: ignore
        gate_circ.add_gate(optype, exprs, unitids)
    else:
        gate_circ.add_gate(optype, unitids)
    auto_rebase_pass({OpType.CX, OpType.U3}).apply(gate_circ)
    RemoveRedundancies().apply(gate_circ)

    return gate_circ


def _write_gate_definition(
    buffer: TextIO,
    n_qubits: int,
    opstr: str,
    optype: OpType,
    header: str,
    include_gate_defs: Set[str],
    params: Optional[List[float]] = None,
) -> None:
    # start writing to stream
    buffer.write("gate " + opstr + " ")
    symbols: Optional[List[Symbol]] = None
    if params:
        # need to add parameters to gate definition
        buffer.write("(")
        symbols = [Symbol("param" + str(index) + "/pi") for index in range(len(params))]  # type: ignore
        symbols_header = [Symbol("param" + str(index)) for index in range(len(params))]  # type: ignore
        for symbol in symbols_header[:-1]:
            buffer.write(symbol.name + ", ")
        buffer.write(symbols_header[-1].name + ") ")

    # add qubits to gate definition
    qubit_args = [Qubit(opstr + "q" + str(index)) for index in list(range(n_qubits))]
    for qb in qubit_args[:-1]:
        buffer.write(str(qb) + ",")
    buffer.write(str(qubit_args[-1]) + " {\n")
    # get rebased circuit for constructing qasm
    gate_circ = _get_gate_circuit(optype, qubit_args, symbols)
    # write circuit to qasm
    circuit_to_qasm_io(gate_circ, buffer, header, include_gate_defs)
    buffer.write("}\n")


def hqs_header(header: str) -> bool:
    return header in ["hqslib1", "hqslib1_dev"]


def circuit_to_qasm_io(
    circ: Circuit,
    stream_out: TextIO,
    header: str = "qelib1",
    include_gate_defs: Optional[Set[str]] = None,
) -> None:
    """A method to generate a qasm text stream from a tket Circuit.
    Note that this will not account for implicit qubit permutations in the Circuit."""
    # Write to a buffer since the output qasm might need additional filtering.
    # e.g. remove unused tket scratch bits.
    buffer = io.StringIO()
    if any(
        circ.n_gates_of_type(typ)
        for typ in (
            OpType.RangePredicate,
            OpType.MultiBit,
            OpType.ExplicitPredicate,
            OpType.ExplicitModifier,
            OpType.SetBits,
            OpType.CopyBits,
            OpType.ClassicalExpBox,
        )
    ) and (not hqs_header(header)):
        raise QASMUnsupportedError(
            "Complex classical gates not supported with qelib1: try converting with "
            "`header=hqslib1`"
        )
    include_module_gates = {"measure", "reset", "barrier"}
    include_module_gates.update(_load_include_module(header, False, True).keys())
    if include_gate_defs is None:
        include_gate_defs = include_module_gates
        include_gate_defs.update(NOPARAM_EXTRA_COMMANDS.keys())
        include_gate_defs.update(PARAM_EXTRA_COMMANDS.keys())
        buffer.write('OPENQASM 2.0;\ninclude "{}.inc";\n\n'.format(header))

        qregs = _retrieve_registers(cast(list[UnitID], circ.qubits), QubitRegister)
        cregs = _retrieve_registers(cast(list[UnitID], circ.bits), BitRegister)
        for reg in qregs.values():
            buffer.write(f"qreg {reg.name}[{reg.size}];\n")
        for bit_reg in cregs.values():
            buffer.write(f"creg {bit_reg.name}[{bit_reg.size}];\n")
    else:
        # gate definition, no header necessary for file
        cregs = {}
        qregs = {}

    added_gate_definitions: Set[str] = set()
    range_preds = dict()
    for command in circ:
        assert isinstance(command, Command)
        checked_op = True
        op = command.op
        args = command.args
        optype, params = _get_optype_and_params(op)
        if optype == OpType.RangePredicate:
            assert isinstance(op, RangePredicateOp)
            range_preds[args[-1]] = (op, args)
            # attach predicate to bit,
            # subsequent conditional will handle it
            continue
        if optype == OpType.Conditional:
            assert isinstance(op, Conditional)
            bits = args[: op.width]
            control_bit = bits[0]
            if control_bit in range_preds:
                # write range predicate in condition
                range_op, range_args = range_preds[control_bit]
                comparator, value = _parse_range(range_op.lower, range_op.upper)
                if op.value == 0 and comparator == "==":
                    comparator = "!="
                if (not hqs_header(header)) and comparator != "==":
                    raise QASMUnsupportedError(
                        "OpenQASM conditions must be on a register's fixed value."
                    )
                bits = range_args[:-1]
                variable: Union[
                    str,
                    UnitID,
                ] = range_args[0].reg_name
            else:
                comparator = "=="
                value = op.value
                if op.width == 1 and hqs_header(header):
                    variable = control_bit
                else:
                    variable = control_bit.reg_name
                    if hqs_header(header) and bits != cregs[variable].to_list():
                        raise QASMUnsupportedError(
                            "hqslib1 QASM conditions must be an entire classical "
                            "register or a single bit"
                        )
            if not hqs_header(header):
                assert isinstance(variable, str)
                if op.width != cregs[variable].size:
                    raise QASMUnsupportedError(
                        "OpenQASM conditions must be an entire classical register"
                    )
                if bits != cregs[variable].to_list():
                    raise QASMUnsupportedError(
                        "OpenQASM conditions must be a single classical register"
                    )
            args = args[op.width :]
            op = op.op
            optype, params = _get_optype_and_params(op)
            if optype != OpType.Phase:
                buffer.write(f"if({variable}{comparator}{value}) ")
        if optype == OpType.Phase:
            # global phase is ignored in QASM
            continue
        if optype == OpType.SetBits:
            assert isinstance(op, SetBitsOp)
            creg_name = args[0].reg_name
            bits, vals = zip(*sorted(zip(args, op.values)))  # type: ignore

            # check if whole register can be set at once
            if bits == tuple(cregs[creg_name].to_list()):
                value = int("".join(map(str, map(int, vals[::-1]))), 2)
                buffer.write(f"{creg_name} = {value};\n")
            else:
                for bit, value in zip(bits, vals):
                    buffer.write(f"{bit} = {int(value)};\n")
            continue
        if optype == OpType.CopyBits:
            assert isinstance(op, CopyBitsOp)
            l_args = args[op.n_inputs :]
            r_args = args[: op.n_inputs]
            l_name = l_args[0].reg_name
            r_name = r_args[0].reg_name

            # check if whole register can be set at once
            if l_args == cregs[l_name].to_list() and r_args == cregs[r_name].to_list():
                buffer.write(f"{l_name} = {r_name};\n")
            else:
                for bit_l, bit_r in zip(l_args, r_args):
                    buffer.write(f"{bit_l} = {bit_r};\n")
            continue
        if optype == OpType.MultiBit:
            assert isinstance(op, MultiBitOp)
            l_args = args[op.n_inputs :]
            op = op.basic_op
            optype = op.type
            registers_involved = [arg.reg_name for arg in args[:2]]
            if len(args) > 2 and args[2].reg_name not in registers_involved:
                # there is a distinct output register
                registers_involved.append(args[2].reg_name)
            args = [cregs[name] for name in registers_involved]  # type: ignore
        if optype in (
            OpType.ExplicitPredicate,
            OpType.ExplicitModifier,
        ):
            # &, ^ and | gates
            opstr = str(op)
            if opstr not in _classical_gatestr_map:
                raise QASMUnsupportedError(f"Classical gate {opstr} not supported.")
            buffer.write(
                f"{args[-1]} = {args[0]} {_classical_gatestr_map[opstr]} {args[1]};\n"
            )
            continue

        if optype == OpType.ClassicalExpBox:
            assert isinstance(op, ClassicalExpBox)
            out_args: list[UnitID] = args[op.get_n_i() :]
            if len(out_args) == 1:
                buffer.write(f"{out_args[0]} = {str(op.get_exp())};\n")
            elif (
                out_args
                == cregs[out_args[0].reg_name].to_list()[: op.get_n_io() + op.get_n_o()]
            ):
                buffer.write(f"{out_args[0].reg_name} = {str(op.get_exp())};\n")
            else:
                raise QASMUnsupportedError(
                    f"ClassicalExpBox only supported"
                    " for writing to a single bit or whole registers."
                )
            continue
        if optype == OpType.WASM:
            assert isinstance(op, WASMOp)
            inputs: List[str] = []
            outputs: List[str] = []
            for reglist, sizes in [
                (inputs, op.input_widths),
                (outputs, op.output_widths),
            ]:
                for in_width in sizes:
                    bits = args[:in_width]
                    args = args[in_width:]
                    regname = bits[0].reg_name
                    if bits != list(cregs[regname]):  # type: ignore
                        QASMUnsupportedError("WASM ops must act on entire registers.")
                    reglist.append(regname)
            if outputs:
                buffer.write(f"{', '.join(outputs)} = ")
            buffer.write(f"{op.func_name}({', '.join(inputs)});\n")
            continue
        if optype == OpType.CustomGate:
            assert isinstance(op, CustomGate)
            if op.gate.name not in include_gate_defs:
                # unroll custom gate
                gate_circ = op.get_circuit()

                if gate_circ.n_gates == 0:
                    raise QASMUnsupportedError(
                        f"CustomGate {op.gate.name} has empty definition."
                        " Empty CustomGates and opaque gates are not supported."
                    )
                gate_circ.rename_units(dict(zip(gate_circ.qubits, args)))
                gate_circ.symbol_substitution(dict(zip(op.gate.args, op.params)))

                circuit_to_qasm_io(gate_circ, buffer, header, include_gate_defs)
                continue
            else:
                opstr = op.gate.name
        elif hqs_header(header) and optype == OpType.ZZPhase:
            # special handling for zzphase
            opstr = "RZZ"
            param = op.params[0]
            # as op.params returns reduced parameters, we can assume
            # that 0 <= param < 4
            if param > 1:
                # first get in to 0 <= param < 2 range
                param = Decimal(str(param)) % Decimal("2")  # type: ignore
                # then flip 1 <= param < 2  range into
                # -1 <= param < 0
                if param > 1:
                    param = -2 + param
            params = [param]  # type: ignore
        elif optype == OpType.Barrier and header == "hqslib1_dev":
            assert isinstance(op, MetaOp)
            if op.data == "":
                opstr = _tk_to_qasm_noparams[optype]
            else:
                opstr = op.data
                checked_op = False

        elif (
            optype in _tk_to_qasm_noparams
            and _tk_to_qasm_noparams[optype] in include_module_gates
        ):
            opstr = _tk_to_qasm_noparams[optype]
        elif (
            optype in _tk_to_qasm_params
            and _tk_to_qasm_params[optype] in include_module_gates
        ):
            opstr = _tk_to_qasm_params[optype]
        elif optype in _tk_to_qasm_extra_noparams:
            opstr = _tk_to_qasm_extra_noparams[optype]
            n_added = len(added_gate_definitions)
            added_gate_definitions.add(opstr)
            if len(added_gate_definitions) != n_added:
                _write_gate_definition(
                    buffer, op.n_qubits, opstr, optype, header, include_gate_defs
                )
        elif optype in _tk_to_qasm_extra_params:
            opstr = _tk_to_qasm_extra_params[optype]
            n_added = len(added_gate_definitions)
            added_gate_definitions.add(opstr)
            if len(added_gate_definitions) != n_added:
                _write_gate_definition(
                    buffer,
                    op.n_qubits,
                    opstr,
                    optype,
                    header,
                    include_gate_defs,
                    params,
                )
        else:
            raise QASMUnsupportedError(
                "Cannot print command of type: {}".format(op.get_name())
            )
        if checked_op and opstr not in include_gate_defs:
            raise QASMUnsupportedError(
                "Gate of type {} is not supported in conversion.".format(opstr)
            )
        buffer.write(opstr)
        if params is not None:
            buffer.write("(")
            for i in range(len(params)):
                reduced = True
                try:
                    p = float(params[i])
                except TypeError:
                    reduced = False
                    p = params[i]
                if i < len(params) - 1:
                    if reduced:
                        buffer.write("{}*pi,".format(p))
                    else:
                        buffer.write("({})*pi,".format(p))

                else:
                    if reduced:
                        buffer.write("{}*pi)".format(p))
                    else:
                        buffer.write("({})*pi)".format(p))
        buffer.write(" ")
        if optype == OpType.Measure:  # assume written to only 1 bit
            buffer.write(
                "{q} -> {c};\n".format(q=args[0].__repr__(), c=args[1].__repr__())
            )
        else:
            for i in range(len(args)):
                buffer.write(args[i].__repr__())
                if i < len(args) - 1:
                    buffer.write(",")
                else:
                    buffer.write(";\n")

    stream_out.write(_filtered_qasm_str(buffer.getvalue()))
