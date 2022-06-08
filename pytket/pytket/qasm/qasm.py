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

from dataclasses import dataclass
import io
import os
import re
import uuid

# TODO: Output custom gates
from importlib import import_module
from itertools import chain, groupby
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

from pytket._tket.circuit import _TEMP_BIT_NAME  # type: ignore
from pytket.circuit import (  # type: ignore
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
    ConstPredicate,
    LogicExp,
    RegEq,
    RegLogicExp,
    RegNeg,
    RegWiseOp,
)
from pytket.qasm.grammar import grammar
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
    "rzz": OpType.ZZPhase,
    "Rz": OpType.Rz,
    "U1q": OpType.PhasedX,
    "crz": OpType.CRz,
    "crx": OpType.CRx,
    "cry": OpType.CRy,
    "cu1": OpType.CU1,
    "cu3": OpType.CU3,
}

_tk_to_qasm_noparams = dict(((item[1], item[0]) for item in NOPARAM_COMMANDS.items()))
_tk_to_qasm_noparams[OpType.CX] = "cx"  # prefer "cx" to "CX"
_tk_to_qasm_params = dict(((item[1], item[0]) for item in PARAM_COMMANDS.items()))
_tk_to_qasm_params[OpType.U3] = "u3"  # prefer "u3" to "U"
_tk_to_qasm_params[OpType.Rz] = "rz"  # prefer "rz" to "Rz"

_classical_gatestr_map = {"AND": "&", "OR": "|", "XOR": "^"}


_all_known_gates = set(NOPARAM_COMMANDS.keys()).union(PARAM_COMMANDS.keys())
_all_string_maps = {
    key: val.name
    for key, val in chain(PARAM_COMMANDS.items(), NOPARAM_COMMANDS.items())
}

unit_regex = re.compile(r"([a-z][a-zA-Z0-9_]*)\[([\d]+)\]")


def _extract_reg(var: Token) -> Tuple[str, int]:
    match = unit_regex.match(var.value)
    if match is None:
        raise QASMParseError(f"Not a valid (qu)bit identifier: {var.value}", var.line)
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

    def _reset_context(self) -> None:
        self.q_registers = {}
        self.c_registers = {}
        self.gate_dict = {}
        self.include = ""
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

    def unroll_all_args(self, args: Iterable[Arg]) -> Iterator[List[UnitID]]:
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
        params = [f"({par})/pi" for par in pars]
        if opstr in self.gate_dict:
            gdef = self.gate_dict[opstr]
            op: Dict[str, Any] = {"type": "CustomGate"}
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
            if optype.startswith("Cn"):
                # n-controlled rotation are only gates supported
                # via this function without fixed signature
                args = list(args)
                op["n_qb"] = len(args)

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
        return args, line  # type: ignore

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

    def cond(self, tree: List[Token]) -> ConstPredicate:
        if tree[1].type == "IARG":
            arg = Bit(*_extract_reg(tree[1]))
        else:
            arg = BitRegister(tree[1].value, self.c_registers[tree[1].value])

        op_enum = BitWiseOp if isinstance(arg, Bit) else RegWiseOp
        comp = cast(
            Type[ConstPredicate],
            LogicExp.factory(
                cast(
                    Union[BitWiseOp, RegWiseOp],
                    op_enum(str(tree[2])),
                )
            ),
        )
        return comp(arg, int(tree[3].value))

    def ifc(self, tree: Sequence) -> Iterable[CommandDict]:
        condition = cast(ConstPredicate, tree[0])

        var, val = condition.args
        condition_bits = []

        if isinstance(var, Bit):
            assert condition.op in (BitWiseOp.EQ, BitWiseOp.NEQ)
            assert val in (0, 1)
            condition_bits = [cast(Bit, var).to_list()]

        else:
            assert isinstance(var, BitRegister)
            reg_bits = next(self.unroll_all_args([cast(BitRegister, var).name]))
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
            com["args"] += chained_uids
            com["op"]["wasm"]["n"] += len(chained_uids)
            com["op"]["wasm"]["no_vec"] = [self.c_registers[reg] for reg in out_args]
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
                # right_reg = cast(BitRegister, exp)
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

        yield {
            "args": list(chain.from_iterable(self.unroll_all_args(params))),
            "op": {
                "type": "WASM",
                "wasm": {
                    "func_name": nam,
                    "n": sum(n_i_vec),
                    "ni_vec": n_i_vec,
                    "no_vec": [],
                    "wasm_uid": str(self.wasm),
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
        def_circ = Circuit.from_dict(circ_dict)

        def_circ.symbol_substitution(symbol_map)
        def_circ.rename_units(rename_map)

        self.gate_dict[gate] = {
            "definition": def_circ.to_dict(),
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
    return Circuit.from_dict(parser.parse(qasm_str))


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
    elif optype == OpType.CustomGate:
        params = op.params
    return (optype, params)


def circuit_to_qasm_io(
    circ: Circuit,
    stream_out: TextIO,
    header: str = "qelib1",
    include_gate_defs: Optional[Set[str]] = None,
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
                OpType.CopyBits,
                OpType.ClassicalExpBox,
            )
        )
        and header != "hqslib1"
    ):
        raise QASMUnsupportedError(
            "Complex classical gates only supported with hqslib1."
        )
    if include_gate_defs is None:
        include_gate_defs = {"measure", "reset", "barrier"}
        include_gate_defs.update(_load_include_module(header, False, True).keys())
        stream_out.write('OPENQASM 2.0;\ninclude "{}.inc";\n\n'.format(header))
        qregs = _retrieve_registers(circ.qubits, QubitRegister)
        cregs = _retrieve_registers(circ.bits, BitRegister)

        for reg in qregs.values():
            stream_out.write(f"qreg {reg.name}[{reg.size}];\n")
        for reg in cregs.values():
            stream_out.write(f"creg {reg.name}[{reg.size}];\n")
    else:
        # gate definition, no header necessary for file
        cregs = {}
        qregs = {}

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
        if optype == OpType.Conditional:
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
            else:
                comparator = "=="
                value = op.value
                if op.width == 1 and header == "hqslib1":
                    variable = control_bit
                else:
                    variable = control_bit.reg_name
            if header != "hqslib1":
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
        if optype == OpType.CopyBits:
            l_args = args[op.n_inputs :]
            r_args = args[: op.n_inputs]
            l_name = l_args[0].reg_name
            r_name = r_args[0].reg_name

            # check if whole register can be set at once
            if l_args == list(cregs[l_name]) or r_args == list(cregs[r_name]):
                stream_out.write(f"{l_name} = {r_name};\n")
            else:
                for bit_l, bit_r in zip(l_args, r_args):
                    stream_out.write(f"{bit_l} = {bit_r};\n")
            continue
        if optype == OpType.MultiBit:
            op = op.basic_op
            optype = op.type
            registers_involved = [arg.reg_name for arg in args[:2]]
            if len(args) > 2 and args[2].reg_name not in registers_involved:
                # there is a distinct output register
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

        if optype == OpType.ClassicalExpBox:
            out_args = args[op.get_n_i() :]
            if (
                out_args
                == list(cregs[out_args[0].reg_name])[: op.get_n_io() + op.get_n_o()]
            ):
                stream_out.write(f"{out_args[0].reg_name} = {str(op.get_exp())};\n")
            elif len(out_args) == 1:
                stream_out.write(f"{out_args[0]} = {str(op.get_exp())};\n")
            else:
                raise QASMUnsupportedError(
                    f"ClassicalExpBox only supported"
                    " for writing to a single bit or whole registers."
                )
            continue
        if optype == OpType.WASM:
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
                    if bits != list(cregs[regname]):
                        QASMUnsupportedError("WASM ops must act on entire registers.")
                    reglist.append(regname)
            if outputs:
                stream_out.write(f"{', '.join(outputs)} = ")
            stream_out.write(f"{op.func_name}({', '.join(inputs)});\n")
            continue
        if optype == OpType.CustomGate:
            if op.gate.name not in include_gate_defs:
                # unroll custom gate
                # TODO when opaque gates are supported, make sure they are not
                # unrolled here
                def_circ = op.get_circuit()

                if def_circ.n_gates == 0:
                    raise QASMUnsupportedError(
                        f"CustomGate {op.gate.name} has empty definition."
                        " Empty CustomGates and opaque gates are not supported."
                    )
                def_circ.rename_units(dict(zip(def_circ.qubits, args)))
                def_circ.symbol_substitution(dict(zip(op.gate.args, op.params)))

                circuit_to_qasm_io(def_circ, stream_out, header, include_gate_defs)
                continue
            else:
                opstr = op.gate.name
        elif header == "hqslib1" and optype == OpType.ZZPhase:
            # special handling for zzphase
            opstr = "RZZ"
            params = op.params
        elif optype in _tk_to_qasm_noparams:
            opstr = _tk_to_qasm_noparams[optype]
        elif optype in _tk_to_qasm_params:
            opstr = _tk_to_qasm_params[optype]
        else:
            raise QASMUnsupportedError(
                "Cannot print command of type: {}".format(op.get_name())
            )
        if opstr not in include_gate_defs:
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
