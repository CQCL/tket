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

"""Functions for decomposing Circuits containing classical expressions
 in to primitive logical operations."""
import copy
from collections.abc import Callable
from heapq import heappop, heappush
from typing import Any, Generic, TypeVar

from pytket._tket.circuit import (
    Circuit,
    ClassicalExpBox,
    ClBitVar,
    ClExpr,
    ClExprOp,
    ClOp,
    ClRegVar,
    Conditional,
    OpType,
    WiredClExpr,
)
from pytket._tket.unit_id import (
    _TEMP_BIT_NAME,
    _TEMP_BIT_REG_BASE,
    _TEMP_REG_SIZE,
    Bit,
    BitRegister,
)
from pytket.circuit.clexpr import check_register_alignments, has_reg_output
from pytket.circuit.logic_exp import (
    BitLogicExp,
    BitWiseOp,
    Constant,
    RegLogicExp,
    RegWiseOp,
    Variable,
    filter_by_type,
)

T = TypeVar("T")


class DecomposeClassicalError(Exception):
    """Error with decomposing classical operations."""


class VarHeap(Generic[T]):
    """A generic heap implementation."""

    def __init__(self) -> None:
        self._heap: list[T] = []
        self._heap_vars: set[T] = set()

    def pop(self) -> T:
        """Pop from top of heap."""
        return heappop(self._heap)

    def push(self, var: T) -> None:
        """Push var to heap."""
        heappush(self._heap, var)
        self._heap_vars.add(var)

    def is_heap_var(self, var: T) -> bool:
        """Check if var was generated from heap."""
        return var in self._heap_vars

    def fresh_var(self) -> T:
        """Generate new variable."""
        raise NotImplementedError


class BitHeap(VarHeap[Bit]):
    """Heap of temporary Bits."""

    def __init__(self, _reg_name: str = _TEMP_BIT_NAME):
        """Initialise new BitHeap.

        :param _reg_name: Name for register of Bits, defaults to _TEMP_BIT_NAME
        :type _reg_name: str, optional
        """

        self.reg_name = _reg_name
        super().__init__()

    @property
    def next_index(self) -> int:
        """Next available bit index, not used by any other heap bit."""
        return max((b.index[0] for b in self._heap_vars), default=-1) + 1

    def fresh_var(self) -> Bit:
        """Return Bit, from heap if available, otherwise create new."""
        if self._heap:
            return self.pop()
        new_bit = Bit(self.reg_name, self.next_index)
        self._heap_vars.add(new_bit)
        return new_bit


class RegHeap(VarHeap[BitRegister]):
    """Heap of temporary BitRegisters."""

    def __init__(self, _reg_name_base: str = _TEMP_BIT_REG_BASE):
        """Initialise new RegHeap.

        :param _reg_name_base: base string for register names, defaults to
            _TEMP_BIT_REG_BASE
        :type _reg_name_base: str, optional
        """
        self._reg_name_base = _reg_name_base
        super().__init__()

    @property
    def next_index(self) -> int:
        """Next available bit index, not used by any other heap register."""
        return (
            max((int(b.name.split("_")[-1]) for b in self._heap_vars), default=-1) + 1
        )

    def fresh_var(self, size: int = _TEMP_REG_SIZE) -> BitRegister:
        """Return BitRegister, from heap if available, otherwise create new.
        Optionally set size of created register."""
        if self._heap:
            return self.pop()
        new_reg = BitRegister(f"{self._reg_name_base}_{self.next_index}", size)
        self._heap_vars.add(new_reg)

        return new_reg


def temp_reg_in_args(args: list[Bit]) -> BitRegister | None:
    """If there are bits from a temporary register in the args, return it."""
    temp_reg_bits = [b for b in args if b.reg_name.startswith(_TEMP_BIT_REG_BASE)]
    if temp_reg_bits:
        return BitRegister(temp_reg_bits[0].reg_name, _TEMP_REG_SIZE)
    return None


VarType = TypeVar("VarType", type[Bit], type[BitRegister])


def int_to_bools(val: Constant, width: int) -> list[bool]:
    # map int to bools via litle endian encoding
    return list(map(bool, map(int, reversed(f"{val:0{width}b}"[-width:]))))


def get_bit_width(x: int) -> int:
    assert x >= 0
    c = 0
    while x:
        x >>= 1
        c += 1
    return c


def _gen_walk(var_type: VarType, newcirc: Circuit, heap: VarHeap) -> Callable[
    [RegLogicExp | BitLogicExp, Variable | None, dict | None],
    Variable,
]:
    """Generate a recursive walk method for decomposing an expression tree."""
    # map operation enum to circuit method
    _method_map = {
        BitWiseOp.AND: newcirc.add_c_and,
        BitWiseOp.OR: newcirc.add_c_or,
        BitWiseOp.XOR: newcirc.add_c_xor,
        RegWiseOp.AND: newcirc.add_c_and_to_registers,
        RegWiseOp.OR: newcirc.add_c_or_to_registers,
        RegWiseOp.XOR: newcirc.add_c_xor_to_registers,
    }
    if var_type is Bit:
        op_type: type[BitWiseOp] | type[RegWiseOp] = BitWiseOp
        exp_type: type[BitLogicExp] | type[RegLogicExp] = BitLogicExp
    else:
        assert var_type is BitRegister
        op_type = RegWiseOp
        exp_type = RegLogicExp

    def add_method(var: Variable) -> None:
        if isinstance(var, Bit):
            newcirc.add_bit(var, reject_dups=False)
        else:
            assert isinstance(var, BitRegister)
            for i in range(var.size):
                newcirc.add_bit(var.__getitem__(i), reject_dups=False)

    # method for setting bits during walk
    def set_bits(var: Variable, val: Constant, kwargs: dict) -> None:
        if isinstance(var, Bit):
            newcirc.add_c_setbits([bool(val)], [var], **kwargs)
        else:
            assert isinstance(var, BitRegister)
            bit_width = get_bit_width(val) if val else 1
            # make sure register size matches constant and add bits to circuit
            assert bit_width <= _TEMP_REG_SIZE
            reg = copy.copy(var)
            reg.size = bit_width
            add_method(reg)
            newcirc.add_c_setreg(val, reg, **kwargs)

    # convert an expression to gates on the circuit
    # and return the variable holding the result
    def recursive_walk(
        exp: RegLogicExp | BitLogicExp,
        targ_bit: Variable | None = None,
        kwargs: dict | None = None,
    ) -> Variable:
        assert isinstance(exp.op, op_type)
        kwargs = kwargs or {}
        # decompose children
        for i, sub_e in filter_by_type(exp.args, exp_type):
            assert isinstance(sub_e, (BitLogicExp, RegLogicExp))
            exp.args[i] = recursive_walk(sub_e, None, kwargs)
        # all args should now be Constant or Variable
        # write Constant to temporary Variable
        for idx, constant in filter_by_type(exp.args, Constant):
            fresh_var = heap.fresh_var()
            add_method(fresh_var)
            set_bits(fresh_var, constant, kwargs)

            exp.args[idx] = fresh_var

        arg1, arg2 = exp.args
        assert isinstance(arg1, var_type) and isinstance(arg2, var_type)
        # if argument is temporary, enable its reuse
        for arg in (arg1, arg2):
            if heap.is_heap_var(arg):
                heap.push(arg)

        if targ_bit is None:
            targ_bit = heap.fresh_var()
        add_method(targ_bit)
        # exp should now be a binary operation on args: List[Bit]

        try:
            _method_map[exp.op](arg1, arg2, targ_bit, **kwargs)  # type: ignore
        except KeyError as e:
            raise DecomposeClassicalError(
                f"{exp.op} cannot be decomposed to TKET primitives."
                " If targetting extended QASM you may not need to decompose."
            ) from e
        return targ_bit

    return recursive_walk


class ClExprDecomposer:
    def __init__(
        self,
        circ: Circuit,
        bit_posn: dict[int, int],
        reg_posn: dict[int, list[int]],
        args: list[Bit],
        bit_heap: BitHeap,
        reg_heap: RegHeap,
        kwargs: dict[str, Any],
    ):
        self.circ: Circuit = circ
        self.bit_posn: dict[int, int] = bit_posn
        self.reg_posn: dict[int, list[int]] = reg_posn
        self.args: list[Bit] = args
        self.bit_heap: BitHeap = bit_heap
        self.reg_heap: RegHeap = reg_heap
        self.kwargs: dict[str, Any] = kwargs
        # Construct maps from int (i.e. ClBitVar) to Bit, and from int (i.e. ClRegVar)
        # to BitRegister:
        self.bit_vars = {i: args[p] for i, p in bit_posn.items()}
        self.reg_vars = {
            i: BitRegister(args[p[0]].reg_name, len(p)) for i, p in reg_posn.items()
        }

    def add_var(self, var: Variable) -> None:
        """Add a Bit or BitRegister to the circuit if not already present."""
        if isinstance(var, Bit):
            self.circ.add_bit(var, reject_dups=False)
        else:
            assert isinstance(var, BitRegister)
            for bit in var.to_list():
                self.circ.add_bit(bit, reject_dups=False)

    def set_bits(self, var: Variable, val: int) -> None:
        """Set the value of a Bit or BitRegister."""
        assert val >= 0
        if isinstance(var, Bit):
            assert val >> 1 == 0
            self.circ.add_c_setbits([bool(val)], [var], **self.kwargs)
        else:
            assert isinstance(var, BitRegister)
            assert val >> var.size == 0
            self.circ.add_c_setreg(val, var, **self.kwargs)

    def decompose_expr(self, expr: ClExpr, out_var: Variable | None) -> Variable:
        """Add the decomposed expression to the circuit and return the Bit or
        BitRegister that contains the result.

        :param expr: the expression to decompose
        :param out_var: where to put the output (if None, create a new scratch location)
        """
        op: ClOp = expr.op
        heap: VarHeap = self.reg_heap if has_reg_output(op) else self.bit_heap

        # Eliminate (recursively) subsidiary expressions from the arguments, and convert
        # all terms to Bit or BitRegister:
        terms: list[Variable] = []
        for arg in expr.args:
            if isinstance(arg, int):
                # Assign to a fresh variable
                fresh_var = heap.fresh_var()
                self.add_var(fresh_var)
                self.set_bits(fresh_var, arg)
                terms.append(fresh_var)
            elif isinstance(arg, ClBitVar):
                terms.append(self.bit_vars[arg.index])
            elif isinstance(arg, ClRegVar):
                terms.append(self.reg_vars[arg.index])
            else:
                assert isinstance(arg, ClExpr)
                terms.append(self.decompose_expr(arg, None))

        # Enable reuse of temporary terms:
        for term in terms:
            if heap.is_heap_var(term):
                heap.push(term)

        if out_var is None:
            out_var = heap.fresh_var()
        self.add_var(out_var)
        match op:
            case ClOp.BitAnd:
                self.circ.add_c_and(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.BitNot:
                self.circ.add_c_not(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.BitOne:
                assert isinstance(out_var, Bit)
                self.circ.add_c_setbits([True], [out_var], **self.kwargs)
            case ClOp.BitOr:
                self.circ.add_c_or(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.BitXor:
                self.circ.add_c_xor(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.BitZero:
                assert isinstance(out_var, Bit)
                self.circ.add_c_setbits([False], [out_var], **self.kwargs)
            case ClOp.RegAnd:
                self.circ.add_c_and_to_registers(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.RegNot:
                self.circ.add_c_not_to_registers(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.RegOne:
                assert isinstance(out_var, BitRegister)
                self.circ.add_c_setbits(
                    [True] * out_var.size, out_var.to_list(), **self.kwargs
                )
            case ClOp.RegOr:
                self.circ.add_c_or_to_registers(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.RegXor:
                self.circ.add_c_xor_to_registers(*terms, out_var, **self.kwargs)  # type: ignore
            case ClOp.RegZero:
                assert isinstance(out_var, BitRegister)
                self.circ.add_c_setbits(
                    [False] * out_var.size, out_var.to_list(), **self.kwargs
                )
            case _:
                raise DecomposeClassicalError(
                    f"{op} cannot be decomposed to TKET primitives."
                )
        return out_var


def _decompose_expressions(circ: Circuit) -> tuple[Circuit, bool]:
    """Rewrite a circuit command-wise, decomposing ClassicalExpBox and ClExprOp."""
    if not check_register_alignments(circ):
        raise DecomposeClassicalError("Circuit contains non-register-aligned ClExprOp.")
    bit_heap = BitHeap()
    reg_heap = RegHeap()
    # add already used heap variables to heaps
    for b in circ.bits:
        if b.reg_name == _TEMP_BIT_NAME:
            bit_heap._heap_vars.add(b)
        elif b.reg_name.startswith(_TEMP_BIT_REG_BASE):
            reg_heap._heap_vars.add(BitRegister(b.reg_name, _TEMP_REG_SIZE))

    newcirc = Circuit(0, name=circ.name)

    for qb in circ.qubits:
        newcirc.add_qubit(qb)
    for cb in circ.bits:
        # lose all temporary bits, add back as required later
        if not (
            cb.reg_name.startswith(_TEMP_BIT_NAME)
            or cb.reg_name.startswith(_TEMP_BIT_REG_BASE)
        ):
            newcirc.add_bit(cb)

    # recursive functions for converting expressions to gates
    bit_recursive_walk = _gen_walk(Bit, newcirc, bit_heap)
    reg_recursive_walk = _gen_walk(BitRegister, newcirc, reg_heap)

    # targets of predicates that need to be relabelled
    replace_targets: dict[Variable, Variable] = dict()
    modified = False
    for command in circ:
        op = command.op
        optype = op.type
        args = command.args
        kwargs = dict()
        if optype == OpType.Conditional:
            assert isinstance(op, Conditional)
            bits = args[: op.width]
            # check if conditional on previously decomposed expression
            if len(bits) == 1 and bits[0] in replace_targets:
                assert isinstance(bits[0], Bit)
                # this op should encode comparison and value
                assert op.value in (0, 1)
                replace_bit = replace_targets[bits[0]]
                # temporary condition bit is available for reuse
                bit_heap.push(replace_bit)  # type: ignore

                # write new conditional op
                kwargs = {"condition_bits": [replace_bit], "condition_value": op.value}
            else:
                kwargs = {"condition_bits": bits, "condition_value": op.value}
            args = args[op.width :]
            op = op.op
            optype = op.type

        if optype == OpType.RangePredicate:
            target = args[-1]
            assert isinstance(target, Bit)
            newcirc.add_bit(target, reject_dups=False)
            temp_reg = temp_reg_in_args(args)  # type: ignore
            # ensure predicate is reading from correct output register
            if temp_reg in replace_targets:
                assert temp_reg is not None
                new_target = replace_targets[temp_reg]
                for i, a in enumerate(args):
                    if a.reg_name == temp_reg.name:
                        args[i] = Bit(new_target.name, a.index[0])  # type: ignore
                # operations conditional on this bit should remain so
                replace_targets[target] = target

        elif optype == OpType.ClassicalExpBox:
            assert isinstance(op, ClassicalExpBox)
            pred_exp = copy.deepcopy(op.get_exp())
            n_out_bits = op.get_n_o() + op.get_n_io()
            # copied as it will be modified in place
            if isinstance(pred_exp, BitLogicExp):
                assert n_out_bits == 1
                target = args[-1]
                assert isinstance(target, Bit)

                comp_bit = bit_recursive_walk(pred_exp, target, kwargs)

                replace_targets[target] = comp_bit

            else:
                assert isinstance(pred_exp, RegLogicExp)
                output_args = args[-n_out_bits:]
                if not all(
                    arg.reg_name == output_args[0].reg_name for arg in output_args
                ):
                    raise DecomposeClassicalError(
                        "Classical Expression must"
                        " only write to one Bit or one Register."
                    )
                out_reg = BitRegister(output_args[0].reg_name, len(output_args))

                comp_reg = reg_recursive_walk(pred_exp, out_reg, kwargs)
                if comp_reg.name != out_reg.name:  # type: ignore
                    replace_targets[out_reg] = comp_reg
            modified = True
            continue

        elif optype == OpType.ClExpr:
            assert isinstance(op, ClExprOp)
            wexpr: WiredClExpr = op.expr
            expr: ClExpr = wexpr.expr
            bit_posn = wexpr.bit_posn
            reg_posn = wexpr.reg_posn
            output_posn = wexpr.output_posn
            assert len(output_posn) > 0
            output0 = args[output_posn[0]]
            assert isinstance(output0, Bit)
            out_var: Variable = (
                BitRegister(output0.reg_name, len(output_posn))
                if has_reg_output(expr.op)
                else output0
            )
            decomposer = ClExprDecomposer(
                newcirc, bit_posn, reg_posn, args, bit_heap, reg_heap, kwargs  # type: ignore
            )
            comp_var = decomposer.decompose_expr(expr, out_var)
            if comp_var != out_var:
                replace_targets[out_var] = comp_var
            modified = True
            continue

        if optype == OpType.Barrier:
            # add_gate doesn't work for metaops
            newcirc.add_barrier(args)
        else:
            for arg in args:
                if (
                    isinstance(arg, Bit)
                    and arg.reg_name != "_w"  # workaround: this shouldn't be type Bit
                    and arg not in newcirc.bits
                ):
                    newcirc.add_bit(arg)
            newcirc.add_gate(op, args, **kwargs)
    return newcirc, modified
