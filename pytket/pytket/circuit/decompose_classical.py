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
from heapq import heappop, heappush
from typing import (
    Callable,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Type,
    TypeVar,
    Union,
    Generic,
)

from pytket._tket.unit_id import (
    _TEMP_BIT_NAME,
    _TEMP_BIT_REG_BASE,
    _TEMP_REG_SIZE,
    BitRegister,
    Bit,
)
from pytket._tket.circuit import Circuit, ClassicalExpBox, Conditional, OpType
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
        self._heap: List[T] = []
        self._heap_vars: Set[T] = set()

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


def temp_reg_in_args(args: List[Bit]) -> Optional[BitRegister]:
    """If there are bits from a temporary register in the args, return it."""
    temp_reg_bits = [b for b in args if b.reg_name.startswith(_TEMP_BIT_REG_BASE)]
    if temp_reg_bits:
        temp_reg = BitRegister(temp_reg_bits[0].reg_name, _TEMP_REG_SIZE)
        return temp_reg
    return None


VarType = TypeVar("VarType", Type[Bit], Type[BitRegister])


def int_to_bools(val: Constant, width: int) -> List[bool]:
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
    [Union[RegLogicExp, BitLogicExp], Optional[Dict]],
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
        op_type: Union[Type[BitWiseOp], Type[RegWiseOp]] = BitWiseOp
        exp_type: Union[Type[BitLogicExp], Type[RegLogicExp]] = BitLogicExp
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
    def set_bits(var: Variable, val: Constant, kwargs: Dict) -> None:
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
        exp: Union[RegLogicExp, BitLogicExp], kwargs: Optional[Dict] = None
    ) -> Variable:
        assert isinstance(exp.op, op_type)
        kwargs = kwargs or {}
        # decompose children
        for i, sub_e in filter_by_type(exp.args, exp_type):
            assert isinstance(sub_e, (BitLogicExp, RegLogicExp))
            exp.args[i] = recursive_walk(sub_e, kwargs)
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
        return targ_bit  # type: ignore

    return recursive_walk


def _decompose_expressions(circ: Circuit) -> Tuple[Circuit, bool]:
    """Rewrite a circuit command-wise, decomposing ClassicalExpBox."""
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
    replace_targets: Dict[Variable, Variable] = dict()
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

                bit_heap.push(target)
                comp_bit = bit_recursive_walk(pred_exp, kwargs)

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

                reg_heap.push(out_reg)
                comp_reg = reg_recursive_walk(pred_exp, kwargs)
                if comp_reg.name != out_reg.name:  # type: ignore
                    replace_targets[out_reg] = comp_reg
            modified = True
            continue
        if optype == OpType.Barrier:
            # add_gate doesn't work for metaops
            newcirc.add_barrier(args)
        else:
            newcirc.add_gate(op, args, **kwargs)
    return newcirc, modified
