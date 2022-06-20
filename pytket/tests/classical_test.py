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

import operator
import re
from typing import Callable, Dict, List, Optional, Tuple, Type, TypeVar, Union, cast
import json
from pathlib import Path

from jsonschema import validate  # type: ignore
from hypothesis import given, settings, strategies
from hypothesis.strategies import SearchStrategy

from pytket import wasm

import pytest

from pytket._tket.circuit import (  # type: ignore
    _TEMP_BIT_NAME,
    _TEMP_BIT_REG_BASE,
    BitRegister,
    QubitRegister,
    Bit,
    Circuit,
    OpType,
    Qubit,
    Conditional,
    Op,
    SetBitsOp,
    MultiBitOp,
    RangePredicateOp,
    ClassicalExpBox,
)
from pytket.circuit.logic_exp import (
    BinaryOp,
    BitLogicExp,
    BitWiseOp,
    ConstPredicate,
    LogicExp,
    RegEq,
    RegGeq,
    RegGt,
    RegLeq,
    RegLogicExp,
    RegLt,
    RegNeq,
    RegPow,
    RegWiseOp,
    UnaryOp,
    reg_eq,
    reg_geq,
    reg_gt,
    reg_leq,
    reg_lt,
    reg_neq,
    if_bit,
    if_not_bit,
)

from pytket.passes import DecomposeClassicalExp  # type: ignore

from strategies import reg_name_regex, binary_digits, uint32  # type: ignore

curr_file_path = Path(__file__).resolve().parent


with open(curr_file_path.parent.parent / "schemas/circuit_v1.json", "r") as f:
    schema = json.load(f)


def print_commands(c: Circuit) -> None:
    print("\n".join(map(str, c)))


def test_c_ops() -> None:
    c = Circuit(0, 4)
    cl_cx_values = [0, 3, 2, 1]  # classical CX transform
    c.add_c_transform(cl_cx_values, [0, 1], "ClCX")
    c.add_c_transform(cl_cx_values, [Bit(1), Bit(2)], "ClCX")
    eq_pred_values = [True, False, False, True]  # test 2 bits for equality
    c.add_c_predicate(eq_pred_values, [0, 1], 2, "EQ")
    c.add_c_predicate(eq_pred_values, [Bit(1), Bit(2)], Bit(3), "EQ")
    and_values = [0, 0, 0, 1]  # binary AND
    c.add_c_modifier(and_values, [1], 2)
    c.add_c_modifier(and_values, [Bit(2)], Bit(3))
    c.add_c_and(1, 2, 3)
    c.add_c_or(Bit(2), Bit(3), Bit(0))
    c.add_c_not(0, 1)
    c.add_c_and(2, 3, 3)
    c.add_c_or(0, 3, 0)
    c.add_c_range_predicate(2, 6, [0, 1, 2], 3)
    c.add_c_setbits([True], [1])
    c.add_c_setbits([False, True], [Bit(2), Bit(3)])
    c.add_c_copybits([0, 1], [2, 3])
    c.add_c_range_predicate(2, 6, [1, 2, 3], 0)
    c0 = c.add_c_register("c0", 3)
    c1 = c.add_c_register("c1", 4)
    c2 = c.add_c_register("c2", 5)
    c.add_c_and_to_registers(c0, c1, c2)
    c.add_c_xor_to_registers(c2, c1, c2)
    c.add_c_not_to_registers(c1, c2)

    c.add_c_setreg(3, c0)
    c.add_c_copyreg(c1, c0)

    cmds = c.get_commands()
    assert len(cmds) == 21
    assert len([cmd for cmd in cmds if cmd.op.get_name() == "AND"]) == 2
    rp_cmds = [cmd for cmd in cmds if cmd.op.type == OpType.RangePredicate]
    assert len(rp_cmds) == 2
    assert rp_cmds[0].op == rp_cmds[1].op
    mb_cmds = [cmd for cmd in cmds if cmd.op.type == OpType.MultiBit]
    assert len(mb_cmds) == 3
    assert (
        mb_cmds[0] != mb_cmds[1]
        and mb_cmds[0] != mb_cmds[2]
        and mb_cmds[1] != mb_cmds[2]
    )
    assert str(cmds[6]) == "CopyBits c1[0], c1[1], c1[2], c0[0], c0[1], c0[2];"
    assert str(cmds[3]) == "SetBits(110) c0[0], c0[1], c0[2];"


def test_wasm() -> None:
    c = Circuit(0, 6)
    c._add_wasm("funcname", "wasmfileuid", [1, 1], [], [Bit(0), Bit(1)])
    c._add_wasm("funcname", "wasmfileuid", [1, 1], [], [Bit(0), Bit(2)])
    c._add_wasm("funcname", "wasmfileuid", [1, 1], [2], [0, 1, 2, 3])
    c._add_wasm("funcname", "wasmfileuid", [1, 1], [2], [0, 1, 2, 4])
    c._add_wasm("funcname", "wasmfileuid", [1], [1, 2], [0, 1, 2, 3])
    c._add_wasm("funcname", "wasmfileuid", [2, 1], [3], [0, 1, 2, 3, 4, 5])

    # the boxes with no output are not counted
    assert c.depth() == 4


def test_wasm_2() -> None:
    c = Circuit(6, 6)
    c0 = c.add_c_register("c0", 3)
    c1 = c.add_c_register("c1", 4)
    c2 = c.add_c_register("c2", 5)

    c._add_wasm("funcname", "wasmfileuid", [c0, c1], [c2])

    # the boxes with no output are not counted
    assert c.depth() == 1


def test_wasm_3() -> None:
    w = wasm.WasmFileHandler("testfile.wasm")

    c = Circuit(0, 6)

    c.add_wasm("funcname", w, [1], [1], [Bit(0), Bit(1)])

    assert c.depth() == 1


def test_wasm_handler() -> None:
    w = wasm.WasmFileHandler("testfile.wasm")

    with pytest.raises(ValueError):
        w2 = wasm.WasmFileHandler("testfile-2.wasm")


def test_add_wasm_to_reg() -> None:
    w = wasm.WasmFileHandler("testfile.wasm")

    c = Circuit(6, 6)
    c0 = c.add_c_register("c0", 3)
    c1 = c.add_c_register("c1", 4)
    c2 = c.add_c_register("c2", 5)

    c.add_wasm_to_reg("funcname", w, [c0, c1], [c2])
    c.add_wasm_to_reg("funcname2", w, [c2], [c2])

    assert c.depth() == 2


def gen_reg(
    name: str, size: int, reg_type: Union[Type[BitRegister], Type[QubitRegister]]
) -> Union[BitRegister, QubitRegister]:
    return reg_type(name, size)


@strategies.composite
def bit_register(
    draw: Callable,
    name: SearchStrategy[str] = strategies.from_regex(reg_name_regex, fullmatch=True),
    size: SearchStrategy[int] = strategies.integers(min_value=2, max_value=32),
) -> BitRegister:
    return cast(
        BitRegister,
        gen_reg(
            draw(name.filter(lambda nm: not nm.startswith("q"))),
            draw(size),
            BitRegister,
        ),
    )


@strategies.composite
def qubit_register(
    draw: Callable,
    name: SearchStrategy[str] = strategies.from_regex(reg_name_regex, fullmatch=True),
    size: SearchStrategy[int] = strategies.integers(min_value=0, max_value=32),
) -> QubitRegister:
    return cast(QubitRegister, gen_reg(draw(name), draw(size), QubitRegister))


@given(
    reg=strategies.one_of(bit_register(), qubit_register()),
    index=strategies.integers(min_value=0, max_value=32),
)
def test_registers(reg: Union[BitRegister, QubitRegister], index: int) -> None:

    unit_type = Qubit if type(reg) is QubitRegister else Bit
    if index < reg.size:
        assert reg[index] == unit_type(reg.name, index)

    assert [x for x in reg] == [unit_type(reg.name, i) for i in range(reg.size)]


@strategies.composite
def bits(
    draw: Callable,
    name: SearchStrategy[str] = strategies.from_regex(reg_name_regex, fullmatch=True),
    index: SearchStrategy[int] = uint32,
) -> Bit:
    return Bit(draw(name.filter(lambda nm: not nm.startswith("q"))), draw(index))


@strategies.composite
def primitive_bit_logic_exps(
    draw: Callable,
    ops: SearchStrategy[BitWiseOp] = strategies.sampled_from(BitWiseOp),
    bits: SearchStrategy[Bit] = bits(),
) -> BitLogicExp:
    op = draw(ops)

    exp_type = LogicExp.factory(op)
    args: List[Bit] = [draw(bits)]
    if issubclass(exp_type, BinaryOp):
        if issubclass(exp_type, ConstPredicate):
            const_compare = draw(binary_digits)
            args.append(const_compare)
        else:
            args.append(draw(bits))

    exp = exp_type(*args)
    assert isinstance(exp, BitLogicExp)
    return exp


def overflow_wrapper(f: Callable[..., int], maxval: int) -> Callable[..., int]:
    def wrapper(*args: int) -> int:
        return f(*args) % maxval

    return wrapper


@given(
    bit_exp=primitive_bit_logic_exps(),
    constants=strategies.tuples(binary_digits, binary_digits),
)
def test_bit_exp(bit_exp: BitLogicExp, constants: Tuple[int, int]) -> None:
    iter_c = iter(constants)
    for inp in bit_exp.all_inputs():
        bit_exp.set_value(inp, next(iter_c))
    op_map: Dict[BitWiseOp, Callable] = {
        BitWiseOp.AND: operator.and_,
        BitWiseOp.OR: operator.or_,
        BitWiseOp.XOR: operator.xor,
        BitWiseOp.NOT: operator.not_,
        BitWiseOp.EQ: operator.eq,
        BitWiseOp.NEQ: operator.ne,
    }
    op_map = {key: overflow_wrapper(val, 2) for key, val in op_map.items()}
    eval_val = bit_exp.eval_vals()

    assert eval_val in (0, 1)

    correct_val = op_map[cast(BitWiseOp, bit_exp.op)](*bit_exp.args)

    if isinstance(correct_val, bool):
        correct_val = int(correct_val)

    assert eval_val == correct_val


@strategies.composite
def primitive_reg_logic_exps(
    draw: Callable,
    ops: SearchStrategy[RegWiseOp] = strategies.sampled_from(RegWiseOp),
    bit_regs: SearchStrategy[BitRegister] = bit_register(),
) -> RegLogicExp:
    op = draw(ops)

    exp_type = LogicExp.factory(op)
    args: List[BitRegister] = [draw(bit_regs)]
    if issubclass(exp_type, BinaryOp):
        if issubclass(
            exp_type,
            (
                RegEq,
                RegNeq,
                RegLt,
                RegGt,
                RegLeq,
                RegGeq,
            ),
        ):
            const_compare = draw(uint32)
            args.append(const_compare)
        else:
            args.append(draw(bit_regs))
    else:
        assert issubclass(exp_type, UnaryOp)
    exp = exp_type(*args)  # type:ignore
    assert isinstance(exp, RegLogicExp)
    return exp


@given(
    reg_exp=primitive_reg_logic_exps(),
    constants=strategies.tuples(
        uint32,
        uint32,
    ),
)
def test_reg_exp(reg_exp: RegLogicExp, constants: Tuple[int, int]) -> None:
    if isinstance(reg_exp, RegPow):
        # to stop massive numbers
        constants = (min(1000, constants[0]), min(constants[1], 3))
    iter_c = iter(constants)
    for inp in reg_exp.all_inputs():
        reg_exp.set_value(inp, next(iter_c))

    op_map: Dict[RegWiseOp, Callable] = {
        RegWiseOp.AND: operator.and_,
        RegWiseOp.OR: operator.or_,
        RegWiseOp.XOR: operator.xor,
        RegWiseOp.NEQ: operator.ne,
        RegWiseOp.EQ: operator.eq,
        RegWiseOp.LT: operator.lt,
        RegWiseOp.GT: operator.gt,
        RegWiseOp.LEQ: operator.le,
        RegWiseOp.GEQ: operator.ge,
    }
    unsupported_ops = {
        RegWiseOp.ADD,
        RegWiseOp.SUB,
        RegWiseOp.MUL,
        RegWiseOp.DIV,
        RegWiseOp.POW,
        RegWiseOp.LSH,
        RegWiseOp.RSH,
        RegWiseOp.NOT,
        RegWiseOp.NEG,
    }
    eval_val = reg_exp.eval_vals()
    op = cast(RegWiseOp, reg_exp.op)

    correct_val = reg_exp if op in unsupported_ops else op_map[op](*reg_exp.args)
    if isinstance(correct_val, bool):
        correct_val = int(correct_val)

    assert eval_val == correct_val


@strategies.composite
def composite_bit_logic_exps(
    draw: Callable,
    bits: SearchStrategy[Bit] = bits(),
    constants: SearchStrategy[int] = binary_digits,
    operators: SearchStrategy[Callable] = strategies.sampled_from(
        [
            operator.and_,
            operator.or_,
            operator.xor,
        ]
    ),
    n_terms: SearchStrategy[int] = strategies.integers(min_value=2, max_value=7),
) -> BitLogicExp:
    terms = draw(n_terms)
    exp: BitLogicExp
    exp = draw(bits)
    for _ in range(terms):
        chosen_operator = draw(operators)
        if chosen_operator is operator.not_:
            exp = chosen_operator(exp)
        else:
            exp = chosen_operator(exp, draw(strategies.one_of(bits, constants)))
    return exp


@strategies.composite
def composite_reg_logic_exps(
    draw: Callable,
    regs: SearchStrategy[BitRegister] = bit_register(),
    constants: SearchStrategy[int] = uint32,
    operators: SearchStrategy[Callable] = strategies.sampled_from(
        [
            operator.and_,
            operator.or_,
            operator.xor,
        ]
    ),
    n_terms: SearchStrategy[int] = strategies.integers(min_value=1, max_value=7),
) -> RegLogicExp:
    terms = draw(n_terms)
    exp = draw(regs)
    used_reg_names = {exp.name}
    for _ in range(terms):
        chosen_operator = draw(operators)
        if chosen_operator is operator.not_:
            exp = chosen_operator(exp)
        else:
            second_operand = draw(
                strategies.one_of(
                    regs.filter(lambda x: x.name not in used_reg_names), constants
                )
            )
            exp = chosen_operator(exp, second_operand)
            if isinstance(second_operand, BitRegister):
                used_reg_names.add(second_operand.name)
    return cast(RegLogicExp, exp)


@strategies.composite
def bit_const_predicates(
    draw: Callable,
    operators: SearchStrategy[Callable] = strategies.sampled_from([if_bit, if_not_bit]),
    exp: SearchStrategy[BitLogicExp] = composite_bit_logic_exps(),
) -> ConstPredicate:
    return cast(ConstPredicate, draw(operators)(draw(exp)))


@strategies.composite
def reg_const_predicates(
    draw: Callable,
    exp: SearchStrategy[RegLogicExp] = composite_reg_logic_exps(),
    operators: SearchStrategy[Callable] = strategies.sampled_from(
        [reg_eq, reg_neq, reg_lt, reg_gt, reg_leq, reg_geq]
    ),
    constants: SearchStrategy[int] = uint32,
) -> ConstPredicate:
    return cast(ConstPredicate, draw(operators)(draw(exp), draw(constants)))


@given(condition=strategies.one_of(bit_const_predicates(), reg_const_predicates()))
@settings(print_blob=True, deadline=None)
def test_const_predicate(condition: ConstPredicate) -> None:
    assert isinstance(condition, ConstPredicate)
    # test serialization round trip here
    # as condition should be the most nested structure
    assert LogicExp.from_dict(condition.to_dict()) == condition
    circ = Circuit()
    qb = Qubit("q_", 1)
    # no bits should start with "q" due to filters
    circ.add_qubit(qb)
    if isinstance(condition, RegLogicExp):
        reg_added = set()
        for inp in condition.all_inputs():
            if inp not in reg_added:
                reg = cast(BitRegister, inp)
                circ.add_c_register(reg)
                reg_added.add(inp)
    else:
        for inp in condition.all_inputs():
            bit = cast(Bit, inp)
            circ.add_bit(bit, reject_dups=False)

    circ.X(qb, condition=condition)
    assert circ.n_gates_of_type(OpType.ClassicalExpBox) == 1
    newcirc = circ.copy()
    DecomposeClassicalExp().apply(newcirc)

    assert newcirc.n_gates_of_type(OpType.ClassicalExpBox) == 0

    assert newcirc.n_gates >= circ.n_gates
    if isinstance(condition, RegLogicExp):
        commands = newcirc.get_commands()
        check_range = commands[-2].op.type == OpType.RangePredicate
        if not check_range:
            print_commands(commands)
            print(condition.args[0])
        assert check_range


def compare_commands_box(
    circ1: Circuit, circ2: Circuit, print_com: bool = False
) -> bool:
    assert circ1.bits == circ2.bits
    assert circ1.qubits == circ2.qubits
    assert circ1.qubits == circ2.qubits
    assert circ1.name == circ2.name

    commands_equal = True
    for c1, c2 in zip(circ1, circ2):
        if print_com:
            print(c1, c2)
        if c1.op.type == OpType.ClassicalExpBox:
            commands_equal &= c1.op.content_equality(c2.op)
            commands_equal &= c1.args == c2.args
        elif c1.op.type == OpType.Conditional:
            commands_equal &= c1.op.value == c2.op.value
            commands_equal &= c1.op.width == c2.op.width
            if c1.op.op.type == OpType.ClassicalExpBox:
                commands_equal &= c1.op.op.content_equality(c2.op.op)
            else:
                commands_equal &= c1.op == c2.op
            commands_equal &= c1.args == c2.args
        else:
            commands_equal &= c1 == c2
        if not commands_equal:
            print(f"Commands not equal\n{c1}\n\n{c2}\n")
            break
    return commands_equal


def test_compare_commands_box() -> None:
    c = Circuit(2)
    bits = [Bit(i) for i in range(2)]
    for b in bits:
        c.add_bit(b)
    c.X(0).X(0)
    c.X(0, condition=(bits[1] & bits[0]))

    assert compare_commands_box(c, c.copy())
    c2 = Circuit(2)
    for b in bits:
        c2.add_bit(b)
    c2.X(0).X(0)
    # swap bits in expression
    c2.X(0, condition=(bits[0] & bits[1]))
    assert not compare_commands_box(c, c2)


def check_serialization_roundtrip(circ: Circuit) -> None:
    circ_dict = circ.to_dict()
    circ_from_dict = Circuit.from_dict(circ_dict)
    validate(instance=circ_dict, schema=schema)
    assert compare_commands_box(circ_from_dict, circ)
    assert circ_from_dict.to_dict() == circ_dict


def test_decomposition_known() -> None:
    bits = [Bit(i) for i in range(10)]
    registers = [BitRegister(c, 3) for c in "abdefghijk"]

    qreg = QubitRegister("q_", 10)
    circ = Circuit()
    conditioned_circ = Circuit()
    decomposed_circ = Circuit()

    for c in (circ, conditioned_circ, decomposed_circ):
        for b in bits:
            c.add_bit(b)
        for br in registers:
            for b in br:
                c.add_bit(b, reject_dups=False)
        c.add_q_register(qreg.name, qreg.size)

    circ.H(qreg[0], condition=bits[0])
    circ.X(qreg[0], condition=if_bit(bits[1]))
    circ.S(qreg[0])
    circ.T(qreg[1], condition=if_not_bit(bits[2]))
    circ.Z(qreg[0], condition=(bits[2] & bits[3]))
    circ.Z(qreg[1], condition=if_not_bit(bits[3] & bits[4]))
    big_exp = bits[4] | bits[5] ^ bits[6] | bits[7] & bits[8]
    # ^ no need for parantheses as python operator precedence
    # will enforce correct precedence in LogicExp
    circ.CX(qreg[0], qreg[1])
    circ.CX(qreg[1], qreg[2], condition=big_exp)

    circ.add_barrier(list(qreg))

    circ.H(qreg[2], condition=reg_eq(registers[0], 3))
    circ.X(qreg[3], condition=reg_lt(registers[1], 6))
    circ.Y(qreg[4], condition=reg_neq(registers[2], 5))
    circ.Z(qreg[5], condition=reg_gt(registers[3], 3))
    circ.S(qreg[6], condition=reg_leq(registers[4], 6))
    circ.T(qreg[7], condition=reg_geq(registers[5], 3))
    big_reg_exp = registers[4] & registers[3] | registers[6] ^ registers[7]
    circ.CX(qreg[3], qreg[4], condition=reg_eq(big_reg_exp, 3))

    circ.add_classicalexpbox_bit(
        bits[4] | bits[5] & bits[3], [bits[0]], condition=bits[1]
    )
    check_serialization_roundtrip(circ)

    temp_bits = BitRegister(_TEMP_BIT_NAME, 32)

    def temp_reg(i: int) -> BitRegister:
        return BitRegister(f"{_TEMP_BIT_REG_BASE}_{i}", 32)

    for b in (temp_bits[i] for i in range(0, 10)):
        conditioned_circ.add_bit(b)

    for t_r in (temp_reg(i) for i in range(0, 1)):
        conditioned_circ.add_c_register(t_r.name, t_r.size)

    # relies on existing interface for adding conditionals
    # may need a more low level interface for that if we decide to get rid of it
    conditioned_circ.H(qreg[0], condition_bits=[bits[0]], condition_value=1)
    conditioned_circ.X(qreg[0], condition_bits=[bits[1]], condition_value=1)
    conditioned_circ.S(qreg[0])
    conditioned_circ.T(qreg[1], condition_bits=[bits[2]], condition_value=0)

    conditioned_circ.add_classicalexpbox_bit((bits[2] & bits[3]), [temp_bits[0]])
    conditioned_circ.Z(qreg[0], condition_bits=[temp_bits[0]], condition_value=1)
    conditioned_circ.add_classicalexpbox_bit((bits[3] & bits[4]), [temp_bits[1]])
    conditioned_circ.Z(qreg[1], condition_bits=[temp_bits[1]], condition_value=0)
    conditioned_circ.CX(qreg[0], qreg[1])
    conditioned_circ.add_classicalexpbox_bit(big_exp, [temp_bits[2]])
    conditioned_circ.CX(
        qreg[1], qreg[2], condition_bits=[temp_bits[2]], condition_value=1
    )

    conditioned_circ.add_barrier(list(qreg))

    conditioned_circ.add_c_range_predicate(3, 3, list(registers[0]), temp_bits[3])
    conditioned_circ.H(qreg[2], condition_bits=[temp_bits[3]], condition_value=1)
    conditioned_circ.add_c_range_predicate(0, 5, list(registers[1]), temp_bits[4])
    conditioned_circ.X(qreg[3], condition_bits=[temp_bits[4]], condition_value=1)
    conditioned_circ.add_c_range_predicate(5, 5, list(registers[2]), temp_bits[5])
    conditioned_circ.Y(qreg[4], condition_bits=[temp_bits[5]], condition_value=0)
    conditioned_circ.add_c_range_predicate(
        4, 4294967295, list(registers[3]), temp_bits[6]
    )
    conditioned_circ.Z(qreg[5], condition_bits=[temp_bits[6]], condition_value=1)
    conditioned_circ.add_c_range_predicate(0, 6, list(registers[4]), temp_bits[7])
    conditioned_circ.S(qreg[6], condition_bits=[temp_bits[7]], condition_value=1)
    conditioned_circ.add_c_range_predicate(
        3, 4294967295, list(registers[5]), temp_bits[8]
    )
    conditioned_circ.T(qreg[7], condition_bits=[temp_bits[8]], condition_value=1)

    temp_reg_bits = list(temp_reg(0))[:3]
    conditioned_circ.add_classicalexpbox_register(big_reg_exp, temp_reg_bits)
    conditioned_circ.add_c_range_predicate(3, 3, temp_reg_bits, temp_bits[9])
    conditioned_circ.CX(
        qreg[3], qreg[4], condition_bits=[temp_bits[9]], condition_value=1
    )
    conditioned_circ.add_classicalexpbox_bit(
        bits[4] | bits[5] & bits[3], [bits[0]], condition=bits[1]
    )

    assert compare_commands_box(circ, conditioned_circ)

    for b in (temp_bits[i] for i in range(0, 11)):
        decomposed_circ.add_bit(b)

    decomposed_circ.add_c_register(BitRegister(f"{_TEMP_BIT_REG_BASE}_0", 3))
    decomposed_circ.add_c_register(BitRegister(f"{_TEMP_BIT_REG_BASE}_1", 32))

    decomposed_circ.H(qreg[0], condition_bits=[bits[0]], condition_value=1)
    decomposed_circ.X(qreg[0], condition_bits=[bits[1]], condition_value=1)
    decomposed_circ.S(qreg[0])
    decomposed_circ.T(qreg[1], condition_bits=[bits[2]], condition_value=0)
    decomposed_circ.add_c_and(bits[2], bits[3], temp_bits[0])
    decomposed_circ.Z(qreg[0], condition_bits=[temp_bits[0]], condition_value=1)
    decomposed_circ.add_c_and(bits[3], bits[4], temp_bits[1])
    decomposed_circ.Z(qreg[1], condition_bits=[temp_bits[1]], condition_value=0)
    decomposed_circ.CX(qreg[0], qreg[1])
    decomposed_circ.add_c_range_predicate(3, 3, list(registers[0]), temp_bits[3])
    decomposed_circ.add_c_range_predicate(0, 5, list(registers[1]), temp_bits[4])
    decomposed_circ.add_c_range_predicate(5, 5, list(registers[2]), temp_bits[5])
    decomposed_circ.add_c_range_predicate(
        4, 4294967295, list(registers[3]), temp_bits[6]
    )
    decomposed_circ.add_c_range_predicate(0, 6, list(registers[4]), temp_bits[7])
    decomposed_circ.add_c_range_predicate(
        3, 4294967295, list(registers[5]), temp_bits[8]
    )

    decomposed_circ.add_c_xor(bits[5], bits[6], temp_bits[2])
    decomposed_circ.add_c_and(bits[7], bits[8], temp_bits[10])
    decomposed_circ.add_c_or(bits[4], temp_bits[2], temp_bits[2])
    decomposed_circ.add_c_or(temp_bits[10], temp_bits[2], temp_bits[2])
    decomposed_circ.CX(
        qreg[1], qreg[2], condition_bits=[temp_bits[2]], condition_value=1
    )

    decomposed_circ.add_barrier(list(qreg))

    decomposed_circ.H(qreg[2], condition_bits=[temp_bits[3]], condition_value=1)
    decomposed_circ.X(qreg[3], condition_bits=[temp_bits[4]], condition_value=1)
    decomposed_circ.Y(qreg[4], condition_bits=[temp_bits[5]], condition_value=0)
    decomposed_circ.Z(qreg[5], condition_bits=[temp_bits[6]], condition_value=1)
    decomposed_circ.S(qreg[6], condition_bits=[temp_bits[7]], condition_value=1)
    decomposed_circ.T(qreg[7], condition_bits=[temp_bits[8]], condition_value=1)

    decomposed_circ.add_c_and_to_registers(registers[4], registers[3], temp_reg(0))
    decomposed_circ.add_c_xor_to_registers(registers[6], registers[7], temp_reg(1))
    decomposed_circ.add_c_or_to_registers(
        temp_reg(0), BitRegister(temp_reg(1).name, 3), temp_reg(0)
    )
    decomposed_circ.add_c_range_predicate(3, 3, list(temp_reg(0))[:3], temp_bits[9])
    decomposed_circ.CX(
        qreg[3], qreg[4], condition_bits=[temp_bits[9]], condition_value=1
    )
    decomposed_circ.add_c_and(
        bits[5], bits[3], bits[0], condition_bits=[bits[1]], condition_value=1
    )
    decomposed_circ.add_c_or(
        bits[4], bits[0], bits[0], condition_bits=[bits[1]], condition_value=1
    )
    check_serialization_roundtrip(decomposed_circ)
    circ_copy = circ.copy()

    DecomposeClassicalExp().apply(circ_copy)
    assert circ_copy == decomposed_circ


def test_conditional() -> None:
    c = Circuit(1, 2)
    c.H(0, condition_bits=[0, 1], condition_value=3)
    op = c.get_commands()[0].op
    cond_op = Conditional(Op.create(OpType.H), 2, 3)
    assert op == cond_op


def test_classical_ops() -> None:
    set_bits = SetBitsOp([True, True, False])
    multi_bit = MultiBitOp(set_bits, 2)
    range_predicate = RangePredicateOp(6, 27, 27)
    c = Circuit(0, 7)
    c.add_gate(multi_bit, [0, 1, 2, 3, 4, 5])
    c.add_gate(range_predicate, [0, 1, 2, 3, 4, 5, 6])
    exp = Bit(2) & Bit(3)
    c.add_classicalexpbox_bit(exp, [Bit(4)])
    cmds = c.get_commands()
    assert cmds[0].op.type == OpType.MultiBit
    assert cmds[1].op.type == OpType.RangePredicate
    ceb = ClassicalExpBox(2, 0, 1, exp)
    op2 = cmds[2].op
    assert ceb.get_exp() == op2.get_exp()


def test_add_expbox_bug() -> None:
    # previously a bug where if IO args weren't
    # at the back of the input iterator, the op signature
    # and wiring were incorrect
    c = Circuit()
    b = c.add_c_register("b", 2)
    c.add_classicalexpbox_bit(b[0] & b[1], [b[0]])
    com = c.get_commands()[0]

    assert com.op.get_n_i() == 1
    assert com.op.get_n_io() == 1
    assert com.op.get_n_o() == 0

    assert com.args == [b[1], b[0]]

    b1 = c.add_c_register("b1", 2)
    c.add_classicalexpbox_register(b | b1, list(b))
    com = c.get_commands()[1]

    assert com.op.get_n_i() == 2
    assert com.op.get_n_io() == 2
    assert com.op.get_n_o() == 0

    assert com.args == [b1[0], b1[1], b[0], b[1]]


def test_conditional_classicals() -> None:
    c = Circuit()
    b = c.add_c_register("b", 2)
    c.add_c_and(b[0], b[1], b[1], condition=b[0])
    assert str(c.get_commands()[0]) == "IF ([b[0]] == 1) THEN AND b[0], b[1];"


def test_conditional_wasm() -> None:
    c = Circuit(0, 6)
    b = c.add_c_register("b", 2)
    c._add_wasm("funcname", "wasmfileuid", [1, 1], [], [Bit(0), Bit(1)], condition=b[0])

    assert c.depth() == 1
    assert str(c.get_commands()[0]) == "IF ([b[0]] == 1) THEN WASM c[0], c[1];"


def test_conditional_wasm_ii() -> None:
    c = Circuit(0, 6)
    b = c.add_c_register("b", 2)
    c._add_wasm("funcname", "wasmfileuid", [b], [], condition=b[0])

    assert c.depth() == 1
    assert str(c.get_commands()[0]) == "IF ([b[0]] == 1) THEN WASM b[0], b[1];"


def test_conditional_wasm_iii() -> None:
    w = wasm.WasmFileHandler("testfile.wasm")

    c = Circuit(6, 6)
    c0 = c.add_c_register("c0", 3)
    c1 = c.add_c_register("c1", 4)
    c2 = c.add_c_register("c2", 5)

    b = c.add_c_register("b", 2)

    c.add_wasm_to_reg("funcname", w, [c0, c1], [c2], condition=b[0])
    c.add_wasm_to_reg("funcname2", w, [c2], [c2], condition=b[1])

    assert c.depth() == 2
    assert (
        str(c.get_commands()[0])
        == "IF ([b[0]] == 1) THEN WASM c0[0], c0[1], c0[2], c1[0], c1[1], c1[2], c1[3], c2[0], c2[1], c2[2], c2[3], c2[4];"
    )
    assert (
        str(c.get_commands()[1])
        == "IF ([b[1]] == 1) THEN WASM c2[0], c2[1], c2[2], c2[3], c2[4], c2[0], c2[1], c2[2], c2[3], c2[4];"
    )


def test_conditional_wasm_iv() -> None:
    w = wasm.WasmFileHandler("testfile.wasm")
    c = Circuit(0, 6)
    b = c.add_c_register("controlreg", 2)
    c.add_wasm("funcname", w, [1], [1], [Bit(0), Bit(1)], condition=b[0])

    assert c.depth() == 1
    assert str(c.get_commands()[0]) == "IF ([controlreg[0]] == 1) THEN WASM c[0], c[1];"


def test_arithmetic_ops() -> None:
    circ = Circuit()
    a = circ.add_c_register("a", 3)
    b = circ.add_c_register("b", 3)
    c = circ.add_c_register("c", 3)

    circ.add_classicalexpbox_register(a + b // c, a)
    circ.add_classicalexpbox_register(b << 2, c)
    circ.add_classicalexpbox_register(c >> 2, b)
    circ.add_classicalexpbox_register(a**c - b, a)

    commands = circ.get_commands()
    assert all(com.op.type == OpType.ClassicalExpBox for com in commands)

    assert commands[0].args == list(b) + list(c) + list(a)
    assert commands[1].args == list(b) + list(c)
    assert commands[2].args == list(c) + list(b)
    assert commands[3].args == list(b) + list(c) + list(a)

    assert str(commands[0].op.get_exp()) == "(a + (b / c))"
    assert str(commands[1].op.get_exp()) == "(b << 2)"
    assert str(commands[2].op.get_exp()) == "(c >> 2)"
    assert str(commands[3].op.get_exp()) == "((a ** c) - b)"
