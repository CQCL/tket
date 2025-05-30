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

"""Enable adding of gates with conditions on Bit or BitRegister expressions."""

from pytket._tket.unit_id import _TEMP_BIT_NAME, _TEMP_BIT_REG_BASE
from pytket.circuit import Bit, BitRegister, Circuit
from pytket.circuit.clexpr import wired_clexpr_from_logic_exp
from pytket.circuit.logic_exp import (
    BitLogicExp,
    Constant,
    PredicateExp,
    RegEq,
    RegGeq,
    RegGt,
    RegLeq,
    RegLogicExp,
    RegLt,
    RegNeq,
)


class NonConstError(Exception):
    """A custom exception class for non constant predicate argument."""


def _add_condition(  # noqa: PLR0912
    circ: Circuit, condition: PredicateExp | Bit | BitLogicExp
) -> tuple[Bit, bool]:
    """Add a condition expression to a circuit using classical expression boxes,
    rangepredicates and conditionals. Return predicate bit and value of said bit.
    """
    if isinstance(condition, Bit):
        return condition, True
    if isinstance(condition, PredicateExp):
        pred_exp, pred_val = condition.args
        # PredicateExp constructor should ensure arg order
        if not isinstance(pred_val, Constant):
            raise NonConstError(
                "Condition expressions must be of type `PredicateExp`\
                with a constant second operand."
            )
    elif isinstance(condition, BitLogicExp):
        pred_val = 1
        pred_exp = condition
    else:
        raise ValueError(
            f"Condition {condition} must be of type Bit, BitLogicExp or PredicateExp"
        )

    next_index = (
        max(
            (bit.index[0] for bit in circ.bits if bit.reg_name == _TEMP_BIT_NAME),
            default=-1,
        )
        + 1
    )
    if isinstance(pred_exp, Bit):
        return pred_exp, bool(pred_val)

    # the resulting condition (a boolean) will be written to this
    # scratch bit
    condition_bit = Bit(_TEMP_BIT_NAME, next_index)
    circ.add_bit(condition_bit)

    if isinstance(pred_exp, BitLogicExp):
        wexpr, args = wired_clexpr_from_logic_exp(pred_exp, [condition_bit])
        circ.add_clexpr(wexpr, args)
        return condition_bit, bool(pred_val)

    assert isinstance(pred_exp, RegLogicExp | BitRegister)
    if isinstance(pred_exp, RegLogicExp):
        inps = pred_exp.all_inputs_ordered()
        reg_sizes: list[int] = []
        for reg in inps:
            assert isinstance(reg, BitRegister)
            reg_sizes.append(reg.size)
        min_reg_size = min(reg_sizes)
        existing_reg_names = {
            bit.reg_name
            for bit in circ.bits
            if bit.reg_name.startswith(_TEMP_BIT_REG_BASE)
        }
        existing_reg_indices = (
            int(r_name.split("_")[-1]) for r_name in existing_reg_names
        )
        next_index = max(existing_reg_indices, default=-1) + 1
        temp_reg = BitRegister(f"{_TEMP_BIT_REG_BASE}_{next_index}", min_reg_size)
        circ.add_c_register(temp_reg)
        target_bits = temp_reg.to_list()
        wexpr, args = wired_clexpr_from_logic_exp(pred_exp, target_bits)
        circ.add_clexpr(wexpr, args)
    elif isinstance(pred_exp, BitRegister):
        target_bits = pred_exp.to_list()

    minval = 0
    maxval = (1 << 64) - 1
    if isinstance(condition, RegLt):
        maxval = pred_val - 1
    elif isinstance(condition, RegGt):
        minval = pred_val + 1
    if isinstance(condition, RegLeq | RegEq | RegNeq):
        maxval = pred_val
    if isinstance(condition, RegGeq | RegEq | RegNeq):
        minval = pred_val

    circ.add_c_range_predicate(minval, maxval, target_bits, condition_bit)
    condition_value = not isinstance(condition, RegNeq)
    return condition_bit, condition_value
