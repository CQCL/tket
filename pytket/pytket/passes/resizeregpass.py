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

from collections.abc import Callable

from pytket.circuit import Bit, Circuit
from pytket.unit_id import _TEMP_BIT_NAME

from .._tket.passes import BasePass, CustomPass

MAX_C_REG_WIDTH = 32


def _is_scratch(bit: Bit) -> bool:
    reg_name = bit.reg_name
    return bool(reg_name == _TEMP_BIT_NAME) or reg_name.startswith(f"{_TEMP_BIT_NAME}_")


def _gen_scratch_transformation(max_size: int) -> Callable[[Circuit], Circuit]:
    def t(circuit: Circuit) -> Circuit:
        # Find all scratch bits
        scratch_bits = list(filter(_is_scratch, circuit.bits))
        # If the total number of scratch bits exceeds the max width, rename them
        if len(scratch_bits) > max_size:
            bits_map = {}
            for i, bit in enumerate(scratch_bits):
                bits_map[bit] = Bit(f"{_TEMP_BIT_NAME}_{i // max_size}", i % max_size)
            circuit.rename_units(bits_map)  # type: ignore
        return circuit

    return t


def scratch_reg_resize_pass(max_size: int = MAX_C_REG_WIDTH) -> BasePass:
    """Create a pass that breaks up the internal scratch bit registers into smaller
    registers.

    :param max_size: desired maximum size of scratch bit registers
    :return: a pass to break up the scratch registers
    """
    return CustomPass(
        _gen_scratch_transformation(max_size), label="resize scratch bits"
    )
