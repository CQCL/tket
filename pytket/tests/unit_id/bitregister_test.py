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
from pytket import Bit
from pytket.unit_id import BitRegister


def test_iteration() -> None:
    reg = BitRegister("test", 4)
    test_current = 0
    for bit in reg:
        assert isinstance(bit, Bit)
        assert bit.reg_name == "test"
        assert bit.index[0] == test_current
        test_current += 1  # noqa: SIM113
