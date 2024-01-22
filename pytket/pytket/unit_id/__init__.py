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

from pytket._tket.unit_id import *
from pytket._tket.unit_id import Bit, BitRegister, Qubit, QubitRegister
from pytket._tket.unit_id import _TEMP_BIT_NAME
from pytket._tket.unit_id import _TEMP_BIT_REG_BASE


def _bitregister_next(self: BitRegister) -> Bit:
    if self._current < self.size:
        result = self[self._current]
        self._current += 1
        return result
    else:
        raise StopIteration


def _qubitregister_next(self: QubitRegister) -> Qubit:
    if self._current < self.size:
        result = self[self._current]
        self._current += 1
        return result
    else:
        raise StopIteration


setattr(BitRegister, "__next__", _bitregister_next)
setattr(QubitRegister, "__next__", _qubitregister_next)
