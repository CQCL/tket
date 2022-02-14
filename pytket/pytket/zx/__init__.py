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

"""
Python API module for interfacing with tket c++ implemented functions and modules.
Exports class ZXDiagram and associated tools.
"""

from pytket._tket.zx import *  # type: ignore
from .tensor_eval import (
    tensor_from_quantum_diagram,
    tensor_from_mixed_diagram,
    unitary_from_classical_diagram,
    unitary_from_quantum_diagram,
    density_matrix_from_cptp_diagram,
    fix_boundaries_to_binary_states,
    fix_inputs_to_binary_state,
    fix_outputs_to_binary_state,
)
