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

"""Utility functions for performing high-level procedures in pytket"""

from .expectations import (
    expectation_from_shots,
    expectation_from_counts,
    get_pauli_expectation_value,
    get_operator_expectation_value,
)
from .measurements import append_pauli_measurement
from .prepare import prepare_circuit
from .results import (
    counts_from_shot_table,
    probs_from_counts,
    probs_from_state,
    permute_qubits_in_statevector,
    permute_basis_indexing,
    permute_rows_cols_in_unitary,
    compare_statevectors,
    compare_unitaries,
)
from .term_sequence import gen_term_sequence_circuit
from .operators import QubitPauliOperator
from .outcomearray import OutcomeArray, readout_counts
from .graph import Graph
from .symbolic import circuit_to_symbolic_unitary, circuit_apply_symbolic_statevector
from .distribution import (
    ProbabilityDistribution,
    EmpiricalDistribution,
    convex_combination,
)
from .stats import gate_counts
