// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <unordered_set>
#include <vector>

#include "OpType.hpp"

namespace tket {

/** Set of operation types */
typedef std::unordered_set<OpType> OpTypeSet;

/** Vector of operation types */
typedef std::vector<OpType> OpTypeVector;

/** Set of all elementary gates */
const OpTypeSet &all_gate_types();

/** Set of all gates over more than one qubit */
const OpTypeSet &all_multi_qubit_types();

/** Set of all gates over a single qubit */
const OpTypeSet &all_single_qubit_types();

/** Set of all single-qubit gates that can be expressed as TK1 */
const OpTypeSet &all_single_qubit_unitary_types();

/** Set of all measurement and reset gates */
const OpTypeSet &all_projective_types();

/** Test for initial, final and barrier "ops" */
bool is_metaop_type(OpType optype);

/** Test for input or creation quantum "ops" */
bool is_initial_q_type(OpType optype);

/** Test for output or discard quantum "ops" */
bool is_final_q_type(OpType optype);

/** Test for input, creation, output or discard quantum "ops" */
bool is_boundary_q_type(OpType optype);

/** Test for input or output for classical "ops" */
bool is_boundary_c_type(OpType optype);

/** Test for elementary gates */
bool is_gate_type(OpType optype);

/** Test for boxes (complex packaged operations) */
bool is_box_type(OpType optype);

/** Test for flowops (just for high-level control flow) */
bool is_flowop_type(OpType optype);

/** Test for rotations (including controlled rotations) */
bool is_rotation_type(OpType optype);

/** Test for rotations around Pauli axes */
bool is_parameterised_pauli_rotation_type(OpType optype);

/** Test for gates over more than one qubit */
bool is_multi_qubit_type(OpType optype);

/** Test for gates over a single qubit */
bool is_single_qubit_type(OpType optype);

/** Test for single-qubit gates that can be expressed as TK1 */
bool is_single_qubit_unitary_type(OpType optype);

/** Test for non-invertible operations */
bool is_oneway_type(OpType optype);

/** Test for Clifford operations */
bool is_clifford_type(OpType optype);

/** Test for measurement and reset gates */
bool is_projective_type(OpType optype);

/** Test for purely classical gates derived from ClassicalOp */
bool is_classical_type(OpType optype);

/** Whether a given operation type belongs to a given set */
bool find_in_set(const OpType &val, const OpTypeSet &set);

}  // namespace tket
