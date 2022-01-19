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

#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Utils/Expression.hpp"

namespace tket {

/**
 * Get an operation with a given type, single parameter and qubit count
 *
 * @param chosen_type operation type
 * @param param operation parameter
 * @param n_qubits number of qubits (only necessary for gates and metaops
 *                 with variable quantum arity)
 */
Op_ptr get_op_ptr(OpType chosen_type, const Expr &param, unsigned n_qubits = 0);

/**
 * Get an operation from a type, vector of parameters and qubit count
 *
 * @param chosen_type operation type
 * @param params operation parameters
 * @param n_qubits number of qubits (only necessary for gates and metaops
 *                 with variable quantum arity)
 */
Op_ptr get_op_ptr(
    OpType chosen_type, const std::vector<Expr> &params = {},
    unsigned n_qubits = 0);

}  // namespace tket
