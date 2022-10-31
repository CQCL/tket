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

#include "Circuit/Circuit.hpp"

namespace tket {

/**
 * Replace a multi-qubit operation with an equivalent circuit using TK2 gates
 *
 * @param op operation
 *
 * @return equivalent circuit
 *
 * @pre \p op is a multi-qubit operation
 *
 * @post The only multi-qubit gates in the replacement circuit are TK2
 */
Circuit TK2_circ_from_multiq(const Op_ptr op);

/**
 * Replace a multi-qubit operation with an equivalent circuit using CX gates
 *
 * @param op operation
 *
 * @return equivalent circuit
 *
 * @pre \p op is a multi-qubit operation
 *
 * @post The only multi-qubit gates in the replacement circuit are CX
 */
Circuit CX_circ_from_multiq(const Op_ptr op);

/**
 * Replace an operation with an equivalent circuit using CX, Rx and Rz
 *
 * @param op operation
 *
 * @return equivalent circuit
 */
Circuit CX_ZX_circ_from_op(const Op_ptr op);

/**
 * Replace CnRy, CnX, CnZ, CnY with 2-qubit gates and single qubit gates
 *
 * @param op operation
 * @param two_q_type whether rebase 2-q gates to CX or TK2
 *
 * @return equivalent circuit
 */
Circuit multi_controlled_to_2q(
    const Op_ptr op, const std::optional<OpType>& two_q_type = std::nullopt);

}  // namespace tket
