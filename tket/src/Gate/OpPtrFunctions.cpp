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

#include "OpPtrFunctions.hpp"

#include "Gate.hpp"
#include "Ops/MetaOp.hpp"
#include "SymTable.hpp"

namespace tket {

Op_ptr get_op_ptr(OpType chosen_type, const Expr& param, unsigned n_qubits) {
  return get_op_ptr(chosen_type, std::vector<Expr>{param}, n_qubits);
}

Op_ptr get_op_ptr(
    OpType chosen_type, const std::vector<Expr>& params, unsigned n_qubits) {
  if (is_gate_type(chosen_type)) {
    SymTable::register_symbols(expr_free_symbols(params));
    return std::make_shared<const Gate>(chosen_type, params, n_qubits);
  } else {
    return std::make_shared<const MetaOp>(chosen_type);
  }
}

}  // namespace tket
