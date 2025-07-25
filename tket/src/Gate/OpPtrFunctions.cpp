// Copyright Quantinuum
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

#include "tket/Gate/OpPtrFunctions.hpp"

#include <mutex>

#include "tket/Gate/Gate.hpp"
#include "tket/Gate/SymTable.hpp"
#include "tket/OpType/OpTypeFunctions.hpp"
#include "tket/Ops/BarrierOp.hpp"
#include "tket/Ops/MetaOp.hpp"

namespace tket {

// Mutex to protect the global symbol table.
static std::mutex mtx;

Op_ptr get_op_ptr(OpType chosen_type, const Expr& param, unsigned n_qubits) {
  return get_op_ptr(chosen_type, std::vector<Expr>{param}, n_qubits);
}

Op_ptr get_op_ptr(
    OpType chosen_type, const std::vector<Expr>& params, unsigned n_qubits) {
  std::lock_guard<std::mutex> lock(mtx);
  if (is_gate_type(chosen_type)) {
    SymTable::register_symbols(expr_free_symbols(params));
    return std::make_shared<const Gate>(chosen_type, params, n_qubits);
  } else if (is_barrier_type(chosen_type)) {
    return std::make_shared<const BarrierOp>();
  } else {
    return std::make_shared<const MetaOp>(chosen_type);
  }
}

}  // namespace tket
