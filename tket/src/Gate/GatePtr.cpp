// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "GatePtr.hpp"

#include <memory>

#include "Gate.hpp"
#include "OpType/OpTypeInfo.hpp"
#include "Ops/OpPtr.hpp"

namespace tket {

Gate_ptr as_gate_ptr(Op_ptr op) {
  Gate_ptr gp = std::dynamic_pointer_cast<const Gate>(op);
  if (!gp) {
    throw BadOpType("Operation is not a gate", op->get_type());
  }
  return gp;
}

}  // namespace tket
