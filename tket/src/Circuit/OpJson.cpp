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

#include "Boxes.hpp"
#include "Conditional.hpp"
#include "Gate/Gate.hpp"
#include "OpType/OpType.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Ops/ClassicalOps.hpp"
#include "Ops/MetaOp.hpp"
#include "Ops/OpPtr.hpp"
#include "Utils/Json.hpp"

namespace tket {

void from_json(const nlohmann::json& j, Op_ptr& op) {
  OpType optype = j.at("type").get<OpType>();
  if (is_metaop_type(optype)) {
    op = MetaOp::deserialize(j);
  } else if (is_box_type(optype)) {
    op = Box::deserialize(j);
  } else if (optype == OpType::Conditional) {
    op = Conditional::deserialize(j);
  } else if (optype == OpType::WASM) {
    op = WASMOp::deserialize(j);
  } else if (is_classical_type(optype)) {
    op = ClassicalOp::deserialize(j);
  } else if (is_gate_type(optype)) {
    op = Gate::deserialize(j);
  } else {
    throw JsonError(
        "Deserialization not yet implemented for " +
        optypeinfo().at(optype).name);
  }
}

}  // namespace tket
