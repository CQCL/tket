// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <tkassert/Assert.hpp>

#include "tket/Circuit/Command.hpp"
#include "tket/Utils/Json.hpp"

namespace tket {

void to_json(nlohmann::json& j, const Command& com) {
  const Op_ptr op = com.get_op_ptr();
  const std::optional<std::string> opgroup = com.get_opgroup();
  j["op"] = op;
  if (opgroup) {
    j["opgroup"] = opgroup.value();
  }

  const op_signature_t& sig = op->get_signature();
  const unit_vector_t& args = com.get_args();

  if (sig.size() != args.size()) {
    JsonError("Number of args does not match signature of op.");
  }

  nlohmann::json j_args;

  for (std::size_t i = 0; i < sig.size(); i++) {
    switch (sig[i]) {
      case EdgeType::WASM: {
        const WasmState& wb = static_cast<const WasmState&>(args[i]);
        j_args.push_back(wb);
        break;
      }
      case EdgeType::Quantum: {
        const Qubit& qb = static_cast<const Qubit&>(args[i]);
        j_args.push_back(qb);
        break;
      }
      case EdgeType::Classical:
      case EdgeType::Boolean: {
        const Bit& cb = static_cast<const Bit&>(args[i]);
        j_args.push_back(cb);
        break;
      }
      default: {
        TKET_ASSERT(!"command to json found invalid edge type in signature");
      }
    }
  }

  j["args"] = j_args;
}
void from_json(const nlohmann::json& j, Command& com) {
  const auto op = j.at("op").get<Op_ptr>();
  std::optional<std::string> opgroup;
  if (j.contains("opgroup")) {
    opgroup = j.at("opgroup").get<std::string>();
  }
  const op_signature_t& sig = op->get_signature();

  const nlohmann::json& j_args = j.at("args");
  if (sig.size() != j_args.size()) {
    JsonError("Number of args does not match signature of op.");
  }
  unit_vector_t args;
  for (std::size_t i = 0; i < sig.size(); i++) {
    switch (sig[i]) {
      case EdgeType::WASM: {
        args.push_back(j_args[i].get<WasmState>());
        break;
      }
      case EdgeType::Quantum: {
        args.push_back(j_args[i].get<Qubit>());
        break;
      }
      case EdgeType::Classical:
      case EdgeType::Boolean: {
        args.push_back(j_args[i].get<Bit>());
        break;
      }
      default: {
        TKET_ASSERT(!"command from json found invalid edge type in signature");
      }
    }
  }
  com = Command(op, args, opgroup);
}

}  // namespace tket
