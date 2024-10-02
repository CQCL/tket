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

#include "tket/Ops/BarrierOp.hpp"

#include "tket/OpType/EdgeType.hpp"
#include "tket/OpType/OpType.hpp"

namespace tket {

BarrierOp::BarrierOp(op_signature_t signature, const std::string& _data)
    : Op(OpType::Barrier), signature_(signature), data_(_data) {}

Op_ptr BarrierOp::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return Op_ptr();
}

SymSet BarrierOp::free_symbols() const { return {}; }

unsigned BarrierOp::n_qubits() const {
  return std::count(signature_.begin(), signature_.end(), EdgeType::Quantum);
}

op_signature_t BarrierOp::get_signature() const { return signature_; }

nlohmann::json BarrierOp::serialize() const {
  nlohmann::json j;
  j["type"] = get_type();
  j["signature"] = get_signature();
  j["data"] = get_data();
  return j;
}

Op_ptr BarrierOp::deserialize(const nlohmann::json& j) {
  op_signature_t sig = j.at("signature").get<op_signature_t>();
  std::string data;
  try {
    data = j.at("data").get<std::string>();
  } catch (const nlohmann::json::out_of_range& e) {
    data = "";
  }
  return std::make_shared<BarrierOp>(sig, data);
}

bool BarrierOp::is_clifford() const { return true; }

BarrierOp::~BarrierOp() {}

bool BarrierOp::is_equal(const Op& op_other) const {
  const BarrierOp& other = dynamic_cast<const BarrierOp&>(op_other);
  return (
      get_signature() == other.get_signature() &&
      get_data() == other.get_data());
}

}  // namespace tket
