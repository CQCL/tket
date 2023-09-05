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

#include "tket/Ops/MetaOp.hpp"

#include <memory>
#include <typeinfo>

#include "tket/OpType/EdgeType.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/Utils/Json.hpp"

namespace tket {

MetaOp::MetaOp(OpType type, op_signature_t signature, const std::string& _data)
    : Op(type), signature_(signature), data_(_data) {
  if (!is_metaop_type(type)) throw BadOpType(type);
}

Op_ptr MetaOp::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return std::make_shared<MetaOp>(get_type(), get_signature(), get_data());
}

SymSet MetaOp::free_symbols() const { return {}; }

unsigned MetaOp::n_qubits() const {
  OptUInt n = desc_.n_qubits();
  if (n == any) {
    return std::count(signature_.begin(), signature_.end(), EdgeType::Quantum);
  } else {
    return n.value();
  }
}

op_signature_t MetaOp::get_signature() const {
  std::optional<op_signature_t> sig = desc_.signature();
  if (sig)
    return *sig;
  else
    return signature_;
}

nlohmann::json MetaOp::serialize() const {
  nlohmann::json j;
  j["type"] = get_type();
  j["signature"] = get_signature();
  j["data"] = get_data();
  return j;
}

Op_ptr MetaOp::deserialize(const nlohmann::json& j) {
  OpType optype = j.at("type").get<OpType>();
  op_signature_t sig = j.at("signature").get<op_signature_t>();
  std::string data;
  try {
    data = j.at("data").get<std::string>();
  } catch (const nlohmann::json::out_of_range& e) {
    data = "";
  }
  return std::make_shared<MetaOp>(optype, sig, data);
}

bool MetaOp::is_clifford() const { return true; }

MetaOp::~MetaOp() {}

bool MetaOp::is_equal(const Op& op_other) const {
  const MetaOp& other = dynamic_cast<const MetaOp&>(op_other);
  return (get_signature() == other.get_signature());
}

MetaOp::MetaOp() : Op(OpType::Barrier) {}

}  // namespace tket
