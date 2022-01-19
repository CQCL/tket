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

#include "Conditional.hpp"

#include "OpType/OpType.hpp"

namespace tket {

Conditional::Conditional(const Op_ptr& op, unsigned width, unsigned value)
    : Op(OpType::Conditional), op_(op), width_(width), value_(value) {}

Conditional::Conditional(const Conditional& other)
    : Op(other), op_(other.op_), width_(other.width_), value_(other.value_) {}

Conditional::~Conditional() {}

Op_ptr Conditional::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  return std::make_shared<Conditional>(
      op_->symbol_substitution(sub_map), width_, value_);
}

SymSet Conditional::free_symbols() const { return op_->free_symbols(); }

bool Conditional::is_equal(const Op& op_other) const {
  const Conditional& other = dynamic_cast<const Conditional&>(op_other);

  return *op_ == *other.get_op() && width_ == other.get_width() &&
         value_ == other.get_value();
}

unsigned Conditional::n_qubits() const { return op_->n_qubits(); }

op_signature_t Conditional::get_signature() const {
  op_signature_t signature(width_, EdgeType::Boolean);
  op_signature_t inner_sig = op_->get_signature();
  signature.insert(signature.end(), inner_sig.begin(), inner_sig.end());
  return signature;
}

nlohmann::json Conditional::serialize() const {
  nlohmann::json j;
  nlohmann::json j_cond;
  j_cond["op"] = get_op();
  j_cond["width"] = get_width();
  j_cond["value"] = get_value();
  j["type"] = OpType::Conditional;
  j["conditional"] = j_cond;
  return j;
}

Op_ptr Conditional::deserialize(const nlohmann::json& j) {
  nlohmann::json j_cond = j.at("conditional");
  Op_ptr cond_op = j_cond.at("op").get<Op_ptr>();
  return std::make_shared<Conditional>(
      cond_op, j_cond.at("width").get<unsigned>(),
      j_cond.at("value").get<unsigned>());
}

std::string Conditional::get_command_str(const unit_vector_t& args) const {
  std::stringstream out;
  out << "IF ([";
  if (width_ > 0) {
    out << args.at(0).repr();
    for (unsigned i = 1; i < width_; ++i) {
      out << ", " << args.at(i).repr();
    }
  }
  out << "] == " << value_;
  unit_vector_t inner_args(args.begin() + width_, args.end());
  out << ") THEN " << op_->get_command_str(inner_args);
  return out.str();
}

Op_ptr Conditional::dagger() const {
  const Op_ptr inner_dagger = op_->dagger();
  return std::make_shared<Conditional>(inner_dagger, width_, value_);
}

Op_ptr Conditional::get_op() const { return op_; }

unsigned Conditional::get_width() const { return width_; }

unsigned Conditional::get_value() const { return value_; }

Conditional::Conditional()
    : Op(OpType::Conditional), op_(), width_(0), value_(0) {}

}  // namespace tket
