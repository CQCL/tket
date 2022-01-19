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

#include "FlowOp.hpp"

namespace tket {

FlowOp::FlowOp(OpType type, std::optional<std::string> label)
    : Op(type), label_(label) {
  if (!is_flowop_type(type)) {
    throw NotValid();
  }
}

Op_ptr FlowOp::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return Op_ptr();
}

SymSet FlowOp::free_symbols() const { return {}; }

std::optional<std::string> FlowOp::get_label() const { return label_; }

FlowOp::~FlowOp() {}

FlowOp::FlowOp() : Op(OpType::Stop) {}

bool FlowOp::is_equal(const Op& op_other) const {
  const FlowOp& other = dynamic_cast<const FlowOp&>(op_other);
  return (this->get_label() == other.get_label());
}

op_signature_t FlowOp::get_signature() const {
  std::optional<op_signature_t> sig = desc_.signature();
  if (sig)
    return *sig;
  else
    throw NotValid();
}

std::string FlowOp::get_name(bool latex) const {
  std::stringstream name;
  if (latex) {
    name << get_desc().latex() << "(";
  } else {
    name << get_desc().name();
  }
  if (type_ != OpType::Stop) {
    name << " " << *label_;
  }
  return name.str();
}

}  // namespace tket
