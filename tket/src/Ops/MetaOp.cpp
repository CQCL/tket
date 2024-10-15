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

#include "tket/Ops/MetaOp.hpp"

#include <memory>

#include "tket/OpType/EdgeType.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/OpType/OpTypeFunctions.hpp"

namespace tket {

MetaOp::MetaOp(OpType type, op_signature_t signature, const std::string& _data)
    : Op(type), signature_(signature), data_(_data) {
  if (!is_metaop_type(type)) throw BadOpType(type);
}

Op_ptr MetaOp::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return std::make_shared<MetaOp>(*this);
}

SymSet MetaOp::free_symbols() const { return {}; }

op_signature_t MetaOp::get_signature() const {
  std::optional<op_signature_t> sig = desc_.signature();
  if (sig)
    return *sig;
  else
    return signature_;
}

bool MetaOp::is_clifford() const { return false; }

MetaOp::~MetaOp() {}

}  // namespace tket
