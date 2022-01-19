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

#include "Op.hpp"

#include <sstream>

#include "Ops/OpPtr.hpp"
#include "Utils/Json.hpp"

namespace tket {

std::string Op::get_name(bool latex) const {
  if (latex) {
    return get_desc().latex();
  } else {
    return get_desc().name();
  }
}

std::string Op::get_command_str(const unit_vector_t& args) const {
  std::stringstream out;
  out << get_name();
  if (!args.empty()) {
    out << " " << args[0].repr();
    for (unsigned i = 1; i < args.size(); i++) {
      out << ", " << args[i].repr();
    }
  }
  out << ";";
  return out.str();
}

std::ostream& operator<<(std::ostream& os, Op const& operation) {
  return os << operation.get_name();
}

void to_json(nlohmann::json& j, const Op_ptr& op) { j = op->serialize(); }

}  // namespace tket
