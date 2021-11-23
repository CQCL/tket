// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Ops/OpPtr.hpp"
#include "Utils/Json.hpp"

namespace tket {

std::ostream& operator<<(std::ostream& os, Op const& operation) {
  return os << operation.get_name();
}

void to_json(nlohmann::json& j, const Op_ptr& op) { j = op->serialize(); }

}  // namespace tket
