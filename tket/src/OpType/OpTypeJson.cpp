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

#include <map>
#include <string>

#include "OpType.hpp"
#include "OpTypeInfo.hpp"

namespace tket {

// map from OpType name to OpType
// Relies on unique OpType names.
const std::map<std::string, OpType>& name_to_optype() {
  auto mapfiller = []() {
    std::map<std::string, OpType> fill_map;
    for (const auto& info : optypeinfo()) {
      fill_map.insert({info.second.name, info.first});
    }
    return fill_map;
  };
  static auto final_map =
      std::make_unique<const std::map<std::string, OpType>>(mapfiller());
  return *final_map;
}

void to_json(nlohmann::json& j, const OpType& type) {
  j = optypeinfo().at(type).name;
}

void from_json(const nlohmann::json& j, OpType& type) {
  const std::string name = j.get<std::string>();
  const auto result = name_to_optype().find(name);
  if (result == name_to_optype().end()) {
    throw JsonError("No OpType with name " + name);
  }

  type = result->second;
}

}  // namespace tket
