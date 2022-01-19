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

#include "OpJsonFactory.hpp"

#include "OpType/OpTypeInfo.hpp"
#include "Ops/Op.hpp"
#include "Utils/Json.hpp"

namespace tket {
std::map<OpType, OpJsonFactory::JsonConstruct>& OpJsonFactory::c_methods_() {
  static std::unique_ptr<std::map<OpType, OpJsonFactory::JsonConstruct>>
      methods =
          std::make_unique<std::map<OpType, OpJsonFactory::JsonConstruct>>();
  return *methods;
}
std::map<OpType, OpJsonFactory::JsonProduce>& OpJsonFactory::p_methods_() {
  static std::unique_ptr<std::map<OpType, OpJsonFactory::JsonProduce>> methods =
      std::make_unique<std::map<OpType, OpJsonFactory::JsonProduce>>();
  return *methods;
}

bool OpJsonFactory::register_method(
    const OpType& type, JsonConstruct create_method,
    JsonProduce produce_method) {
  if (auto it = c_methods_().find(type); it == c_methods_().end()) {
    c_methods_()[type] = create_method;
    p_methods_()[type] = produce_method;
    return true;
  }
  return false;
}

Op_ptr OpJsonFactory::from_json(const nlohmann::json& j) {
  const auto type = j.at("type").get<OpType>();
  if (auto it = c_methods_().find(type); it != c_methods_().end()) {
    return it->second(j);
  }
  throw JsonError(
      "No from_json conversion for type " + optypeinfo().at(type).name);
}

nlohmann::json OpJsonFactory::to_json(const Op_ptr& op) {
  const OpType& type = op->get_type();
  if (auto it = p_methods_().find(type); it != p_methods_().end()) {
    return it->second(op);
  }
  throw JsonError(
      "No to_json conversion registered for type: " +
      optypeinfo().at(type).name);
}

}  // namespace tket
