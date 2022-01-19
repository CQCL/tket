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

#pragma once

#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Utils/Json.hpp"

/**
 * When an OpType needs custom JSON conversion methods (as is the case for box
 * types), it should be implemented as a dedicated class and implement its own
 * @p from_json and @p to_json methods, and register them using this macro.
 */
#define REGISTER_OPFACTORY(type, opclass)                     \
  bool registered_##opclass = OpJsonFactory::register_method( \
      OpType::type, opclass::from_json, opclass::to_json);

namespace tket {

class OpJsonFactory {
 public:
  using JsonConstruct = Op_ptr (*)(const nlohmann::json &j);
  using JsonProduce = nlohmann::json (*)(const Op_ptr &j);

  OpJsonFactory() = delete;

  // register conversion methods to type
  static bool register_method(
      const OpType &type, JsonConstruct create_method,
      JsonProduce produce_method);

  // abstract interfaces for converting op types
  static Op_ptr from_json(const nlohmann::json &j);
  static nlohmann::json to_json(const Op_ptr &op);

 private:
  // map from optype to method which can produce JSON for that Op
  static std::map<OpType, JsonConstruct> &c_methods_();
  // map from optype to method which can construct that Op from JSON
  static std::map<OpType, JsonProduce> &p_methods_();
};

}  // namespace tket
