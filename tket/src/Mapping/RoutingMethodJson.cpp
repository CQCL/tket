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

#include "Mapping/RoutingMethodJson.hpp"

#include "Mapping/LexiLabelling.hpp"

namespace tket {

void to_json(nlohmann::json& j, const RoutingMethod& rm) { j = rm.serialize(); }

void from_json(const nlohmann::json& /*j*/, RoutingMethod& rm) {
  rm = RoutingMethod();
}

void to_json(nlohmann::json& j, const std::vector<RoutingMethodPtr>& rmp_v) {
  for (const auto& r : rmp_v) {
    j.push_back(*r);
  }
}

void from_json(const nlohmann::json& j, std::vector<RoutingMethodPtr>& rmp_v) {
  for (const auto& c : j) {
    std::string name = c.at("name").get<std::string>();
    if (name == "LexiLabellingMethod") {
      rmp_v.push_back(std::make_shared<LexiLabellingMethod>(
          LexiLabellingMethod::deserialize(c)));
    } else if (name == "LexiRouteRoutingMethod") {
      rmp_v.push_back(std::make_shared<LexiRouteRoutingMethod>(
          LexiRouteRoutingMethod::deserialize(c)));
    } else if (name == "RoutingMethod") {
      rmp_v.push_back(std::make_shared<RoutingMethod>());
    } else if (name == "AASRouteRoutingMethod") {
      rmp_v.push_back(std::make_shared<AASRouteRoutingMethod>(
          AASRouteRoutingMethod::deserialize(c)));
    } else if (name == "AASLabellingMethod") {
      rmp_v.push_back(std::make_shared<AASLabellingMethod>(
          AASLabellingMethod::deserialize(c)));
    } else if (name == "MultiGateReorderRoutingMethod") {
      rmp_v.push_back(std::make_shared<MultiGateReorderRoutingMethod>(
          MultiGateReorderRoutingMethod::deserialize(c)));
    } else if (name == "BoxDecompositionRoutingMethod") {
      rmp_v.push_back(std::make_shared<BoxDecompositionRoutingMethod>(
          BoxDecompositionRoutingMethod::deserialize(c)));
    } else {
      std::logic_error(
          "Deserialization for given RoutingMethod not supported.");
    }
  }
}

}  // namespace tket
