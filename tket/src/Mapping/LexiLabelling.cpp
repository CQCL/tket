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
#include "Mapping/LexiLabelling.hpp"

namespace tket {

std::pair<bool, unit_map_t> LexiLabellingMethod::routing_method(
    MappingFrontier_ptr& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  LexiRoute lr(architecture, mapping_frontier);
  return {lr.solve_labelling(), {}};
}

nlohmann::json LexiLabellingMethod::serialize() const {
  nlohmann::json j;
  j["name"] = "LexiLabellingMethod";
  return j;
}

LexiLabellingMethod LexiLabellingMethod::deserialize(
    const nlohmann::json& /*j*/) {
  return LexiLabellingMethod();
}

}  // namespace tket
