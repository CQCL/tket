// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Circuit/ResourceData.hpp"

#include "tket/Utils/Json.hpp"

namespace tket {

bool ResourceData::operator==(const ResourceData& other) const {
  return OpTypeCount == other.OpTypeCount && GateDepth == other.GateDepth &&
         OpTypeDepth == other.OpTypeDepth &&
         TwoQubitGateDepth == other.TwoQubitGateDepth;
}

void to_json(nlohmann::json& j, const ResourceData& data) {
  j["op_type_count"] = data.OpTypeCount;
  j["gate_depth"] = data.GateDepth;
  j["op_type_depth"] = data.OpTypeDepth;
  j["two_qubit_gate_depth"] = data.TwoQubitGateDepth;
}

void from_json(const nlohmann::json& j, ResourceData& data) {
  data.OpTypeCount =
      j.at("op_type_count").get<std::map<OpType, ResourceBounds<unsigned>>>();
  data.GateDepth = j.at("gate_depth").get<ResourceBounds<unsigned>>();
  data.OpTypeDepth =
      j.at("op_type_depth").get<std::map<OpType, ResourceBounds<unsigned>>>();
  data.TwoQubitGateDepth =
      j.at("two_qubit_gate_depth").get<ResourceBounds<unsigned>>();
}

}  // namespace tket
