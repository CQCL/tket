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

bool LexiLabellingMethod::check_method(
    const std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  std::shared_ptr<unit_frontier_t> frontier_edges =
      frontier_convert_vertport_to_edge(
          mapping_frontier->circuit_, mapping_frontier->quantum_boundary);
  CutFrontier next_cut = mapping_frontier->circuit_.next_q_cut(frontier_edges);

  for (const Vertex& vert : *next_cut.slice) {
    EdgeVec ev = mapping_frontier->circuit_.get_in_edges_of_type(
        vert, EdgeType::Quantum);
    // lexilabelling can't support dynamic labelling of >2 qubit gates
    if (ev.size() > 2) {
      return false;
    }
    for (const Edge& e : ev) {
      for (const std::pair<UnitID, Edge>& pair :
           frontier_edges->get<TagKey>()) {
        if (pair.second == e) {
          if (!architecture->node_exists(Node(pair.first))) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

unit_map_t LexiLabellingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  LexiRoute lr(architecture, mapping_frontier);
  lr.solve_labelling();
  return {};
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
