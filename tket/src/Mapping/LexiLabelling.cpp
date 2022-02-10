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
  std::set<Node> already_checked;
  // *it = {UnitID, {Vertex, Port}}
  for (auto it = mapping_frontier->quantum_boundary->get<TagKey>().begin();
       it != mapping_frontier->quantum_boundary->get<TagKey>().end(); ++it) {
    Edge e0 = mapping_frontier->circuit_.get_nth_out_edge(
        it->second.first, it->second.second);
    Vertex v0 = mapping_frontier->circuit_.target(e0);
    Node node = Node(it->first);
    // i.e. skip already checked vertices
    if (already_checked.find(node) == already_checked.end()) {
      already_checked.insert(node);
      // for counter<n_edges check, to avoid iterating through many many qubits
      int n_edges =
          mapping_frontier->circuit_.n_in_edges_of_type(v0, EdgeType::Quantum);
      int counter = 1;  // 1 edge
      auto jt = it;
      ++jt;
      while (jt != mapping_frontier->quantum_boundary->get<TagKey>().end() &&
             counter < n_edges) {
        Edge e1 = mapping_frontier->circuit_.get_nth_out_edge(
            jt->second.first, jt->second.second);
        Vertex v1 = mapping_frontier->circuit_.target(e1);
        if (v0 == v1) {
          counter++;
          // confirms that there is at least one multi-qubit gate in the first
          // layer which is assigned to some Qubit not in the architecture
          if (!architecture->node_exists(Node(it->first)) ||
              !architecture->node_exists(Node(jt->first))) {
            return true;
          }
        }
        ++jt;
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
