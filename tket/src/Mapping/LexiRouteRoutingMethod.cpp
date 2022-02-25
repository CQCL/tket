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

#include "Mapping/LexiRouteRoutingMethod.hpp"

namespace tket{

LexiRouteRoutingMethod::LexiRouteRoutingMethod(unsigned _max_depth)
: max_depth_(_max_depth){};

    
std::pair<bool, unit_map_t> LexiRouteRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  LexiRoute lr(architecture, mapping_frontier);
  return {lr.solve(this->max_depth_), {}};
}

unsigned LexiRouteRoutingMethod::get_max_depth() const {
  return this->max_depth_;
}

nlohmann::json LexiRouteRoutingMethod::serialize() const {
  nlohmann::json j;
  j["depth"] = this->get_max_depth();
  j["name"] = "LexiRouteRoutingMethod";
  return j;
}

LexiRouteRoutingMethod LexiRouteRoutingMethod::deserialize(
    const nlohmann::json& j) {
  return LexiRouteRoutingMethod(j.at("depth").get<unsigned>());
}
    
}  // namespace tket


// bool LexiRouteRoutingMethod::check_method(
//     const std::shared_ptr<MappingFrontier>& mapping_frontier,
//     const ArchitecturePtr& architecture) const {
//   std::shared_ptr<unit_frontier_t> frontier_edges =
//       frontier_convert_vertport_to_edge(
//           mapping_frontier->circuit_, mapping_frontier->quantum_boundary);
//   CutFrontier next_cut = mapping_frontier->circuit_.next_q_cut(frontier_edges);
//   for (const Vertex& vert : *next_cut.slice) {
//     Op_ptr op = mapping_frontier->circuit_.get_Op_ptr_from_Vertex(vert);
//     // can't work wih box ops, or gates with more than 2 qubits that aren't a
//     // BRIDGE

//     if ((mapping_frontier->circuit_.n_in_edges_of_type(
//              vert, EdgeType::Quantum) > 2 &&
//          op->get_type() != OpType::BRIDGE) ||
//         (op->get_desc().is_box() || (op->get_type() == OpType::Conditional &&
//                                      static_cast<const Conditional&>(*op)
//                                          .get_op()
//                                          ->get_desc()
//                                          .is_box()))) {
//       return false;
//     } else {
//       // second check that all input UnitID are actually in architecture
//       for (const Edge& e : mapping_frontier->circuit_.get_in_edges_of_type(
//                vert, EdgeType::Quantum)) {
//         for (const std::pair<UnitID, Edge>& pair :
//              frontier_edges->get<TagKey>()) {
//           if (pair.second == e) {
//             if (!architecture->node_exists(Node(pair.first))) {
//               return false;
//             }
//           }
//         }
//       }
//     }
//   }
//   return true;
// }