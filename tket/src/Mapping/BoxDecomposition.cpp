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
#include "Mapping/BoxDecomposition.hpp"

#include "Mapping/MappingFrontier.hpp"

namespace tket {

BoxDecomposition::BoxDecomposition(
    const ArchitecturePtr &_architecture,
    MappingFrontier_ptr &_mapping_frontier)
    : architecture_(_architecture), mapping_frontier_(_mapping_frontier) {}

bool BoxDecomposition::solve() {
  // Box type vertices are later removed from DAG
  VertexList bin;
  bool modified = false;
  std::shared_ptr<unit_frontier_t> frontier_edges =
      frontier_convert_vertport_to_edge(
          this->mapping_frontier_->circuit_,
          this->mapping_frontier_->linear_boundary);
  CutFrontier next_cut =
      this->mapping_frontier_->circuit_.next_q_cut(frontier_edges);
  for (Vertex &vert : *next_cut.slice) {
    Op_ptr op = this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(vert);
    if (op->get_desc().is_box() ||
        (op->get_type() == OpType::Conditional &&
         static_cast<const Conditional &>(*op).get_op()->get_desc().is_box())) {
      if (this->mapping_frontier_->circuit_.substitute_box_vertex(
              vert, Circuit::VertexDeletion::No)) {
        modified = true;
        bin.push_back(vert);
      }
    }
  }
  if (!modified) {
    return false;
  }
  // Delete vertices
  this->mapping_frontier_->circuit_.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return true;
}

BoxDecompositionRoutingMethod::BoxDecompositionRoutingMethod(){};

std::pair<bool, unit_map_t> BoxDecompositionRoutingMethod::routing_method(
    MappingFrontier_ptr &mapping_frontier,
    const ArchitecturePtr &architecture) const {
  BoxDecomposition bd(architecture, mapping_frontier);
  bool modified = bd.solve();
  return {modified, {}};
}

nlohmann::json BoxDecompositionRoutingMethod::serialize() const {
  nlohmann::json j;
  j["name"] = "BoxDecompositionRoutingMethod";
  return j;
}

BoxDecompositionRoutingMethod BoxDecompositionRoutingMethod::deserialize(
    const nlohmann::json & /*j*/) {
  return BoxDecompositionRoutingMethod();
}

}  // namespace tket