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

#include "WeightSubgrMono/Searching/WeightUpdater.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

std::optional<WeightUpdater::Result> WeightUpdater::operator()(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    const PossibleAssignments& possible_assignments,
    const std::vector<std::pair<VertexWSM, VertexWSM>>& assignments,
    std::size_t number_of_p_vertices_previously_processed_in_this_node,
    WeightWSM current_weight, WeightWSM max_weight,
    std::set<VertexWSM>& unassigned_neighbour_vertices) const {
  m_p_vertices_seen.clear();
  Result result;
  result.scalar_product = current_weight;
  result.total_extra_p_edge_weights = 0;

  for (auto ii = number_of_p_vertices_previously_processed_in_this_node;
       ii < assignments.size(); ++ii) {
    const VertexWSM& pv = assignments[ii].first;
    const VertexWSM& tv = assignments[ii].second;

    TKET_ASSERT(m_p_vertices_seen.insert(pv).second);

    // Look for assigned neighbours.
    for (const auto& entry : pattern_ndata.get_neighbours_and_weights(pv)) {
      const VertexWSM& other_pv = entry.first;
      const auto& other_domain = possible_assignments.at(other_pv);
      switch (other_domain.size()) {
        case 0:
          return {};
        case 1: {
          // We have an assigned edge.
          if (m_p_vertices_seen.count(other_pv) != 0) {
            // We've already seen both vertices,
            // so the edge must already have been added.
            break;
          }
          const VertexWSM& other_tv = *other_domain.cbegin();
          const auto t_edge_weight_opt =
              target_ndata.get_edge_weight_opt(tv, other_tv);
          if (!t_edge_weight_opt) {
            return {};
          }
          result.scalar_product += entry.second * t_edge_weight_opt.value();
          if (result.scalar_product > max_weight) {
            return {};
          }
          result.total_extra_p_edge_weights += entry.second;
        } break;
        default:
          // An unassigned vertex, which is ALSO adjacent to an assigned one.
          unassigned_neighbour_vertices.insert(other_pv);
      }
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
