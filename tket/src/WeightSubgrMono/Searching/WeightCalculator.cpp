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

#include "WeightSubgrMono/Searching/WeightCalculator.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "WeightSubgrMono/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

std::optional<WeightCalculator::Result> WeightCalculator::operator()(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    const DomainsAccessor& accessor,
    std::size_t number_of_processed_assignments,
    WeightWSM max_scalar_product) const {
  m_p_vertices_seen.clear();

  Result result;
  result.scalar_product = accessor.get_scalar_product();
  result.total_extra_p_edge_weights = 0;

  const auto& assignments = accessor.get_new_assignments();

  for (auto ii = number_of_processed_assignments; ii < assignments.size();
       ++ii) {
    const VertexWSM& pv = assignments[ii].first;
    const VertexWSM& tv = assignments[ii].second;

    TKET_ASSERT(m_p_vertices_seen.insert(pv).second);

    // Look for assigned neighbours.
    for (const auto& entry : pattern_ndata.get_neighbours_and_weights(pv)) {
      const VertexWSM& other_pv = entry.first;
      const auto& other_domain = accessor.get_domain(other_pv);
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
          if (result.scalar_product > max_scalar_product) {
            return {};
          }
          result.total_extra_p_edge_weights += entry.second;
          break;
        }
        default:
          break;
      }
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
