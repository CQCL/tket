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

#include "tkwsm/Searching/WeightCalculator.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

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

  const std::vector<std::pair<VertexWSM, VertexWSM>>& assignments =
      accessor.get_new_assignments();

  for (auto ii = number_of_processed_assignments; ii < assignments.size();
       ++ii) {
    const VertexWSM& pv = assignments[ii].first;
    const VertexWSM& tv = assignments[ii].second;

    TKET_ASSERT(m_p_vertices_seen.insert(pv).second);

    // Look for assigned neighbours.
    for (const std::pair<VertexWSM, WeightWSM>& entry :
         pattern_ndata.get_neighbours_and_weights(pv)) {
      const VertexWSM& other_pv = entry.first;

      const auto& other_domain = accessor.get_domain(other_pv);
      const auto first_tv = other_domain.find_first();
      if (first_tv >= other_domain.size()) {
        // Empty domain; a nogood!
        return {};
      }
      const auto second_tv = other_domain.find_next(first_tv);
      if (second_tv >= other_domain.size() &&
          m_p_vertices_seen.count(other_pv) == 0) {
        // The domain has size exactly 1,
        // so we have an assigned edge.
        // But we also haven't seen both vertices yet,
        // so the edge cannot have been added already.
        const auto t_edge_weight_opt =
            target_ndata.get_edge_weight_opt(tv, first_tv);
        if (!t_edge_weight_opt) {
          return {};
        }
        result.scalar_product += entry.second * t_edge_weight_opt.value();
        if (result.scalar_product > max_scalar_product) {
          return {};
        }
        result.total_extra_p_edge_weights += entry.second;
      }
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
