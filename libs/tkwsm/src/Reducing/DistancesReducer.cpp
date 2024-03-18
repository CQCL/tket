// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/Reducing/DistancesReducer.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/FilterUtils.hpp"
#include "tkwsm/GraphTheoretic/NearNeighboursData.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DistancesReducer::DistancesReducer(
    NearNeighboursData& pattern_near_ndata,
    NearNeighboursData& target_near_ndata, unsigned distance)
    : m_pattern_near_ndata(pattern_near_ndata),
      m_target_near_ndata(target_near_ndata),
      m_distance(distance) {
  TKET_ASSERT(m_distance > 0);
}

bool DistancesReducer::check(std::pair<VertexWSM, VertexWSM> assignment) {
  return FilterUtils::compatible_sorted_degree_counts(
      m_pattern_near_ndata.get_degree_counts_at_exact_distance(
          assignment.first, m_distance),
      m_target_near_ndata.get_degree_counts_up_to_distance(
          assignment.second, m_distance));
}

ReductionResult DistancesReducer::reduce(
    std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
    boost::dynamic_bitset<>& work_bitset) {
  auto result = ReductionResult::SUCCESS;

  const boost::dynamic_bitset<>& pv_at_distance_d =
      m_pattern_near_ndata.get_vertices_at_exact_distance(
          assignment.first, m_distance);

  for (auto new_pv = pv_at_distance_d.find_first();
       new_pv < pv_at_distance_d.size();
       new_pv = pv_at_distance_d.find_next(new_pv)) {
    work_bitset = m_target_near_ndata.get_vertices_up_to_distance(
        assignment.second, m_distance);

    switch (accessor.intersect_domain_with_swap(new_pv, work_bitset)
                .reduction_result) {
      case ReductionResult::NEW_ASSIGNMENTS:
        // Normally we would return now,
        // BUT the assignment PV->TV hasn't been fully processed
        // until ALL PV' with  dist(PV, PV')=d  have been considered.
        // Thus, continue, to finish off this assignment,
        // rather than trying to store partial data
        // and resuming later.
        result = ReductionResult::NEW_ASSIGNMENTS;
        break;
      case ReductionResult::NOGOOD:
        return ReductionResult::NOGOOD;
      case ReductionResult::SUCCESS:
        break;
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
