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

#include "tkwsm/Reducing/DistancesReducer.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Common/SetIntersection.hpp"
#include "tkwsm/GraphTheoretic/FilterUtils.hpp"
#include "tkwsm/GraphTheoretic/NearNeighboursData.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

#include "WeightSubgrMono/Common/TemporaryRefactorCode.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DistancesReducer::DistancesReducer(
    NearNeighboursData& pattern_near_ndata, const NeighboursData& target_ndata,
    NearNeighboursData& target_near_ndata, unsigned distance)
    : m_pattern_near_ndata(pattern_near_ndata),
      m_target_ndata(target_ndata),
      m_target_near_ndata(target_near_ndata),
      m_distance(distance) {
  TKET_ASSERT(m_distance > 1);
}

bool DistancesReducer::check(std::pair<VertexWSM, VertexWSM> assignment) {
  return FilterUtils::compatible_sorted_degree_counts(
      m_pattern_near_ndata.get_degree_counts(assignment.first, m_distance),
      m_target_near_ndata.get_degree_counts(assignment.second, m_distance));
}


ReductionResult DistancesReducer::reduce(
    std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
    std::set<VertexWSM>&) {
  const auto& pv_at_distance_d = m_pattern_near_ndata.get_vertices_at_distance(
      assignment.first, m_distance);

  if (pv_at_distance_d.empty()) {
    // Nothing to check!
    return ReductionResult::SUCCESS;
  }

  TemporaryRefactorCode();
  boost::dynamic_bitset<> work_bitset;
  
  const std::vector<std::pair<VertexWSM, WeightWSM>>& tv_neighbours =
      m_target_ndata.get_neighbours_and_weights(assignment.second);

  auto result = ReductionResult::SUCCESS;

  for (VertexWSM new_pv : pv_at_distance_d) {
    {
      const auto& current_domain = accessor.get_domain(assignment.first);
      if (other_vertex_reduction_can_be_skipped_by_symmetry(
              current_domain, accessor, assignment.first, new_pv)) {
        continue;
      }
      work_bitset.resize(current_domain.size());
    }
    work_bitset.reset();
    // This is for d=1, obviously a special case.
    for(const auto& entry : tv_neighbours) {
      // TODO: store bitsets!
      TKET_ASSERT(entry.first < work_bitset.size());
      TKET_ASSERT(!work_bitset.test_set(entry.first));
    }

    for (unsigned jj = 2; jj <= m_distance; ++jj) {
      const auto& new_tv_list =
          m_target_near_ndata.get_vertices_at_distance(assignment.second, jj);
      
      if(new_tv_list.empty()) {
        break;
      }
      for(auto tv : new_tv_list) {
        TKET_ASSERT(tv<work_bitset.size());
        TKET_ASSERT(!work_bitset.test_set(tv));
      }
    }

    switch (accessor.intersect_domain_with_swap(new_pv, work_bitset).reduction_result) {
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
