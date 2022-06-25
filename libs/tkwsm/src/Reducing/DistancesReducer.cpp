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

#include "WeightSubgrMono/Reducing/DistancesReducer.hpp"

#include <tkassert/Assert.hpp>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NearNeighboursData.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "WeightSubgrMono/Searching/DomainsAccessor.hpp"

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

  // Ensure that we never need to resize again.
  m_work_vectors_list.resize(m_distance - 1);
}

bool DistancesReducer::check(std::pair<VertexWSM, VertexWSM> assignment) {
  return FilterUtils::compatible_sorted_degree_counts(
      m_pattern_near_ndata.get_degree_counts(assignment.first, m_distance),
      m_target_near_ndata.get_degree_counts(assignment.second, m_distance));
}

ReductionResult DistancesReducer::reduce(
    std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
    std::set<VertexWSM>& work_set) {
  const auto& pv_at_distance_d = m_pattern_near_ndata.get_vertices_at_distance(
      assignment.first, m_distance);

  if (pv_at_distance_d.empty()) {
    // Nothing to check!
    return ReductionResult::SUCCESS;
  }

  // To save space, we don't explicitly construct the union
  // {v : dist(TV,v) <= j}, which could be large.
  // For each PV' with  dist(PV, PV')=d,  we intersect Dom(PV')
  // with the disjoint UNION of  A(1), A(2), ..., A(d),
  // where A(j) = {v : dist(TV,v)=j}.
  // We build it up by combining
  //    Dom(PV') intersect A(j)   (which are disjoint), directly.
  const std::vector<std::pair<VertexWSM, WeightWSM>>& tv_neighbours =
      m_target_ndata.get_neighbours_and_weights(assignment.second);

  auto result = ReductionResult::SUCCESS;

  for (VertexWSM new_pv : pv_at_distance_d) {
    const std::set<VertexWSM>& current_domain = accessor.get_domain(new_pv);
    if (other_vertex_reduction_can_be_skipped_by_symmetry(
            current_domain, accessor, assignment.first, new_pv)) {
      continue;
    }
    // This is for d=1, obviously a special case.
    fill_intersection_ignoring_second_elements(
        current_domain, tv_neighbours, work_set);

    auto size_so_far = work_set.size();
    unsigned next_work_vector_index = 0;

    // Add all the other intersected sets.
    for (unsigned jj = 2; jj <= m_distance; ++jj) {
      const auto& new_tv_list =
          m_target_near_ndata.get_vertices_at_distance(assignment.second, jj);

      // The constructor already resized the vectors_list,
      // to ensure a valid index.
      std::vector<VertexWSM>& work_vector =
          m_work_vectors_list[next_work_vector_index];

      fill_intersection(current_domain, new_tv_list, work_vector);
      if (!work_vector.empty()) {
        // If it was empty, just reuse for the next distance.
        size_so_far += work_vector.size();
        ++next_work_vector_index;
      }
    }
    if (size_so_far == current_domain.size()) {
      // No change!
      continue;
    }

    if (size_so_far == 0) {
      return ReductionResult::NOGOOD;
    }

    // Now, take the union of all the disjoint sets.
    for (unsigned index = 0; index < next_work_vector_index; ++index) {
      for (VertexWSM tv : m_work_vectors_list[index]) {
        work_set.insert(tv);
      }
    }
    // Now "m_work_vector" is the final intersection,
    // although NOT sorted; but that doesn't matter.
    switch (accessor.overwrite_domain_with_set_swap(new_pv, work_set)) {
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
        // We should already have handled nogoods.
        TKET_ASSERT(false);
      case ReductionResult::SUCCESS:
        break;
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
