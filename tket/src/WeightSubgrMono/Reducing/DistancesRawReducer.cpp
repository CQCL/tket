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

#include "WeightSubgrMono/Reducing/DistancesRawReducer.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/GraphTheoretic/NearNeighboursData.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DistancesRawReducer::DistancesRawReducer(
    const NeighboursData& pattern_ndata, NearNeighboursData& pattern_near_ndata,
    const NeighboursData& target_ndata, NearNeighboursData& target_near_ndata)
    : m_pattern_ndata(pattern_ndata),
      m_pattern_near_ndata(pattern_near_ndata),
      m_target_ndata(target_ndata),
      m_target_near_ndata(target_near_ndata) {}

bool DistancesRawReducer::check(
    const std::pair<VertexWSM, VertexWSM>& assignment,
    unsigned distance) const {
  // We must have   #{u : dist(pv,u) = j } <= #{v : dist(tv,v) <= j}
  // for each j; but ALSO
  //  #{u : dist(pv,u) <= j } <= #{v : dist(tv,v) <= j }
  // The second clearly implies the first.
  auto pv_count_at_lower_distances =
      m_pattern_ndata.get_degree(assignment.first);
  auto tv_count_at_lower_distances =
      m_target_ndata.get_degree(assignment.second);
  if (pv_count_at_lower_distances > tv_count_at_lower_distances) {
    return false;
  }
  for (unsigned ii = 2; ii <= distance; ++ii) {
    pv_count_at_lower_distances +=
        m_pattern_near_ndata.get_vertices_at_distance(assignment.first, ii)
            .size();

    tv_count_at_lower_distances +=
        m_target_near_ndata.get_vertices_at_distance(assignment.second, ii)
            .size();
    if (pv_count_at_lower_distances > tv_count_at_lower_distances) {
      return false;
    }
  }
  return true;
}

DistancesRawReducer::Result DistancesRawReducer::reduce_neighbours(
    const std::pair<VertexWSM, VertexWSM>& assignment, NodeWSM& node) {
  const auto& pv_neighbours =
      m_pattern_ndata.get_neighbours_and_weights(assignment.first);
  const auto& tv_neighbours =
      m_target_ndata.get_neighbours_and_weights(assignment.second);
  if (pv_neighbours.size() > tv_neighbours.size()) {
    return Result::IMPOSSIBLE_ASSIGNMENT;
  }
  for (const auto& pv_entry : pv_neighbours) {
    const VertexWSM& pv_other = pv_entry.first;
    const auto& current_domain = node.get_possible_assignments().at(pv_other);
    fill_intersection_ignoring_second_elements(
        current_domain, tv_neighbours, m_work_set);

    if (current_domain.size() == m_work_set.size()) {
      // No change.
      continue;
    }
    if (m_work_set.empty()) {
      return Result::IMPOSSIBLE_NODE;
    }
    node.overwrite_domain(pv_other, m_work_set);
  }
  return Result::SUCCESS;
}

DistancesRawReducer::Result DistancesRawReducer::operator()(
    const std::pair<VertexWSM, VertexWSM>& assignment, NodeWSM& node,
    unsigned distance) {
  TKET_ASSERT(distance > 0);
  if (distance == 1) {
    return reduce_neighbours(assignment, node);
  }
  const auto& pv_at_distance_d =
      m_pattern_near_ndata.get_vertices_at_distance(assignment.first, distance);

  if (pv_at_distance_d.empty()) {
    // Nothing to check!
    return Result::SUCCESS;
  }

  // Now, for each PV' with  dist(PV, PV')=d,  we must intersect Dom(PV')
  // with the disjoint UNION of  A(1), A(2), ..., A(d),
  // where A(j) = { u : dist(tv, u)=j }.
  // To save space, we don't explicitly construct this union,
  // which could be large; we just build it up by combining
  //    Dom(PV') intersect A(j)   (which are disjoint) directly.
  for (VertexWSM new_pv : pv_at_distance_d) {
    const auto& current_domain = node.get_possible_assignments().at(new_pv);
    const auto& tv_neighbours =
        m_target_ndata.get_neighbours_and_weights(assignment.second);
    fill_intersection_ignoring_second_elements(
        current_domain, tv_neighbours, m_work_set);
    m_work_vector = {m_work_set.cbegin(), m_work_set.cend()};

    // Now, add all the other intersected sets.
    for (unsigned jj = 2; jj <= distance; ++jj) {
      const auto& new_tv_list =
          m_target_near_ndata.get_vertices_at_distance(assignment.second, jj);

      fill_intersection(current_domain, new_tv_list, m_work_set);

      // Append: the different sets are automatically disjoint, of course.
      for (VertexWSM new_tv : m_work_set) {
        m_work_vector.push_back(new_tv);
      }
    }
    // Now "m_work_vector" is the final intersection,
    // although NOT sorted; but that doesn't matter.
    if (m_work_vector.empty()) {
      return Result::IMPOSSIBLE_NODE;
    }
    node.overwrite_domain(new_pv, m_work_vector);
  }
  return Result::SUCCESS;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
