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

#include "WeightSubgrMono/Reducing/NeighboursReducer.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "WeightSubgrMono/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NeighboursReducer::NeighboursReducer(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata)
    : m_pattern_ndata(pattern_ndata), m_target_ndata(target_ndata) {}

bool NeighboursReducer::check(std::pair<VertexWSM, VertexWSM> assignment) {
  return m_pattern_ndata.get_degree(assignment.first) <=
         m_target_ndata.get_degree(assignment.second);
}

ReductionResult NeighboursReducer::reduce(
    std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
    std::set<VertexWSM>& work_set) {
  const std::vector<std::pair<VertexWSM, WeightWSM>>&
      target_neighbours_and_weights =
          m_target_ndata.get_neighbours_and_weights(assignment.second);
  auto result = ReductionResult::SUCCESS;

  for (const std::pair<VertexWSM, WeightWSM>& p_neighbour_entry :
       m_pattern_ndata.get_neighbours_and_weights(assignment.first)) {
    const VertexWSM& p_neighbour = p_neighbour_entry.first;
    const std::set<VertexWSM>& domain = accessor.get_domain(p_neighbour);
    if (other_vertex_reduction_can_be_skipped_by_symmetry(
            domain, accessor, assignment.first, p_neighbour)) {
      continue;
    }
    fill_intersection_ignoring_second_elements(
        domain, target_neighbours_and_weights, work_set);

    switch (accessor.overwrite_domain_with_set_swap(p_neighbour, work_set)) {
      case ReductionResult::SUCCESS:
        break;
      case ReductionResult::NOGOOD:
        return ReductionResult::NOGOOD;
      case ReductionResult::NEW_ASSIGNMENTS:
        // We'll still continue, just to finish this assignment off.
        result = ReductionResult::NEW_ASSIGNMENTS;
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
