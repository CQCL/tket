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

#include "tkwsm/Reducing/NeighboursReducer.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Common/SetIntersection.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

//#include "WeightSubgrMono/Searching/NodesRawData.hpp"

#include "WeightSubgrMono/Common/TemporaryRefactorCode.hpp"

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

  TemporaryRefactorCode();
  boost::dynamic_bitset<> work_bitset;


  for (const std::pair<VertexWSM, WeightWSM>& p_neighbour_entry :
       m_pattern_ndata.get_neighbours_and_weights(assignment.first)) {
    const VertexWSM& p_neighbour = p_neighbour_entry.first;

    const boost::dynamic_bitset<>& bitset_domain = accessor.get_domain(p_neighbour);
    if (other_vertex_reduction_can_be_skipped_by_symmetry(
            bitset_domain, accessor, assignment.first, p_neighbour)) {
      continue;
    }

    work_bitset.resize(bitset_domain.size());
    work_bitset.reset();
    for(const auto& tv_weight_pair : target_neighbours_and_weights) {
      TKET_ASSERT(!work_bitset.test_set(tv_weight_pair.first));
    }

    switch (accessor.intersect_domain_with_swap(p_neighbour, work_bitset).reduction_result) {
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
