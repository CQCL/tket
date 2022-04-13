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

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/GraphTheoretic/NearNeighboursData.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "WeightSubgrMono/Searching/SearchBranch.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DistancesReducer::DistancesReducer(
    const NeighboursData& pattern_ndata, NearNeighboursData& pattern_near_ndata,
    const NeighboursData& target_ndata, NearNeighboursData& target_near_ndata)
    : m_raw_reducer(
          pattern_ndata, pattern_near_ndata, target_ndata, target_near_ndata) {}

bool DistancesReducer::check(
    const std::pair<VertexWSM, VertexWSM>& assignment,
    unsigned distance) const {
  return m_raw_reducer.check(assignment, distance);
}

void DistancesReducer::reset(unsigned distance_value) {
  m_data.resize(distance_value);
  std::fill(m_data.begin(), m_data.end(), 0);
}

DistancesReducer::Result DistancesReducer::operator()(NodeWSM& node) {
  // May change, but only indirectly via "node".
  const auto& new_assignments = node.get_new_assignments();
  const auto size = new_assignments.size();
  Result result;

  for (unsigned ii = 0; ii < m_data.size(); ++ii) {
    auto& index = m_data[ii];
    while (index < size) {
      const auto reduce_result =
          m_raw_reducer(new_assignments[index], node, ii + 1);
      ++index;
      switch (reduce_result) {
        case DistancesRawReducer::Result::IMPOSSIBLE_ASSIGNMENT:
          result.impossible_assignment = new_assignments[index];
          result.new_assignments_created = false;
          result.nogood_found = true;
          return result;

        case DistancesRawReducer::Result::IMPOSSIBLE_NODE:
          result.new_assignments_created = false;
          result.nogood_found = true;
          return result;

        case DistancesRawReducer::Result::SUCCESS:
          break;
      }
      if (size != new_assignments.size()) {
        result.new_assignments_created = true;
        result.nogood_found = false;
        return result;
      }
    }
  }
  result.new_assignments_created = false;
  result.nogood_found = false;
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
