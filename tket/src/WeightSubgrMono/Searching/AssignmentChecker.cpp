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

#include "WeightSubgrMono/Searching/AssignmentChecker.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Reducing/DerivedGraphsReducer.hpp"
#include "WeightSubgrMono/Reducing/DistancesReducer.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

AssignmentChecker::AssignmentChecker(
    DerivedGraphsReducer& derived_graphs_reducer,
    DistancesReducer& distances_reducer)
    : m_derived_graphs_reducer(derived_graphs_reducer),
      m_distances_reducer(distances_reducer) {}

bool AssignmentChecker::operator()(
    const std::pair<VertexWSM, VertexWSM>& assignment,
    unsigned distances_reducer_max_dist) {
  auto& acceptable_tv_set = m_checked_assignments[assignment.first];
  if (acceptable_tv_set.count(assignment.second) != 0) {
    return true;
  }
  // We haven't checked before, we must do so now.
  const bool acceptable =
      m_derived_graphs_reducer.check(assignment) &&
      m_distances_reducer.check(assignment, distances_reducer_max_dist);
  if (acceptable) {
    acceptable_tv_set.insert(assignment.second);
  }
  return acceptable;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
