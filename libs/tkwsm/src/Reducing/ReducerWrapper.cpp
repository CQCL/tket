// Copyright Quantinuum
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

#include "tkwsm/Reducing/ReducerWrapper.hpp"

#include "tkwsm/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool ReducerInterface::other_vertex_reduction_can_be_skipped_by_symmetry(
    const boost::dynamic_bitset<>& other_domain,
    const DomainsAccessor& accessor, VertexWSM this_vertex,
    VertexWSM other_vertex) {
  // If this other pv was already assigned in a previous node
  // (i.e., its domain was the same as now), then the reducer
  // already reduced this domain when that node was reduced
  // (and, we could only have reached this current node by moving down,
  // so our current domain is a subset of that one).
  //
  // Otherwise, if pv1, pv2 both had domains reduced to size 1 in
  // the current node, we only need to perform one reduction.
  // We could decide by knowing which vertex was assigned earlier,
  // but that information is not available (although it could be deduced
  // using some labour).
  // Instead, we decide by using vertex numbers.
  return other_domain.count() == 1 &&
         (!accessor.domain_created_in_current_node(other_vertex) ||
          other_vertex < this_vertex);
}

ReducerWrapper::ReducerWrapper(ReducerInterface& reducer_interface)
    : m_reducer(reducer_interface) {}

void ReducerWrapper::clear() { m_number_of_processed_assignments = 0; }

bool ReducerWrapper::check(std::pair<VertexWSM, VertexWSM> assignment) {
  return m_reducer.check(assignment);
}

ReductionResult ReducerWrapper::reduce(
    DomainsAccessor& accessor, boost::dynamic_bitset<>& work_set) {
  auto result = ReductionResult::SUCCESS;
  for (const std::vector<std::pair<VertexWSM, VertexWSM>>& new_assignments =
           accessor.get_new_assignments();
       m_number_of_processed_assignments < new_assignments.size();) {
    result = m_reducer.reduce(
        new_assignments[m_number_of_processed_assignments], accessor, work_set);

    // We should increment this BEFORE breaking out of the loop.
    // We've finished processing this particular assignment,
    // so don't want to waste time reprocessing it
    // next time "reduce" is called.
    ++m_number_of_processed_assignments;
    if (result != ReductionResult::SUCCESS) {
      break;
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
