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

#include "WeightSubgrMono/Reducing/DerivedGraphsReducer.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "WeightSubgrMono/Searching/NodeWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphsReducer::DerivedGraphsReducer(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata)
    : m_derived_pattern_graphs(
          pattern_ndata, m_calculator, m_storage, m_counts_storage),
      m_derived_target_graphs(
          target_ndata, m_calculator, m_storage, m_counts_storage),
      m_number_of_assignments_processed(0) {}

bool DerivedGraphsReducer::check(
    const std::pair<VertexWSM, VertexWSM>& assignment) {
  const auto pattern_vdata =
      m_derived_pattern_graphs.get_data(assignment.first);
  const auto target_vdata = m_derived_target_graphs.get_data(assignment.second);

  return pattern_vdata.triangle_count <= target_vdata.triangle_count &&
         pattern_vdata.d2_neighbours->size() <=
             target_vdata.d2_neighbours->size() &&
         pattern_vdata.d3_neighbours->size() <=
             target_vdata.d3_neighbours->size() &&

         // It happens that the sorted degree sequences filter algorithm
         // is identical to the algorithm here (the pattern edge weight must be
         // <= the target edge weight).
         FilterUtils::compatible_sorted_degree_sequences(
             *pattern_vdata.d2_sorted_counts_iter,
             *target_vdata.d2_sorted_counts_iter) &&
         FilterUtils::compatible_sorted_degree_sequences(
             *pattern_vdata.d3_sorted_counts_iter,
             *target_vdata.d3_sorted_counts_iter);
}


DerivedGraphsReducer::ReductionResult DerivedGraphsReducer::reduce(
    const DerivedGraphStructs::NeighboursAndCounts& pattern_neighbours,
    const DerivedGraphStructs::NeighboursAndCounts& target_neighbours,
    NodeWSM& node) {
  // If we do create a new assignment, we will continue
  // so that this assignment at least is fully processed.
  bool found_new_assignment = false;

  // We will not change this directly;
  // but it may change indirectly as we reduce the node.
  // The reference is still valid, however.
  const auto& domains = node.get_possible_assignments();

  for (const auto& p_entry : pattern_neighbours) {
    const VertexWSM& pv = p_entry.first;
    // This is the edge weight in the derived graph,
    // from the (unknown) root PV to this new PV.
    const auto& p_count = p_entry.second;
    const auto& domain = domains.at(pv);

    fill_intersection(
        domain, target_neighbours, m_new_domain,
        [](const std::pair<VertexWSM, DerivedGraphStructs::Count>& pair) {
          return pair.first;
        },

        // To get the next pair (tv2, count)  after tv1,
        // we add the extra condition that count >= p-count.
        // (Since, edges in the derived p-graph must map to edges with equal
        // or greater counts in the target graph).
        // Luckily, this fits in exactly with the general framework.
        [p_count](VertexWSM tv) { return std::make_pair(tv, p_count); });

    if (m_new_domain.size() == domain.size()) {
      continue;
    }
    if (m_new_domain.empty()) {
      return ReductionResult::NOGOOD;
    }
    if (m_new_domain.size() == 1) {
      found_new_assignment = true;
    }
    node.overwrite_domain(pv, m_new_domain);
  }
  if (found_new_assignment) {
    return ReductionResult::NEW_ASSIGNMENT;
  }
  return ReductionResult::SUCCESS;
}

DerivedGraphsReducer::ReductionResult DerivedGraphsReducer::reduce(
    const std::pair<VertexWSM, VertexWSM>& assignment, NodeWSM& node) {
  const auto pattern_vdata =
      m_derived_pattern_graphs.get_data(assignment.first);
  const auto target_vdata = m_derived_target_graphs.get_data(assignment.second);

  // All the d2, d3 pattern neighbours must map to target neighbours.
  // Furthermore, the edge weights must increase in the target graph.
  // Even if a new assignment is found, let's reduce with all graphs
  // before returning. Thus, we know at least that the assignment is fully
  // processed.
  const auto d2_result =
      reduce(*pattern_vdata.d2_neighbours, *target_vdata.d2_neighbours, node);

  if (d2_result == ReductionResult::NOGOOD) {
    return ReductionResult::NOGOOD;
  }
  const auto d3_result =
      reduce(*pattern_vdata.d3_neighbours, *target_vdata.d3_neighbours, node);

  if (d3_result == ReductionResult::NOGOOD) {
    return ReductionResult::NOGOOD;
  }
  if (d2_result == ReductionResult::SUCCESS &&
      d3_result == ReductionResult::SUCCESS) {
    return ReductionResult::SUCCESS;
  }
  return ReductionResult::NEW_ASSIGNMENT;
}

void DerivedGraphsReducer::clear() { m_number_of_assignments_processed = 0; }

DerivedGraphsReducer::ReductionResult DerivedGraphsReducer::reduce(
    NodeWSM& node) {
  // The assignments are stored in a fixed place in "node",
  // separately from the rest of the node data,
  // so are not invalidated.
  const auto& assignments = node.get_new_assignments();
  while (m_number_of_assignments_processed < assignments.size()) {
    const auto result =
        reduce(assignments[m_number_of_assignments_processed], node);
    ++m_number_of_assignments_processed;
    if (result != ReductionResult::SUCCESS) {
      return result;
    }
  }
  return ReductionResult::SUCCESS;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
