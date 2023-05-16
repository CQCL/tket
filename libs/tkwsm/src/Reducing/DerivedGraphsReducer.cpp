// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tkwsm/Reducing/DerivedGraphsReducer.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/GraphTheoretic/FilterUtils.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphsReducer::DerivedGraphsReducer(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata)
    : m_derived_pattern_graphs(pattern_ndata, m_calculator),
      m_derived_target_graphs(target_ndata, m_calculator) {}

bool DerivedGraphsReducer::check(std::pair<VertexWSM, VertexWSM> assignment) {
  const DerivedGraphs::VertexData pattern_vdata =
      m_derived_pattern_graphs.get_data(assignment.first);
  const DerivedGraphs::VertexData target_vdata =
      m_derived_target_graphs.get_data(assignment.second);

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

ReductionResult DerivedGraphsReducer::reduce_with_derived_data(
    const DerivedGraphStructs::NeighboursAndCounts&
        pattern_derived_neighbours_data,
    const DerivedGraphStructs::NeighboursAndCounts&
        target_derived_neighbours_data,
    VertexWSM root_pattern_vertex, DomainsAccessor& accessor,
    boost::dynamic_bitset<>& work_bitset) {
  // If we do create a new assignment, we will continue
  // so that this assignment at least is fully processed.
  // This is simpler than trying to split the work up
  // and resume halfway through (although, really, that's what we should do
  // to squeeze the maximum performance from lazy evaluation).
  bool found_new_assignment = false;

  for (const std::pair<VertexWSM, DerivedGraphStructs::Count>& p_entry :
       pattern_derived_neighbours_data) {
    // This PV is a neighbour of the root PV, in some derived graph.
    const VertexWSM& pv = p_entry.first;
    {
      const auto& domain = accessor.get_domain(pv);
      if (other_vertex_reduction_can_be_skipped_by_symmetry(
              domain, accessor, root_pattern_vertex, pv)) {
        continue;
      }
      work_bitset.resize(domain.size());
      work_bitset.reset();
    }
    // This is the edge weight of PV--(root PV) in the derived graph.
    const DerivedGraphStructs::Count& p_count = p_entry.second;

    for (const auto& entry : target_derived_neighbours_data) {
      if (entry.second >= p_count) {
        TKET_ASSERT(!work_bitset.test_set(entry.first));
      }
    }

    switch (
        accessor.intersect_domain_with_swap(pv, work_bitset).reduction_result) {
      case ReductionResult::NOGOOD:
        return ReductionResult::NOGOOD;
      case ReductionResult::NEW_ASSIGNMENTS:
        found_new_assignment = true;
        break;
      case ReductionResult::SUCCESS:
        break;
    }
  }
  if (found_new_assignment) {
    return ReductionResult::NEW_ASSIGNMENTS;
  }
  return ReductionResult::SUCCESS;
}

ReductionResult DerivedGraphsReducer::reduce(
    std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
    boost::dynamic_bitset<>& work_set) {
  const auto pattern_vdata =
      m_derived_pattern_graphs.get_data(assignment.first);
  const auto target_vdata = m_derived_target_graphs.get_data(assignment.second);

  // All the d2, d3 pattern neighbours must map to target neighbours.
  // Furthermore, the edge weights must increase in the target graph.
  // Even if a new assignment is found, let's reduce with all graphs
  // before returning. Thus, we know at least that the assignment is fully
  // processed.
  const auto d2_result = reduce_with_derived_data(
      *pattern_vdata.d2_neighbours, *target_vdata.d2_neighbours,
      assignment.first, accessor, work_set);

  if (d2_result == ReductionResult::NOGOOD) {
    return ReductionResult::NOGOOD;
  }
  const auto d3_result = reduce_with_derived_data(
      *pattern_vdata.d3_neighbours, *target_vdata.d3_neighbours,
      assignment.first, accessor, work_set);

  if (d3_result == ReductionResult::NOGOOD) {
    return ReductionResult::NOGOOD;
  }
  if (d2_result == ReductionResult::NEW_ASSIGNMENTS ||
      d3_result == ReductionResult::NEW_ASSIGNMENTS) {
    return ReductionResult::NEW_ASSIGNMENTS;
  }
  return ReductionResult::SUCCESS;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
