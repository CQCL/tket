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

#include "WeightSubgrMono/Reducing/CloseVerticesFilter.hpp"

#include <cmath>
#include <sstream>

#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// In experiments, this appeared to be the best value
// (in combination with all other reductions).
CloseVerticesFilter::CloseVerticesFilter() : m_num_levels(2) {}

bool CloseVerticesFilter::reduce(
    const FixedData& fixed_data, Assignments& assignments,
    std::size_t number_of_assignments_previously_processed_in_this_node,
    SearchNodeWrapper& node_wrapper) {
  if (fixed_data.target_is_complete) {
    return true;
  }
  const auto& chosen_assignments = node_wrapper.get().chosen_assignments;

  while (number_of_assignments_previously_processed_in_this_node <
         chosen_assignments.size()) {
    const auto& new_assignment = chosen_assignments
        [number_of_assignments_previously_processed_in_this_node];

    ++number_of_assignments_previously_processed_in_this_node;
    const auto& new_pv = new_assignment.first;
    const auto& new_tv = new_assignment.second;

    // Get the existing data for PV, TV,
    // calculating it if this is the first usage.
    auto& p_data = m_pattern_data[new_pv];
    if (p_data.empty()) {
      // There are no isolated vertices,
      // so empty data MUST mean uninitialised.
      m_close_neighbours_calculator(
          m_num_levels, p_data, fixed_data.pattern_neighbours_data, new_pv);
    }
    auto& t_data = m_target_data[new_tv];
    if (t_data.empty()) {
      m_close_neighbours_calculator(
          m_num_levels, t_data, fixed_data.target_neighbours_data, new_tv);
    }

    // NOTE: m_pattern_data, m_target_data are separate objects,
    // and we only ask for a single element from each.
    // Therefore the p_data, t_data references ARE valid
    // within this iteration of the loop,
    // even if map reallocation took place.
    const auto& domains_map = node_wrapper.get().pattern_v_to_possible_target_v;

    for (unsigned ii = 0; ii < p_data.size(); ++ii) {
      // Every p-vertex at this distance must correspond to a t-vertex
      // at an equal or smaller distance in the target graph
      // (since subgraph isomorphisms decrease distances).
      for (VertexWSM p_vertex : p_data[ii]) {
        const auto citer = domains_map.find(p_vertex);
        if (citer == domains_map.cend()) {
          // This pv is NOT unassigned, so it should already
          // have been assigned elsewhere...
          // does this contradict its assignment?
          const VertexWSM tv = assignments.at(p_vertex);
          if (t_vertex_is_suitable(t_data, tv, ii)) {
            continue;
          }
          return false;
        }

        // Now, we must INTERSECT this domain with the UNION
        // of all the t-sets at level <= ii.
        // Is there a really fast way to do this?
        m_reduced_domain.clear();
        const auto& domain = citer->second;

        for (VertexWSM tv : domain) {
          if (t_vertex_is_suitable(t_data, tv, ii)) {
            // Note: automatically sorted.
            m_reduced_domain.push_back(tv);
          }
        }
        if (m_reduced_domain.empty()) {
          return false;
        }
        node_wrapper.overwrite_domain(m_reduced_domain, p_vertex, assignments);
      }
    }
  }
  return true;
}

bool CloseVerticesFilter::t_vertex_is_suitable(
    const CloseNeighboursCalculator::CloseVertices& close_target_vertices,
    VertexWSM original_target_vertex, unsigned current_p_level) {
  for (unsigned jj = 0;
       jj <= current_p_level && jj < close_target_vertices.size(); ++jj) {
    if (std::binary_search(
            close_target_vertices[jj].cbegin(),
            close_target_vertices[jj].cend(), original_target_vertex)) {
      return true;
    }
  }
  return false;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
