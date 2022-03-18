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
#include "WeightSubgrMono/GraphTheoretic/DerivedGraphs.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

#include <algorithm>

namespace tket {
namespace WeightedSubgraphMonomorphism {

// As soon as we get a new assignment PV->TV, it means that
// each neighbour of PV (in some derived graph) must be either
// unassigned, or assigned to neighbours of TV.
// This merely checks, though, it does not reduce anything.
static bool check_assigned_edges_filter(
      const DerivedGraphStructs::NeighboursAndCounts& pattern_neighbours,
      const DerivedGraphStructs::NeighboursAndCounts& target_neighbours,
      const Assignments& assignments) {
  // Could do fancy back-and-forth iterators
  for(const auto& p_entry : pattern_neighbours) {
    const auto assigned_tv_opt = get_optional_value(assignments, p_entry.first);
    if(!assigned_tv_opt) {
      continue;
    }
    const auto assigned_tv = assigned_tv_opt.value();
    // Now, (pv origin)--pv is an edge in the pattern graph,
    // so (tv origin)--tv MUST be an edge in the target graph.

    const auto target_citer = std::lower_bound(
        target_neighbours.cbegin(),
        target_neighbours.cend(),
        // Also, the edge weight must be at least as great.
        // It's lexicographic ordering, so want (x,y) >= (TV, p-weight).
        // Either x=TV [so that y >= p-weight automatically],
        // or x>TV. (Or "x" doesn't exist).
        std::make_pair(assigned_tv, p_entry.second));
    if(target_citer == target_neighbours.cend() ||
          target_citer->first != assigned_tv) {
      return false;
    }
  }
  return true;
}


static bool reduce_domains_of_neighbours(
      const DerivedGraphStructs::NeighboursAndCounts& pattern_neighbours,
      const DerivedGraphStructs::NeighboursAndCounts& target_neighbours,
      SearchNodeWrapper& node_wrapper,
      Assignments& assignments,
      std::vector<VertexWSM>& reduced_domain) {
  const auto& domains_map = node_wrapper.get().pattern_v_to_possible_target_v;
  // Could do fancy back-and-forth iterators
  for(const auto& p_entry : pattern_neighbours) {
    const auto citer = domains_map.find(p_entry.first);
    if(citer == domains_map.cend()) {
      continue;
    }
    const auto& domain = citer->second;
    reduced_domain.clear();

    // Every PV neighbour (in the derived graph) must map to a TV neighbour.
    // Hence we must INTERSECT the existing domain with this set.
    for(const auto& t_entry : target_neighbours) {
      if(domain.count(t_entry.first) != 0 &&
            t_entry.second >= p_entry.second) {
        reduced_domain.push_back(t_entry.first);
      }
    }
    if(reduced_domain.empty()) {
      return false;
    }
    if(reduced_domain.size() != domain.size()) {
      node_wrapper.overwrite_domain(reduced_domain, p_entry.first, assignments);
    }
  }
  return true;
}


static bool check_and_reduce(const DerivedGraphStructs::NeighboursAndCounts& pattern_neighbours,
      const DerivedGraphStructs::NeighboursAndCounts& target_neighbours,
      SearchNodeWrapper& node_wrapper,
      Assignments& assignments,
      std::vector<VertexWSM>& reduced_domain) {
  return check_assigned_edges_filter(pattern_neighbours, target_neighbours, assignments) &&
      reduce_domains_of_neighbours(pattern_neighbours, target_neighbours,
        node_wrapper, assignments, reduced_domain);
}


bool DerivedGraphsReducer::reduce_domains(
      const FixedData& fixed_data, Assignments& assignments,
      std::size_t number_of_assignments_previously_processed_in_this_node,
      SearchNodeWrapper& node_wrapper,
      DerivedGraphs& derived_pattern_graphs,
      DerivedGraphs& derived_target_graphs) {
  if (fixed_data.target_is_complete) {
    return true;
  }
  
  const auto& chosen_assignments = node_wrapper.get().chosen_assignments;
  const auto& domains_map = node_wrapper.get().pattern_v_to_possible_target_v;

  while (number_of_assignments_previously_processed_in_this_node <
         chosen_assignments.size()) {
    const auto& new_assignment = chosen_assignments
        [number_of_assignments_previously_processed_in_this_node];

    ++number_of_assignments_previously_processed_in_this_node;
    const auto& new_pv = new_assignment.first;
    const auto& new_tv = new_assignment.second;

    if(!check_and_reduce(
            derived_pattern_graphs.d2_graph.get_neighbours(new_pv),
            derived_target_graphs.d2_graph.get_neighbours(new_tv),
            node_wrapper, assignments, m_reduced_domain) 
            
            ||
            
          !check_and_reduce(
            derived_pattern_graphs.d3_graph.get_neighbours(new_pv),
            derived_target_graphs.d3_graph.get_neighbours(new_tv),
            node_wrapper, assignments, m_reduced_domain)
            ) {
      return false;
    }
  }
  return true;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
