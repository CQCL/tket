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
#include "WeightSubgrMono/GraphTheoretic/DerivedGraphsContainer.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// We have PV->TV, and PV is at distance 2 in a derived pattern graph
// from some central pattern v.
// Check if TV is at distance <= 2 in the target graph
// fomr the corresponding central target v.
// Return false if not, so the assignment is impossible.
static bool check_assigned_value(VertexWSM pv,
        const std::vector<VertexWSM>& target_distance_two_vertices,
        const DerivedGraphsCalculator::NeighboursAndCounts& target_distance_one_vertices,
        const Assignments& assignments) {

  const VertexWSM tv = assignments.at(pv);

  const auto expected_tv_citer = std::lower_bound(
      target_distance_one_vertices.cbegin(),
      target_distance_one_vertices.cend(),
      std::make_pair(tv, std::size_t(1)));

  const bool tv_is_at_distance_one =
      expected_tv_citer != target_distance_one_vertices.cend() &&
      expected_tv_citer->first == tv;

  // We want to return TRUE if and only if this TV is at distance 1 or 2.
  return tv_is_at_distance_one ||
        std::binary_search(
              target_distance_two_vertices.cbegin(),
              target_distance_two_vertices.cend(), tv);
}


// Any pv2 at distance two (in a derived graph) from some central
// pv1 must map to a tv2 at distance <= 2 from tv1 in the target graph.
// This does the actual domain reduction.
static bool reduce_domains_with_derived_distance_two_vertices(
        const std::vector<VertexWSM>& pattern_distance_two_vertices,
        const std::vector<VertexWSM>& target_distance_two_vertices,
        const DerivedGraphsCalculator::NeighboursAndCounts& target_distance_one_vertices,
        Assignments& assignments,
        SearchNodeWrapper& node_wrapper,
        std::vector<VertexWSM>& reduced_domain_work_vector) {

  if(pattern_distance_two_vertices.size() > target_distance_two_vertices.size() +
          target_distance_one_vertices.size()) {
    return false;
  }
  const auto& domains_map = node_wrapper.get().pattern_v_to_possible_target_v;          
  for(auto pv : pattern_distance_two_vertices) {
    const auto& domain_citer = domains_map.find(pv);
    // PV is at distance 2 from the root, so the corresponding TV
    // is at distance 1 or 2 in the derived target graph.
    if (domain_citer == domains_map.cend()) {
      // PV has no domain, so must be assigned.
      if(!check_assigned_value(pv, target_distance_two_vertices,
              target_distance_one_vertices, assignments)) {
        return false;
      }
      continue;
    }

    // This PV (at distance 2 from the root pv in the derived graph)
    // is not yet assigned.
    const auto& current_domain = domain_citer->second;
    reduced_domain_work_vector.clear();
    for(auto tv : target_distance_two_vertices) {
      if(current_domain.count(tv) != 0) {
        reduced_domain_work_vector.push_back(tv);
      }
    }
    // We currently have all distance 2 TV in reduced_domain_work_vector,
    // we must add the distance 1 vertices also.
    // These are necessarily distinct.
    for(const auto& entry : target_distance_one_vertices) {
      const auto& tv = entry.first;
      if(current_domain.count(tv) != 0) {
        reduced_domain_work_vector.push_back(tv);
      }
    }
    // Now, reduce the domain.
    if(reduced_domain_work_vector.empty()) {
      return false;
    }
    if(reduced_domain_work_vector.size() != current_domain.size()) {
      // We don't need reduced_domain_work_vector to be sorted.
      node_wrapper.overwrite_domain(reduced_domain_work_vector, pv, assignments);
    }
  }
  return true;
}


// "func" knows how to find neighbours in the derived graph,
// returning a permanent reference.
template<class GetDerivedNeighboursFunc>
static const std::vector<VertexWSM>& get_derived_distance_two_vertices(
      VertexWSM v,
      std::map<VertexWSM, std::vector<VertexWSM>>& map_to_fill,
      const GetDerivedNeighboursFunc& func,
      std::set<VertexWSM>& work_set) {
  {
    const auto iter = map_to_fill.find(v);
    if(iter != map_to_fill.end()) {
      return iter->second;
    }
  }
  const DerivedGraphsCalculator::NeighboursAndCounts& derived_neighbours = func(v);
  work_set.clear();
  for(const auto& entry : derived_neighbours) {
    const auto& v1 = entry.first;
    for(const auto& inner_entry : func(v1)) {
      // v2 is at distance <= 2 from v.
      const auto& v2 = inner_entry.first;
      if(v2 == v) {
        continue;
      }
      // Ensure that v2 is NOT at distance 1 from v.
      const auto inner_citer = std::lower_bound(
          derived_neighbours.cbegin(),
          derived_neighbours.cend(),
          std::make_pair(v2, std::size_t(1)));
      if(inner_citer == derived_neighbours.cend() ||
            inner_citer->first != v2) {
        // NOT at distance 1.
        work_set.insert(v2);
      }
    }
  }
  auto& vertices = map_to_fill[v];
  vertices = { work_set.cbegin(), work_set.cend() };
  return vertices;
}


template<class GetNeighboursFromVertexDataFunc>
static bool reduce_with_derived_distance_two_vertices(
        const FixedData& fixed_data, Assignments& assignments,
        SearchNodeWrapper& node_wrapper,
        DerivedGraphsContainer& derived_graphs,
        VertexWSM pv, VertexWSM tv,
        std::map<VertexWSM, std::vector<VertexWSM>>& pattern_graph_map,
        std::map<VertexWSM, std::vector<VertexWSM>>& target_graph_map,
        std::set<VertexWSM>& work_set,
        std::vector<VertexWSM>& reduced_domain_work_vector,
        const GetNeighboursFromVertexDataFunc& func) {
        
  const auto& pv_at_distance_two = get_derived_distance_two_vertices(
      pv, pattern_graph_map,
      [&fixed_data, &derived_graphs, &func](VertexWSM v) -> const DerivedGraphsCalculator::NeighboursAndCounts& {
        return func(derived_graphs.get_pattern_v_data_permanent_reference(
            v, fixed_data.pattern_neighbours_data));
      },
      work_set);

  if(pv_at_distance_two.empty()) {
    return true;
  }

  const auto get_target_neighbours_func = [&fixed_data, &derived_graphs, &func](VertexWSM v) -> const DerivedGraphsCalculator::NeighboursAndCounts& {
    return func(derived_graphs.get_target_v_data_permanent_reference(
        v, fixed_data.target_neighbours_data));
  };

  const auto& tv_at_distance_two = get_derived_distance_two_vertices(
      tv, target_graph_map,
      get_target_neighbours_func, work_set);

  return reduce_domains_with_derived_distance_two_vertices(
        pv_at_distance_two, tv_at_distance_two,
        get_target_neighbours_func(tv),
        assignments, node_wrapper, reduced_domain_work_vector);
}


bool DerivedGraphsReducer::reduce_domains(
      const FixedData& fixed_data, Assignments& assignments,
      std::size_t number_of_assignments_previously_processed_in_this_node,
      SearchNodeWrapper& node_wrapper,
      DerivedGraphsContainer& derived_graphs) {
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
    const auto& pv_data = derived_graphs.get_pattern_v_data_permanent_reference(
        new_pv, fixed_data.pattern_neighbours_data);
    const auto& tv_data = derived_graphs.get_target_v_data_permanent_reference(
        new_tv, fixed_data.target_neighbours_data);

    if(!reduce_with_derived_weighted_neighbours(pv_data.depth_2_neighbours,
            tv_data.depth_2_neighbours, assignments, node_wrapper) ||
          !reduce_with_derived_weighted_neighbours(pv_data.depth_3_neighbours,
            tv_data.depth_3_neighbours, assignments, node_wrapper)) {
      return false;
    }
    // Try more distance derived neighbours.
    if(!reduce_with_derived_distance_two_vertices(fixed_data, assignments,
          node_wrapper, derived_graphs, new_pv, new_tv,
          m_pattern_neighbours_at_distance_two_in_d2,
          m_target_neighbours_at_distance_two_in_d2,
          m_work_set,
          m_reduced_domain,
          [](const DerivedGraphsContainer::VertexData& vdata) -> const DerivedGraphsCalculator::NeighboursAndCounts& {
            return vdata.depth_2_neighbours;
          })) {
      return false;
    }
    if(!reduce_with_derived_distance_two_vertices(fixed_data, assignments,
          node_wrapper, derived_graphs, new_pv, new_tv,
          m_pattern_neighbours_at_distance_two_in_d3,
          m_target_neighbours_at_distance_two_in_d3,
          m_work_set,
          m_reduced_domain,
          [](const DerivedGraphsContainer::VertexData& vdata) -> const DerivedGraphsCalculator::NeighboursAndCounts& {
            return vdata.depth_3_neighbours;
          })) {
      return false;
    }
  }
  return true;
}


bool DerivedGraphsReducer::reduce_with_derived_weighted_neighbours(
          const DerivedGraphsCalculator::NeighboursAndCounts& pattern_neighbours,
          const DerivedGraphsCalculator::NeighboursAndCounts& target_neighbours,
          Assignments& assignments,
          SearchNodeWrapper& node_wrapper) {
  const auto& domains_map = node_wrapper.get().pattern_v_to_possible_target_v;
  for(const auto& pv_entry : pattern_neighbours) {
    const auto& pv_other = pv_entry.first;
    const auto& p_edge_weight = pv_entry.second;
    const auto& domain_citer = domains_map.find(pv_other);
    if (domain_citer == domains_map.cend()) {
      // If it has no domain, it must be assigned.
      const VertexWSM tv_other = assignments.at(pv_other);
      const auto expected_tv_citer = std::lower_bound(
          target_neighbours.cbegin(), target_neighbours.cend(),
          std::make_pair(tv_other, p_edge_weight));
      if(expected_tv_citer == target_neighbours.cend() ||
            expected_tv_citer->first != tv_other) {
        return false;
      }
      continue;
    }
    const auto& current_domain = domain_citer->second;
    m_reduced_domain.clear();
    for(const auto& entry : target_neighbours) {
      if(entry.second >= p_edge_weight &&
            current_domain.count(entry.first) != 0) {
        m_reduced_domain.push_back(entry.first);
      }
    }
    if(m_reduced_domain.empty()) {
      return false;
    }
    if(m_reduced_domain.size() != current_domain.size()) {
      node_wrapper.overwrite_domain(m_reduced_domain, pv_other, assignments);
    }
  }
  return true;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
