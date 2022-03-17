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

#include "WeightSubgrMono/Searching/VariableOrdering.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchNode.hpp"
#include "WeightSubgrMono/Searching/SharedData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

void VariableOrdering::fill_pattern_vertices_with_smallest_domain(
    const SearchNode& node, const Assignments& assignments,
    SharedData& shared_data) {
  m_pattern_vertices_with_smallest_domain.clear();
  std::size_t min_domain_size;
  set_maximum(min_domain_size);
  bool vertices_adjacent_to_assigned = false;

  for (const auto& entry : node.pattern_v_to_possible_target_v) {
    const auto& p_vertex = entry.first;
    TKET_ASSERT(assignments.count(p_vertex) == 0);
    const auto domain_size = entry.second.size();
    TKET_ASSERT(domain_size >= 2);
    if (vertices_adjacent_to_assigned) {
      // We'll only accept if this domain is adjacent,
      // AND has small enough domain.
      if (domain_size > min_domain_size) {
        continue;
      }
      if (!shared_data.fixed_data.pattern_neighbours_data.is_adjacent_to_assigned_pv(p_vertex, assignments)) {
        continue;
      }
      if (domain_size < min_domain_size) {
        // Strictly better.
        m_pattern_vertices_with_smallest_domain.clear();
        min_domain_size = domain_size;
      }
      m_pattern_vertices_with_smallest_domain.push_back(p_vertex);
      continue;
    }
    // We currently have NO vertices adjacent to an assigned one.
    if (shared_data.fixed_data.pattern_neighbours_data.is_adjacent_to_assigned_pv(p_vertex, assignments)) {
      // It's strictly better than anything so far.
      m_pattern_vertices_with_smallest_domain.clear();
      min_domain_size = domain_size;
      m_pattern_vertices_with_smallest_domain.push_back(p_vertex);
      vertices_adjacent_to_assigned = true;
      continue;
    }
    // This is also not adjacent.
    if (domain_size > min_domain_size) {
      continue;
    }
    // Now, we definitely accept.
    if (domain_size < min_domain_size) {
      // Strictly better.
      m_pattern_vertices_with_smallest_domain.clear();
      min_domain_size = domain_size;
    }
    m_pattern_vertices_with_smallest_domain.push_back(p_vertex);
  }
}

VertexWSM VariableOrdering::choose_next_variable(
    const SearchNode& node, const Assignments& assignments,
    SharedData& shared_data) {
  fill_pattern_vertices_with_smallest_domain(node, assignments, shared_data);
  TKET_ASSERT(!m_pattern_vertices_with_smallest_domain.empty());
  return shared_data.rng.get_element(m_pattern_vertices_with_smallest_domain);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
