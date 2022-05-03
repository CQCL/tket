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

#include "Utils/RNG.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool VariableOrdering::check_candidate(
    VertexWSM pv, std::size_t domain_size, std::size_t& min_domain_size,
    std::set<VertexWSM>& current_node_unassigned_vertices_to_overwrite,
    bool write_unassigned_vertices) {
  if (domain_size == 0) {
    return false;
  }
  if (write_unassigned_vertices && domain_size > 1) {
    current_node_unassigned_vertices_to_overwrite.insert(pv);
  }
  if (domain_size == 1 || domain_size > min_domain_size) {
    return true;
  }
  if (domain_size < min_domain_size) {
    // Strictly better.
    m_pv_list.clear();
    min_domain_size = domain_size;
  }
  m_pv_list.emplace_back(pv);
  return true;
}

VariableOrdering::Result VariableOrdering::get_variable(
    const DomainsAccessor& accessor, RNG& rng,
    std::set<VertexWSM>& current_node_unassigned_vertices_to_overwrite) {
  m_pv_list.clear();
  std::size_t min_domain_size;
  set_maximum(min_domain_size);
  Result result;

  for (VertexWSM pv : accessor.get_candidate_vertices_for_assignment()) {
    const auto domain_size = accessor.get_domain(pv).size();
    if (!check_candidate(
            pv, domain_size, min_domain_size,
            // The candidate vertices are NOT necessarily
            // all the unassigned vertices, therefore we DON'T
            // try to fill the set.
            current_node_unassigned_vertices_to_overwrite, false)) {
      result.empty_domain = true;
      return result;
    }
  }
  if (m_pv_list.empty()) {
    // No candidates are unassigned, so look through ALL vertices.
    const bool write_unassigned_vertices =
        current_node_unassigned_vertices_to_overwrite.empty();
    for (VertexWSM pv : accessor.get_unassigned_pattern_vertices_superset()) {
      const auto domain_size = accessor.get_domain(pv).size();
      if (!check_candidate(
              pv, domain_size, min_domain_size,
              current_node_unassigned_vertices_to_overwrite,
              write_unassigned_vertices)) {
        result.empty_domain = true;
        return result;
      }
    }
  }
  result.empty_domain = false;
  if (m_pv_list.empty()) {
    return result;
  }
  // We have a vertex!
  const auto index = rng.get_size_t(m_pv_list.size() - 1);
  result.variable_opt = m_pv_list[index];
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
