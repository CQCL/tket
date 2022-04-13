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

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool VariableOrdering::check_candidate(
    VertexWSM pv, std::size_t domain_size, std::size_t& min_domain_size) {
  if (domain_size == 0) {
    return false;
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
    const PossibleAssignments& possible_assignments,
    const std::set<VertexWSM>& candidate_vertices, RNG& rng) {
  m_pv_list.clear();
  std::size_t min_domain_size;
  set_maximum(min_domain_size);
  Result result;

  for (VertexWSM pv : candidate_vertices) {
    const auto domain_size = possible_assignments.at(pv).size();
    if (!check_candidate(pv, domain_size, min_domain_size)) {
      result.empty_domain = true;
      return result;
    }
  }
  if (m_pv_list.empty()) {
    // No candidates are unassigned, so look through ALL vertices.
    for (const auto& entry : possible_assignments) {
      const VertexWSM& pv = entry.first;
      const auto domain_size = entry.second.size();
      if (!check_candidate(pv, domain_size, min_domain_size)) {
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
