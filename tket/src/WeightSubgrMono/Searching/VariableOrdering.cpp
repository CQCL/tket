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

VariableOrdering::Result VariableOrdering::get_variable(
    const DomainsAccessor& accessor, RNG& rng) {
  m_pv_list.clear();
  std::size_t min_domain_size;
  set_maximum(min_domain_size);
  Result result;

  // NOTE: it might seem natural to give priority to PV adjacent to
  // newly assigned vertices, but tests showed that this is actually
  // a bad idea; it's faster and simpler just to consider all
  // unassigned vertices equally.
  for (VertexWSM pv : accessor.get_unassigned_pattern_vertices_superset()) {
    const auto domain_size = accessor.get_domain(pv).size();
    if (domain_size == 0) {
      result.empty_domain = true;
      return result;
    }
    if (domain_size == 1 || domain_size > min_domain_size) {
      continue;
    }
    if (domain_size < min_domain_size) {
      // Strictly better.
      m_pv_list.clear();
      min_domain_size = domain_size;
    }
    m_pv_list.emplace_back(pv);
  }
  result.empty_domain = false;
  if (m_pv_list.empty()) {
    return result;
  }
  // We have an unassigned vertex!
  const auto index = rng.get_size_t(m_pv_list.size() - 1);
  result.variable_opt = m_pv_list[index];
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
