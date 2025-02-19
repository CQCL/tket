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

#include "tkwsm/Searching/VariableOrdering.hpp"

#include <algorithm>
#include <tkrng/RNG.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

VariableOrdering::Result VariableOrdering::get_variable(
    DomainsAccessor& accessor, RNG& rng) {
  m_pv_list.clear();
  std::size_t min_domain_size;
  set_maximum(min_domain_size);
  Result result;
  m_work_vector.clear();

  // NOTE: it might seem natural to give priority to PV adjacent to
  // newly assigned vertices, but tests showed that this is actually
  // a bad idea; it's faster and simpler just to consider all
  // unassigned vertices equally.
  for (VertexWSM pv : accessor.get_unassigned_pattern_vertices_superset()) {
    const auto domain_size = accessor.get_domain_size(pv);
    if (domain_size == 0) {
      result.empty_domain = true;
      return result;
    }
    if (domain_size == 1) {
      continue;
    }
    // The point is, we're guaranteed to consider every
    // unassigned PV (and maybe some assigned ones).
    // When we move down the search tree, the set of unassigned
    // vertices decreases; thus we can save time by only searching
    // this subset of PV.
    // When we move up and return to a node,
    // as we reduce further some vertices become assigned, BUT that
    // doesn't bother us; it just means we have a few extra PV
    // to search through next time.
    // This is quicker than trying to maintain a strictly accurate
    // list of unassigned PV (which would need a std::set).
    m_work_vector.push_back(pv);

    if (domain_size > min_domain_size) {
      continue;
    }
    if (domain_size < min_domain_size) {
      // Strictly better.
      m_pv_list.clear();
      min_domain_size = domain_size;
    }
    m_pv_list.emplace_back(pv);
  }

  // NOTE: what we're doing is, building up the new vector
  // (the pattern_vertices_superset) in m_work_vector.
  // Then when we reach this point, we just write the data back
  // into pattern_vertices_superset (by a vector swap, which is O(1) time).
  // This is practically just as efficient as overwriting
  // pattern_vertices_superset directly, but simpler because
  // we need to copy and manipulate the old data to produce the new data,
  // so we don't have to do it all within the SAME vector object
  // (although we could do).
  accessor.get_unassigned_pattern_vertices_superset_to_overwrite().swap(
      m_work_vector);

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
