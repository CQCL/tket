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

#include "TokenSwapping/PartialMappingLookup.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

const ExactMappingLookup::Result& PartialMappingLookup::operator()(
    const VertexMapping& desired_mapping, const vector<Swap>& edges,
    const std::set<size_t>& vertices_with_tokens_at_start,
    unsigned max_number_of_swaps) {
  const auto& exact_mapping_result =
      m_exact_mapping_lookup(desired_mapping, edges, max_number_of_swaps);

  if (exact_mapping_result.success && exact_mapping_result.swaps.empty()) {
    return exact_mapping_result;
  }

  // Are there any empty vertices?
  m_empty_source_vertices.clear();
  m_empty_target_vertices.clear();
  for (const auto& entry : desired_mapping) {
    if (vertices_with_tokens_at_start.count(entry.first) == 0) {
      m_empty_source_vertices.push_back(entry.first);
      m_empty_target_vertices.push_back(entry.second);
    }
  }
  if (m_empty_source_vertices.size() <= 1) {
    // Only an exact lookup is needed (or possible).
    return exact_mapping_result;
  }

  // There are some empty vertices at the start.
  // These END UP at empty target vertices
  // (which, of course, might be completely different!)
  // For next_permutation, let's permute the empty SOURCE vertices.
  // They are already sorted, thus already at the first permutation
  // in the ordering, because they came from the keys of desired_mapping.
  {
    const bool next_permutation = std::next_permutation(
        m_empty_source_vertices.begin(), m_empty_source_vertices.end());
    TKET_ASSERT(next_permutation);
  }
  m_altered_mapping = desired_mapping;

  for (unsigned perm_count = 0;;) {
    for (unsigned ii = 0; ii < m_empty_source_vertices.size(); ++ii) {
      m_altered_mapping[m_empty_source_vertices[ii]] =
          m_empty_target_vertices[ii];
    }
    const auto& exact_map_result_for_permuted_vertices =
        m_exact_mapping_lookup.improve_upon_existing_result(
            m_altered_mapping, edges, max_number_of_swaps);

    if (exact_map_result_for_permuted_vertices.success &&
        exact_map_result_for_permuted_vertices.swaps.empty()) {
      return exact_map_result_for_permuted_vertices;
    }
    ++perm_count;
    if (perm_count >= m_parameters.max_number_of_empty_vertex_permutations ||
        !std::next_permutation(
            m_empty_source_vertices.begin(), m_empty_source_vertices.end())) {
      return exact_map_result_for_permuted_vertices;
    }
  }
}

PartialMappingLookup::Parameters::Parameters()
    : max_number_of_empty_vertex_permutations(10) {}

PartialMappingLookup::Parameters& PartialMappingLookup::get_parameters() {
  return m_parameters;
}

}  // namespace tsa_internal
}  // namespace tket
