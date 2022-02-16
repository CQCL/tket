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

#include "TokenSwapping/ExactMappingLookup.hpp"

#include <algorithm>

#include "TokenSwapping/FilteredSwapSequences.hpp"
#include "TokenSwapping/GeneralFunctions.hpp"
#include "TokenSwapping/SwapConversion.hpp"
#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

const ExactMappingLookup::Result& ExactMappingLookup::operator()(
    const VertexMapping& desired_mapping, const vector<Swap>& edges,
    unsigned max_number_of_swaps) {
  m_result.success = false;
  m_result.too_many_vertices = desired_mapping.size() > 6;
  m_result.swaps.clear();
  if (m_result.too_many_vertices) {
    return m_result;
  }
  return improve_upon_existing_result(
      desired_mapping, edges, max_number_of_swaps);
}

const ExactMappingLookup::Result&
ExactMappingLookup::improve_upon_existing_result(
    const VertexMapping& desired_mapping, const vector<Swap>& edges,
    unsigned max_number_of_swaps) {
  max_number_of_swaps = std::min(max_number_of_swaps, 16u);
  const auto& relabelling = m_relabeller(desired_mapping);

  if (relabelling.identity) {
    // This beats whatever was there before,
    // whether or not it was successful.
    m_result.success = true;
    m_result.too_many_vertices = false;
    m_result.swaps.clear();
    return m_result;
  }
  if (relabelling.too_many_vertices) {
    // We cannot get a new result, so just return the existing one, whether or
    // not it succeeded.
    if (!m_result.success) {
      m_result.too_many_vertices = true;
    }
    return m_result;
  }
  TKET_ASSERT(relabelling.permutation_hash != 0);
  {
    const bool size_match = relabelling.new_to_old_vertices.size() ==
                            relabelling.old_to_new_vertices.size();
    TKET_ASSERT(size_match);
  }
  TKET_ASSERT(relabelling.new_to_old_vertices.size() >= 2);

  fill_result_from_table(relabelling, edges, max_number_of_swaps);
  return m_result;
}

void ExactMappingLookup::fill_result_from_table(
    const CanonicalRelabelling::Result& relabelling_result,
    const vector<Swap>& old_edges, unsigned max_number_of_swaps) {
  if (m_result.success) {
    if (m_result.swaps.empty()) {
      return;
    }
    max_number_of_swaps =
        std::min<unsigned>(max_number_of_swaps, m_result.swaps.size() - 1);
    if (max_number_of_swaps == 0) {
      return;
    }
  } else {
    m_result.swaps.clear();
  }
  SwapConversion::EdgesBitset new_edges_bitset = 0;

  for (auto old_edge : old_edges) {
    const auto new_v1_opt = get_optional_value(
        relabelling_result.old_to_new_vertices, old_edge.first);
    if (!new_v1_opt) {
      continue;
    }
    const auto new_v2_opt = get_optional_value(
        relabelling_result.old_to_new_vertices, old_edge.second);
    if (!new_v2_opt) {
      continue;
    }
    const auto new_v1 = new_v1_opt.value();
    const auto new_v2 = new_v2_opt.value();
    TKET_ASSERT(new_v1 <= 5);
    TKET_ASSERT(new_v2 <= 5);
    new_edges_bitset |= SwapConversion::get_edges_bitset(
        SwapConversion::get_hash_from_swap(get_swap(new_v1, new_v2)));
  }

  const FilteredSwapSequences::SingleSequenceData table_result(
      relabelling_result.permutation_hash, new_edges_bitset,
      max_number_of_swaps);

  TKET_ASSERT(table_result.number_of_swaps > 0);
  if (table_result.number_of_swaps > max_number_of_swaps) {
    // No result in the table.
    return;
  }
  TKET_ASSERT(table_result.edges_bitset != 0);
  TKET_ASSERT(table_result.swaps_code > 0);

  m_result.success = true;
  m_result.swaps.clear();
  auto swaps_code_copy = table_result.swaps_code;
  while (swaps_code_copy != 0) {
    const auto& new_swap =
        SwapConversion::get_swap_from_hash(swaps_code_copy & 0xF);
    swaps_code_copy >>= 4;
    m_result.swaps.push_back(get_swap(
        relabelling_result.new_to_old_vertices.at(new_swap.first),
        relabelling_result.new_to_old_vertices.at(new_swap.second)));
  }
  TKET_ASSERT(m_result.swaps.size() <= 16);
}

}  // namespace tsa_internal
}  // namespace tket
