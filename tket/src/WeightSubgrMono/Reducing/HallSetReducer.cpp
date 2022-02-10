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

#include "WeightSubgrMono/Reducing/HallSetReducer.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Searching/SearchBranch.hpp"
#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// We put smaller domains at the TOP (back) of a vector, so we can pop_back
// the smaller domains.
bool HallSetReducer::VariableData::operator<(const VariableData& other) const {
  return domain_size > other.domain_size ||
         (domain_size == other.domain_size && vertex < other.vertex);
}

bool HallSetReducer::reduce(
    SearchNodeWrapper& search_node_wrapper, SearchBranch& branch) const {
  m_domain_sizes_and_vertices.clear();

  // Initial fill of the data.
  for (const auto& entry :
       search_node_wrapper.get().pattern_v_to_possible_target_v) {
    m_domain_sizes_and_vertices.emplace_back();
    m_domain_sizes_and_vertices.back().vertex = entry.first;
    m_domain_sizes_and_vertices.back().domain_size = entry.second.size();

    if (m_domain_sizes_and_vertices.back().domain_size == 0) {
      return false;
    }
  }
  std::sort(
      m_domain_sizes_and_vertices.begin(), m_domain_sizes_and_vertices.end());

  for (;;) {
    const auto old_size = m_domain_sizes_and_vertices.size();
    if (!find_and_remove_top_hall_set_block(search_node_wrapper, branch)) {
      return false;
    }
    if (old_size == m_domain_sizes_and_vertices.size()) {
      break;
    }
    TKET_ASSERT(old_size > m_domain_sizes_and_vertices.size());
  }
  return true;
}

bool HallSetReducer::find_and_remove_top_hall_set_block(
    SearchNodeWrapper& search_node_wrapper, SearchBranch& branch) const {
  m_combined_domains.clear();
  std::size_t hall_set_size = 0;
  bool hall_set_found = false;

  for (auto citer = m_domain_sizes_and_vertices.crbegin();
       citer != m_domain_sizes_and_vertices.crend(); ++citer) {
    ++hall_set_size;
    const auto& domains =
        search_node_wrapper.get().pattern_v_to_possible_target_v;

    const auto domain_citer = domains.find(citer->vertex);
    if (domain_citer == domains.cend()) {
      return false;
    }
    for (auto tv : domain_citer->second) {
      m_combined_domains.insert(tv);
    }
    if (m_combined_domains.size() < hall_set_size) {
      return false;
    }
    if (m_combined_domains.size() == hall_set_size) {
      hall_set_found = true;
      break;
    }
  }
  if (!hall_set_found) {
    return true;
  }
  return remove_top_hall_set_block(hall_set_size, search_node_wrapper, branch);
}

bool HallSetReducer::remove_top_hall_set_block(
    std::size_t hall_set_size, SearchNodeWrapper& search_node_wrapper,
    SearchBranch& branch) const {
  m_domain_sizes_and_vertices.resize(
      m_domain_sizes_and_vertices.size() - hall_set_size);
  bool needs_reorder = false;
  for (auto& entry : m_domain_sizes_and_vertices) {
    const auto new_size = search_node_wrapper.remove_elements_from_domain(
        entry.vertex, m_combined_domains, branch.get_assignments_mutable());
    if (new_size == 0) {
      return false;
    }
    if (entry.domain_size != new_size) {
      needs_reorder = true;
      entry.domain_size = new_size;
    }
  }
  if (needs_reorder) {
    std::sort(
        m_domain_sizes_and_vertices.begin(), m_domain_sizes_and_vertices.end());
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
