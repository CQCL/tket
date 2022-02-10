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

#include "WeightSubgrMono/DomainInitialising/DistanceCounts.hpp"

#include <algorithm>

#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

void DistanceCounts::initialise(
    const NeighboursData& neighbours_data, VertexWSM v) {
  m_counts.clear();
  m_vertices_seen.clear();
  m_current_frontier.clear();
  m_vertices_seen.insert(v);
  m_current_frontier.insert(v);
  push_back(neighbours_data);
}

bool DistanceCounts::push_back(const NeighboursData& neighbours_data) {
  if (!m_counts.empty() && m_counts.back() == 0) {
    return false;
  }
  // This will be the new frontier.
  m_work_set.clear();
  for (auto frontier_v : m_current_frontier) {
    const auto& new_neighbours =
        neighbours_data.get_neighbours_and_weights(frontier_v);
    for (const auto& entry : new_neighbours) {
      const auto& new_v = entry.first;
      if (m_vertices_seen.count(new_v) == 0) {
        m_vertices_seen.insert(new_v);
        m_work_set.insert(new_v);
      }
    }
  }
  m_counts.push_back(m_work_set.size());
  if (m_counts.back() == 0) {
    return false;
  }
  m_work_set.swap(m_current_frontier);
  return true;
}

const std::vector<std::size_t>& DistanceCounts::get_counts() const {
  return m_counts;
}

std::size_t DistanceCounts::size() const { return m_counts.size(); }

bool DistanceCounts::test_against_target(const DistanceCounts& other) const {
  return test_against_target(m_counts, other.m_counts);
}

bool DistanceCounts::test_against_target(
    const std::vector<std::size_t>& p_counts,
    const std::vector<std::size_t>& t_counts) {
  if (p_counts.empty()) {
    return true;
  }
  if (t_counts.empty()) {
    for (auto count : p_counts) {
      if (count > 0) {
        return false;
      }
    }
    return true;
  }
  unsigned pattern_index = 0;
  unsigned target_index = 0;
  auto remaining_p_count_at_this_level = p_counts[0];
  auto remaining_t_holes_at_this_level = t_counts[0];

  for (;;) {
    if (remaining_p_count_at_this_level <= remaining_t_holes_at_this_level) {
      // We can clear this p-level.
      ++pattern_index;
      if (pattern_index >= p_counts.size()) {
        return true;
      }
      remaining_t_holes_at_this_level -= remaining_p_count_at_this_level;
      remaining_p_count_at_this_level = p_counts[pattern_index];
      if (remaining_p_count_at_this_level == 0) {
        return true;
      }
      if (remaining_t_holes_at_this_level != 0) {
        continue;
      }

      // The t-vertices have also been cleared.
      // But some new p-vertices still need to be cleared.
      ++target_index;
      if (target_index >= t_counts.size()) {
        return false;
      }
      remaining_t_holes_at_this_level = t_counts[target_index];
      continue;
    }

    // Here, there are strictly more p than can be paired off with t;
    // so the t-vertices are cleared instead.
    remaining_p_count_at_this_level -= remaining_t_holes_at_this_level;
    ++target_index;
    if (target_index >= t_counts.size() || target_index > pattern_index) {
      break;
    }
    remaining_t_holes_at_this_level = t_counts[target_index];
  }
  return false;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
