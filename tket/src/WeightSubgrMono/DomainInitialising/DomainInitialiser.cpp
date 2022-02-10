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

#include "WeightSubgrMono/DomainInitialising/DomainInitialiser.hpp"

#include <algorithm>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/DomainInitialising/DistanceCounts.hpp"
#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DomainInitialiser::Parameters::Parameters() : max_path_length(10) {}

bool DomainInitialiser::full_initialisation(
    PossibleAssignments& possible_assignments,
    const std::vector<VertexWSM>& pattern_vertices,
    const NeighboursData& pattern_neighbours_data,
    const std::vector<VertexWSM>& target_vertices,
    const NeighboursData& target_neighbours_data, const Parameters& params) {
  return degree_sequence_initialisation(
             possible_assignments, pattern_vertices, pattern_neighbours_data,
             target_vertices, target_neighbours_data) &&

         distance_counts_reduction(
             possible_assignments, pattern_neighbours_data,
             target_neighbours_data, params) &&

         triangle_counts_reduction(
             possible_assignments, pattern_neighbours_data,
             target_neighbours_data) &&

         square_counts_reduction(
             possible_assignments, pattern_neighbours_data,
             target_neighbours_data) &&

         alldiff_reduction(possible_assignments);
}

bool DomainInitialiser::degree_sequence_initialisation(
    PossibleAssignments& possible_assignments,
    const std::vector<VertexWSM>& pattern_vertices,
    const NeighboursData& pattern_neighbours_data,
    const std::vector<VertexWSM>& target_vertices,
    const NeighboursData& target_neighbours_data) {
  possible_assignments.clear();
  for (auto pv : pattern_vertices) {
    possible_assignments[pv];
  }
  // Get the target degree sequences.
  std::vector<std::vector<std::size_t>> target_degree_sequences(
      target_vertices.size());
  // Sort by decreasing sequence length; nonstable sort is fine.
  std::vector<std::size_t> sequence_indices(target_vertices.size());
  for (unsigned ii = 0; ii < target_vertices.size(); ++ii) {
    target_degree_sequences[ii] =
        target_neighbours_data.get_sorted_degree_sequence_expensive(
            target_vertices[ii]);
    sequence_indices[ii] = ii;
  }
  std::sort(
      sequence_indices.begin(), sequence_indices.end(),
      [&target_degree_sequences](std::size_t lhs, std::size_t rhs) -> bool {
        return target_degree_sequences[lhs].size() >
               target_degree_sequences[rhs].size();
      });

  std::vector<std::size_t> pattern_sequence;
  for (auto& entry : possible_assignments) {
    pattern_sequence =
        pattern_neighbours_data.get_sorted_degree_sequence_expensive(
            entry.first);
    for (auto index : sequence_indices) {
      const auto& target_sequence = target_degree_sequences[index];
      if (target_sequence.size() < pattern_sequence.size()) {
        break;
      }
      if (FilterUtils::compatible_sorted_degree_sequences(
              pattern_sequence, target_sequence)) {
        entry.second.insert(target_vertices[index]);
      }
    }
    if (entry.second.empty()) {
      return false;
    }
  }
  return true;
}

bool DomainInitialiser::distance_counts_reduction(
    PossibleAssignments& possible_assignments,
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data, const Parameters& params) {
  auto& tv_to_erase = m_work_vector;

  DistanceCounts p_distance_counts;
  std::map<VertexWSM, DistanceCounts> target_distance_counts;

  const auto get_t_counts = [&target_distance_counts, &target_neighbours_data](
                                VertexWSM tv) -> DistanceCounts& {
    const auto iter = target_distance_counts.find(tv);
    if (iter == target_distance_counts.end()) {
      DistanceCounts& counts = target_distance_counts[tv];
      counts.initialise(target_neighbours_data, tv);
      return counts;
    }
    return iter->second;
  };

  for (auto& entry : possible_assignments) {
    p_distance_counts.initialise(pattern_neighbours_data, entry.first);

    // We're going to need the full depth for P vertices,
    // as we will keep trying to get "false": we won't be happy
    // with it returning "true" if we could make it more stringent.
    while (p_distance_counts.size() < params.max_path_length) {
      if (!p_distance_counts.push_back(pattern_neighbours_data)) {
        break;
      }
    }
    tv_to_erase.clear();
    auto& domain = entry.second;
    for (auto tv : domain) {
      DistanceCounts& t_counts = get_t_counts(tv);

      for (;;) {
        if (p_distance_counts.test_against_target(t_counts)) {
          break;
        }
        // We've got a "false" result; but maybe this is just because
        // we haven't grown the T-subgraph enough?
        if (t_counts.size() >= p_distance_counts.size()) {
          // It's already as big as it needs to be,
          // so the false result is correct.
          tv_to_erase.push_back(tv);
          break;
        }
        if (!t_counts.push_back(target_neighbours_data)) {
          // We can't expand any more; so the "false" reading stands.
          tv_to_erase.push_back(tv);
          break;
        }
      }
    }
    if (tv_to_erase.size() == domain.size()) {
      return false;
    }
    for (auto tv : tv_to_erase) {
      domain.erase(tv);
    }
  }
  return true;
}

bool DomainInitialiser::alldiff_reduction(
    PossibleAssignments& possible_assignments) {
  // We only push back to this at the moment when
  // the domain drops down from 2 to 1.
  // Domains cannot grow, only decrease.
  // Therefore, no duplicates.
  m_assigned_vertices.clear();

  for (const auto& entry : possible_assignments) {
    switch (entry.second.size()) {
      case 0:
        return false;
      case 1:
        m_assigned_vertices.push_back(entry.first);
      default:
        break;
    }
  }
  std::size_t n_vertices_processed = 0;

  while (n_vertices_processed < m_assigned_vertices.size()) {
    const auto pv = m_assigned_vertices[n_vertices_processed];
    ++n_vertices_processed;
    const auto& domain = possible_assignments.at(pv);
    if (domain.size() != 1) {
      return false;
    }
    const auto tv = *domain.cbegin();
    for (auto& entry : possible_assignments) {
      if (entry.first == pv) {
        continue;
      }
      if (entry.second.erase(tv) == 1) {
        // We did erase.
        if (entry.second.empty()) {
          return false;
        }
        if (entry.second.size() == 1) {
          m_assigned_vertices.push_back(entry.first);
        }
      }
    }
  }
  std::sort(m_assigned_vertices.begin(), m_assigned_vertices.end());
  return true;
}

const std::vector<VertexWSM>& DomainInitialiser::get_assigned_vertices() const {
  return m_assigned_vertices;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
