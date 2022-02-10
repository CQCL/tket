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

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/DomainInitialising/DomainInitialiser.hpp"
#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace {

struct LengthTwoData {
  // KEY: a vertex which is the end of a length-2 path
  // starting from v.
  // VALUE: the number of distinct middle vertices we can pass through
  //   to reach the end.
  // std::map<VertexWSM, std::size_t> length_two_map;
  PossibleAssignments length_two_detailed_map;

  // TODO: more expensive but more discriminating filters
  // would look at the detailed vertex sets and consider matching.
  //
  // Fills and returns by reference a sorted vector of values.
  void initialise(
      VertexWSM v, const NeighboursData& neighbours_data,
      std::vector<std::size_t>& values) {
    const auto& neighbours = neighbours_data.get_neighbours_and_weights(v);
    length_two_detailed_map.clear();
    for (const auto& entry : neighbours) {
      const auto& middle_v = entry.first;
      const auto& final_edges =
          neighbours_data.get_neighbours_and_weights(middle_v);
      for (const auto& end_entry : final_edges) {
        const auto& end_v = end_entry.first;
        if (end_v == v) {
          continue;
        }
        length_two_detailed_map[end_v].insert(middle_v);
      }
    }
    values.clear();
    values.reserve(length_two_detailed_map.size());
    for (const auto& entry : length_two_detailed_map) {
      values.push_back(entry.second.size());
    }
    std::sort(values.begin(), values.end());
  }
};
}  // namespace

static bool first_embeds_into_second(
    const std::vector<std::size_t>& lhs, const std::vector<std::size_t>& rhs) {
  return FilterUtils::compatible_sorted_degree_sequences(lhs, rhs);
}

// A more thorough reduction would look at the detailed
// midpoint vertices in each length-2 path and try to map them
// using bipartite matching
// (or rather, show that they cannot be mapped).
bool DomainInitialiser::square_counts_reduction(
    PossibleAssignments& possible_assignments,
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data) {
  LengthTwoData calculator;
  std::map<VertexWSM, std::vector<std::size_t>> target_data_map;
  auto& t_vertices_to_erase = m_work_vector;
  std::vector<std::size_t> pattern_vertex_values;

  for (auto& entry : possible_assignments) {
    t_vertices_to_erase.clear();
    calculator.initialise(
        entry.first, pattern_neighbours_data, pattern_vertex_values);
    auto& domain = entry.second;
    for (auto tv : domain) {
      auto& target_values = target_data_map[tv];
      if (target_values.empty()) {
        // It must be uninitialised; isolated vertices have been removed
        calculator.initialise(tv, target_neighbours_data, target_values);
      }
      if (!first_embeds_into_second(pattern_vertex_values, target_values)) {
        t_vertices_to_erase.push_back(tv);
      }
    }
    if (t_vertices_to_erase.size() == domain.size()) {
      return false;
    }
    for (auto tv : t_vertices_to_erase) {
      domain.erase(tv);
    }
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
