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
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static std::size_t get_triangles_count(
    VertexWSM v, const NeighboursData& neighbours_data) {
  std::size_t count = 0;
  const auto& neighbours = neighbours_data.get_neighbours_and_weights(v);
  for (unsigned ii = 0; ii < neighbours.size(); ++ii) {
    for (unsigned jj = ii + 1; jj < neighbours.size(); ++jj) {
      if (neighbours_data.get_edge_weight_opt(
              neighbours[ii].first, neighbours[jj].first)) {
        ++count;
      }
    }
  }
  return count;
}

// A more thorough reduction would look at the detailed vertices
// within the triangles and try to map them
// (or rather, show that they cannot be mapped).
bool DomainInitialiser::triangle_counts_reduction(
    PossibleAssignments& possible_assignments,
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data) {
  auto& t_vertices_to_erase = m_work_vector;
  std::map<VertexWSM, std::size_t> target_triangle_counts;

  for (auto& entry : possible_assignments) {
    t_vertices_to_erase.clear();
    const auto count =
        get_triangles_count(entry.first, pattern_neighbours_data);
    auto& domain = entry.second;
    for (auto tv : domain) {
      const auto target_count_opt =
          get_optional_value(target_triangle_counts, tv);
      if (target_count_opt) {
        if (count <= target_count_opt) {
          continue;
        }
      } else {
        const auto target_count =
            get_triangles_count(tv, target_neighbours_data);
        target_triangle_counts[tv] = target_count;
        if (count <= target_count) {
          continue;
        }
      }
      t_vertices_to_erase.push_back(tv);
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
