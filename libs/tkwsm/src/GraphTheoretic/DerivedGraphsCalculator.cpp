// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/DerivedGraphsCalculator.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

void DerivedGraphsCalculator::fill_mid_vertices_for_length_two_paths(
    const NeighboursData& ndata, VertexWSM v) {
  m_mid_vertices_for_length_two_paths.clear();
  const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours =
      ndata.get_neighbours_and_weights(v);
  for (const std::pair<VertexWSM, WeightWSM>& entry : neighbours) {
    const VertexWSM& v1 = entry.first;
    for (const std::pair<VertexWSM, WeightWSM>& v2_entry :
         ndata.get_neighbours_and_weights(v1)) {
      const VertexWSM& v2 = v2_entry.first;
      if (v2 == v) {
        continue;
      }
      // Notice that v1 occurs in sorted order!
      m_mid_vertices_for_length_two_paths[v2].push_back(v1);
    }
  }
}

void DerivedGraphsCalculator::fill_d2_neighbours_and_counts(
    DerivedGraphStructs::NeighboursAndCounts& depth_2_neighbours_and_counts) {
  depth_2_neighbours_and_counts.clear();
  for (const auto& entry : m_mid_vertices_for_length_two_paths) {
    const VertexWSM& v2 = entry.first;
    const auto& v1_set = entry.second;
    // Notice, the v2 vertices already are in sorted order!
    TKET_ASSERT(is_sorted_and_unique(v1_set));
    depth_2_neighbours_and_counts.emplace_back(v2, v1_set.size());
  }
}

void DerivedGraphsCalculator::fill_d3_neighbours_and_counts_map(
    const NeighboursData& ndata) {
  m_depth_3_neighbours_and_counts_map.clear();
  for (const auto& entry : m_mid_vertices_for_length_two_paths) {
    const VertexWSM& v2 = entry.first;
    const auto& v1_set = entry.second;

    // Now, v--v1--v2--v3 will be either a path, or a triangle (if v3=v).
    // v--v1--v2--v3 is constructed so that each vertex is a neighbour
    // of the previous. Thus it's a valid path as long as v2 != v
    // (already checked and excluded) and v3 != v1.
    // Note that v3 IS allowed to be a neighbour of v,
    // as long as it isn't v1.
    // Let   |{v1}|=N.
    for (const std::pair<VertexWSM, WeightWSM>& v3_entry :
         ndata.get_neighbours_and_weights(v2)) {
      // The weight is simply ignored.
      const VertexWSM& v3 = v3_entry.first;
      if (std::binary_search(v1_set.cbegin(), v1_set.cend(), v3)) {
        // It's a path of the form  v--v1--v2--(v1)'.
        // So we contribute N-1, since  v1 != (v1)' is the only restriction.
        m_depth_3_neighbours_and_counts_map[v3] += v1_set.size() - 1;
      } else {
        // It's a path of the form  v--v1--v2--v3, so contributes the whole N.
        m_depth_3_neighbours_and_counts_map[v3] += v1_set.size();
      }
    }
  }
}

void DerivedGraphsCalculator::fill_remaining_d3_data(
    VertexWSM v, DerivedGraphStructs::Count& triangle_count,
    DerivedGraphStructs::NeighboursAndCounts& depth_3_neighbours_and_counts) {
  depth_3_neighbours_and_counts.clear();
  triangle_count = 0;
  for (const auto& entry : m_depth_3_neighbours_and_counts_map) {
    if (entry.first == v) {
      triangle_count = entry.second;
    } else {
      depth_3_neighbours_and_counts.emplace_back(entry);
    }
  }
}

void DerivedGraphsCalculator::fill(
    const NeighboursData& ndata, VertexWSM v,
    DerivedGraphStructs::Count& triangle_count,
    DerivedGraphStructs::NeighboursAndCounts& depth_2_neighbours_and_counts,
    DerivedGraphStructs::NeighboursAndCounts& depth_3_neighbours_and_counts) {
  fill_mid_vertices_for_length_two_paths(ndata, v);
  fill_d2_neighbours_and_counts(depth_2_neighbours_and_counts);
  fill_d3_neighbours_and_counts_map(ndata);
  fill_remaining_d3_data(v, triangle_count, depth_3_neighbours_and_counts);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
