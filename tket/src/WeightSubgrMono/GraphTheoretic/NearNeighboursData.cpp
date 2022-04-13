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

#include "WeightSubgrMono/GraphTheoretic/NearNeighboursData.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NearNeighboursData::NearNeighboursData(const NeighboursData& ndata)
    : m_ndata(ndata) {}

void NearNeighboursData::fill_counts_vector(
    VertexWSM vv, unsigned max_distance,
    std::vector<std::size_t>& counts_vector) {
  counts_vector.reserve(max_distance);
  counts_vector.clear();
  counts_vector.push_back(m_ndata.get_degree(vv));
  for (unsigned distance = 2; distance <= max_distance; ++distance) {
    const auto& vertices = get_vertices_at_distance(vv, distance);
    counts_vector.push_back(vertices.size());
    if (counts_vector.back() == 0) {
      // Once a zero occurs, no point filling in more zeros.
      break;
    }
  }
}

bool NearNeighboursData::test_against_target(
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
  // We simply push each p-vertex at distance d into a t-vertex
  // of smallest possible distance <= d.
  unsigned pattern_index = 0;
  unsigned target_index = 0;

  // The number of p-vertices we must erase.
  auto remaining_p_count = p_counts[0];

  // How many t-vertex holes remain at the lowest possible level
  // for the p-vertices to fall into.
  auto remaining_t_holes = t_counts[0];

  for (;;) {
    while (remaining_p_count == 0) {
      ++pattern_index;
      if (pattern_index >= p_counts.size()) {
        return true;
      }
      remaining_p_count = p_counts[pattern_index];
    }

    TKET_ASSERT(remaining_p_count > 0);

    // Still some p-vertices need to be cleared.
    // So we MUST get some target holes.
    while (remaining_t_holes == 0) {
      ++target_index;
      if (target_index >= t_counts.size() || target_index > pattern_index) {
        return false;
      }
      remaining_t_holes = t_counts[target_index];
    }

    TKET_ASSERT(remaining_t_holes > 0);
    TKET_ASSERT(target_index <= pattern_index);

    // Now, clear as many p-vertices as we can.
    if (remaining_p_count <= remaining_t_holes) {
      remaining_t_holes -= remaining_p_count;
      remaining_p_count = 0;
    } else {
      // We can only clear some p-vertices.
      remaining_p_count -= remaining_t_holes;
      remaining_t_holes = 0;
    }
  }
  return true;
}

static void fill_neighbours_of_neighbours(
    const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours_and_weights,
    const NeighboursData& ndata, std::set<VertexWSM>& vertices_workset,
    std::vector<VertexWSM>& vertices) {
  vertices_workset.clear();
  for (const auto& entry : neighbours_and_weights) {
    const VertexWSM& vv1 = entry.first;
    for (const auto& other_entry : ndata.get_neighbours_and_weights(vv1)) {
      const VertexWSM& vv2 = other_entry.first;
      if (NeighboursData::binary_search(vv2, neighbours_and_weights)) {
        continue;
      }
      vertices_workset.insert(vv2);
    }
  }
  vertices = {vertices_workset.cbegin(), vertices_workset.cend()};
}

// result[i] is for distance i+2.
// We've already filled in result[0] at least, and resized to the desired size.
// Now fill in all the result[i].
static void fill_more_distant_vertices(
    std::vector<std::vector<VertexWSM>>& result, unsigned old_size,
    const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours_and_weights,
    std::set<VertexWSM>& vertices_workset, const NeighboursData& ndata) {
  TKET_ASSERT(old_size > 0);
  TKET_ASSERT(old_size <= result.size());

  for (unsigned index = old_size; index < result.size(); ++index) {
    vertices_workset.clear();
    for (auto vv_prev : result[index - 1]) {
      for (const auto& entry : ndata.get_neighbours_and_weights(vv_prev)) {
        const VertexWSM& vv_new = entry.first;
        if (NeighboursData::binary_search(vv_new, neighbours_and_weights)) {
          continue;
        }
        // We must also search previous entries.
        bool seen_already = false;
        for (unsigned jj = 0; jj < index; ++jj) {
          if (std::binary_search(
                  result[jj].cbegin(), result[jj].cend(), vv_new)) {
            seen_already = true;
            break;
          }
        }
        if (!seen_already) {
          vertices_workset.insert(vv_new);
        }
      }
    }
    if (vertices_workset.empty()) {
      // no point in filling in more empty vertex lists;
      // they're already empty.
      break;
    }
    result[index] = {vertices_workset.cbegin(), vertices_workset.cend()};
  }
}

const std::vector<VertexWSM>& NearNeighboursData::get_vertices_at_distance(
    VertexWSM vv, unsigned max_distance) {
  TKET_ASSERT(max_distance >= 2);
  const auto index = max_distance - 2;
  auto& results_for_this_vertex = m_data[vv];
  unsigned old_size = results_for_this_vertex.size();

  if (index >= old_size) {
    // "results_for_this_vertex" is too small, so we must fill in more entries.
    // Give it the correct size, prior to filling.
    results_for_this_vertex.resize(index + 1);
    const auto& neighbours_and_weights = m_ndata.get_neighbours_and_weights(vv);

    if (old_size == 0) {
      // A special case. We must fill in distance 2 data.
      fill_neighbours_of_neighbours(
          neighbours_and_weights, m_ndata, m_vertices_workset,
          results_for_this_vertex[0]);
      old_size = 1;
    }
    fill_more_distant_vertices(
        results_for_this_vertex, old_size, neighbours_and_weights,
        m_vertices_workset, m_ndata);
  }
  return results_for_this_vertex.at(index);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
