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

static void fill_neighbours_of_neighbours(
    VertexWSM root_vertex,
    const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours_and_weights,
    const NeighboursData& ndata, std::set<VertexWSM>& vertices_workset,
    std::vector<VertexWSM>& vertices) {
  vertices_workset.clear();
  for (const auto& entry : neighbours_and_weights) {
    const VertexWSM& vv1 = entry.first;
    for (const auto& other_entry : ndata.get_neighbours_and_weights(vv1)) {
      const VertexWSM& vv2 = other_entry.first;
      if (vv2 == root_vertex ||
          NeighboursData::binary_search(vv2, neighbours_and_weights)) {
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
      TKET_ASSERT(result[index].empty());
      for (const auto& entry : ndata.get_neighbours_and_weights(vv_prev)) {
        const VertexWSM& vv_new = entry.first;
        // vv_prev is at distance d from the root vertex;
        // thus neighbours of vv_prev are at distance d-1, d, d+1.
        // So we only need to check two vectors (for d-1,d).
        if (std::binary_search(
                result[index - 1].cbegin(), result[index - 1].cend(), vv_new)) {
          continue;
        }
        if (index == 1) {
          if (NeighboursData::binary_search(vv_new, neighbours_and_weights)) {
            continue;
          }
        } else {
          if (std::binary_search(
                  result[index - 2].cbegin(), result[index - 2].cend(),
                  vv_new)) {
            continue;
          }
        }
        vertices_workset.insert(vv_new);
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
    VertexWSM vv, unsigned distance) {
  TKET_ASSERT(distance >= 2);
  const auto index = distance - 2;
  auto& results_for_this_vertex = m_data[vv].vertices_at_distance;
  unsigned old_size = results_for_this_vertex.size();

  if (index >= old_size) {
    // "results_for_this_vertex" is too small, so we must fill in more entries.
    // Give it the correct size, prior to filling.
    results_for_this_vertex.resize(index + 1);
    const auto& neighbours_and_weights = m_ndata.get_neighbours_and_weights(vv);

    if (old_size == 0) {
      // A special case. We must fill in distance 2 data.
      fill_neighbours_of_neighbours(
          vv, neighbours_and_weights, m_ndata, m_vertices_workset,
          results_for_this_vertex[0]);
      old_size = 1;
    }
    fill_more_distant_vertices(
        results_for_this_vertex, old_size, neighbours_and_weights,
        m_vertices_workset, m_ndata);
  }
  return results_for_this_vertex.at(index);
}

std::size_t NearNeighboursData::get_n_vertices_at_max_distance(
    VertexWSM vv, unsigned max_distance) {
  switch (max_distance) {
    case 0:
      return 0;
    case 1:
      return m_ndata.get_degree(vv);
    default:
      break;
  }
  auto& sizes_list = m_data[vv].n_vertices_at_max_distance;
  const unsigned index = max_distance - 2;
  if (index < sizes_list.size()) {
    return sizes_list[index];
  }
  auto old_size = sizes_list.size();
  sizes_list.resize(index + 1);
  if (old_size == 0) {
    sizes_list[0] =
        m_ndata.get_degree(vv) + get_vertices_at_distance(vv, 2).size();
    ++old_size;
  }
  for (unsigned ii = old_size; ii <= index; ++ii) {
    sizes_list[ii] =
        sizes_list[ii - 1] + get_vertices_at_distance(vv, ii + 2).size();
  }
  return sizes_list[index];
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
