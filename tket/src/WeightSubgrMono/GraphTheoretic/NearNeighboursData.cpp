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

NearNeighboursData::NearNeighboursData(const NeighboursData& ndata, Type type)
    : m_ndata(ndata), m_type(type) {
  m_data.resize(ndata.get_number_of_nonisolated_vertices());
}

static void fill_neighbours_of_neighbours(
    VertexWSM root_vertex,
    const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours_and_weights,
    const NeighboursData& ndata, std::set<VertexWSM>& vertices_workset,
    std::vector<VertexWSM>& vertices) {
  vertices_workset.clear();
  for (const std::pair<VertexWSM, WeightWSM>& entry : neighbours_and_weights) {
    const VertexWSM& vv1 = entry.first;
    for (const std::pair<VertexWSM, WeightWSM>& other_entry :
         ndata.get_neighbours_and_weights(vv1)) {
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
    for (VertexWSM vv_prev : result[index - 1]) {
      TKET_ASSERT(result[index].empty());
      for (const std::pair<VertexWSM, WeightWSM>& entry :
           ndata.get_neighbours_and_weights(vv_prev)) {
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
  const unsigned index = distance - 2;
  std::vector<std::vector<VertexWSM>>& results_for_this_vertex =
      m_data[vv].vertices_at_distance;
  unsigned old_size = results_for_this_vertex.size();

  if (index >= old_size) {
    // "results_for_this_vertex" is too small, so we must fill in more entries.
    // Give it the correct size, prior to filling.
    results_for_this_vertex.resize(index + 1);
    const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours_and_weights =
        m_ndata.get_neighbours_and_weights(vv);

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
  std::vector<std::size_t>& sizes_list = m_data[vv].n_vertices_at_max_distance;
  const unsigned index = max_distance - 2;
  if (index < sizes_list.size()) {
    return sizes_list[index];
  }
  unsigned old_size = sizes_list.size();
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

const FilterUtils::DegreeCounts& NearNeighboursData::get_degree_counts(
    VertexWSM vv, unsigned distance) {
  TKET_ASSERT(distance >= 2);
  std::vector<FilterUtils::DegreeCounts>& counts_list =
      m_data[vv].degree_counts_for_distance;
  const unsigned index = distance - 2;
  if (index < counts_list.size()) {
    return counts_list[index];
  }
  // We don't have data up to this distance, for this vertex.
  // Ensure first that we HAVE calculated the vertices up to this distance.
  get_vertices_at_distance(vv, distance);

  // Notice that the vertices, and degree counts here, are stored separately,
  // but on the SAAME value m_data[vv],
  // so references aren't invalidated.
  const std::vector<std::vector<VertexWSM>>& vertices_lists =
      m_data[vv].vertices_at_distance;
  TKET_ASSERT(index < vertices_lists.size());

  unsigned old_size = counts_list.size();
  counts_list.resize(index + 1);
  m_work_map.clear();
  if (old_size == 0) {
    // We have to fill in distance 2 data as the initial "seed".
    for (VertexWSM distance_2_vertex : vertices_lists[0]) {
      // C++ standard guarantees it will be initialised with 0.
      m_work_map[m_ndata.get_degree(distance_2_vertex)] += 1;
    }
    if (m_type == Type::TARGET_GRAPH) {
      // We ALSO must fill in distance 1 data, i.e. neighbours.
      for (const std::pair<VertexWSM, WeightWSM>& entry :
           m_ndata.get_neighbours_and_weights(vv)) {
        m_work_map[m_ndata.get_degree(entry.first)] += 1;
      }
    }
    counts_list[0].reserve(m_work_map.size());
    for (const std::pair<const std::size_t, std::size_t>& entry : m_work_map) {
      counts_list[0].emplace_back(entry);
    }
    // Simply pretend that we'd already filled it in.
    ++old_size;
  }

  for (unsigned index_to_fill = old_size; index_to_fill <= index;
       ++index_to_fill) {
    switch (m_type) {
      case Type::PATTERN_GRAPH:
        m_work_map.clear();
        break;
      case Type::TARGET_GRAPH:
        // We'll maintain the data within m_work_map and add to it
        // as we increase the distances.
        if (m_work_map.empty()) {
          // Need the initial fill of the previous data.
          const FilterUtils::DegreeCounts& previous_counts =
              counts_list[index_to_fill - 1];
          for (const std::pair<std::size_t, std::size_t>& entry :
               previous_counts) {
            m_work_map.emplace(entry);
          }
        }
        break;
    }
    const std::vector<VertexWSM>& vertices = vertices_lists[index_to_fill];
    for (VertexWSM other_vertex : vertices) {
      m_work_map[m_ndata.get_degree(other_vertex)] += 1;
    }
    // Now dump it out from the map back into the vector, to be stored.
    counts_list[index_to_fill].reserve(m_work_map.size());
    for (const std::pair<const std::size_t, std::size_t>& entry : m_work_map) {
      counts_list[index_to_fill].emplace_back(entry);
    }
    TKET_ASSERT(counts_list[index_to_fill].size() == m_work_map.size());
  }
  return counts_list[index];
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
