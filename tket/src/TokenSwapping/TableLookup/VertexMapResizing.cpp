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

#include "TokenSwapping/VertexMapResizing.hpp"

#include <algorithm>
#include <limits>

#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

VertexMapResizing::VertexMapResizing(NeighboursInterface& neighbours)
    : m_neighbours(neighbours) {}

const vector<size_t>& VertexMapResizing::operator()(size_t vertex) {
  const auto citer = m_cached_neighbours.find(vertex);
  if (citer != m_cached_neighbours.cend()) {
    return citer->second;
  }
  auto& list = m_cached_neighbours[vertex];
  list = m_neighbours(vertex);
  for (auto other_v : list) {
    m_cached_full_edges.insert(get_swap(vertex, other_v));
  }
  return list;
}

const VertexMapResizing::Result& VertexMapResizing::resize_mapping(
    VertexMapping& mapping, unsigned desired_size) {
  m_result.success = false;
  m_result.edges.clear();
  if (mapping.size() > desired_size) {
    for (auto infinite_loop_guard = 1 + mapping.size(); infinite_loop_guard > 0;
         --infinite_loop_guard) {
      const auto old_size = mapping.size();
      remove_vertex(mapping);
      const auto new_size = mapping.size();
      if (new_size <= desired_size) {
        fill_result_edges(mapping);
        m_result.success = true;
        return m_result;
      }
      if (old_size <= new_size) {
        return m_result;
      }
    }
    TKET_ASSERT(!"VertexMapResizing::resize_mapping");
  }
  TKET_ASSERT(mapping.size() <= desired_size);
  bool terminated_correctly = false;
  for (auto infinite_loop_guard = 1 + desired_size; infinite_loop_guard > 0;
       --infinite_loop_guard) {
    const auto old_size = mapping.size();
    if (old_size >= desired_size) {
      terminated_correctly = true;
      break;
    }
    add_vertex(mapping);
    const auto new_size = mapping.size();
    if (old_size == new_size) {
      // Couldn't add a vertex.
      terminated_correctly = true;
      break;
    }
    // Must have added exactly one vertex.
    TKET_ASSERT(old_size + 1 == new_size);
  }
  TKET_ASSERT(terminated_correctly);
  // It's acceptable to have too few vertices,
  // it can still be looked up in the table.
  m_result.success = true;
  fill_result_edges(mapping);
  return m_result;
}

size_t VertexMapResizing::get_edge_count(
    const VertexMapping& mapping, size_t vertex) {
  const auto& neighbours = operator()(vertex);
  return std::count_if(
      neighbours.cbegin(), neighbours.cend(),
      // Note that "neighbours" automatically will not contain "vertex" itself.
      [&mapping](size_t vertex) { return mapping.count(vertex) != 0; });
}

void VertexMapResizing::add_vertex(VertexMapping& mapping) {
  std::set<size_t> new_vertices;

  // Multipass, maybe a bit inefficient, but doesn't matter.
  // After a few calls, it's just map lookup so not so bad.
  for (const auto& existing_vertex_pair : mapping) {
    // A valid mapping should have the same source/target vertices,
    // so don't need to consider .second.
    const auto& neighbours = operator()(existing_vertex_pair.first);
    for (auto vv : neighbours) {
      if (mapping.count(vv) == 0) {
        new_vertices.insert(vv);
      }
    }
  }

  // Now find the new vertex which would add the largest number of new edges.
  size_t maximum_new_edges = 0;
  size_t best_new_vertex = std::numeric_limits<size_t>::max();

  for (auto new_v : new_vertices) {
    const auto edge_count = get_edge_count(mapping, new_v);
    if (edge_count > maximum_new_edges) {
      best_new_vertex = new_v;
      maximum_new_edges = edge_count;
    }
  }
  if (maximum_new_edges > 0) {
    mapping[best_new_vertex] = best_new_vertex;
  }
}

void VertexMapResizing::remove_vertex(VertexMapping& mapping) {
  const auto invalid_number_of_edges = std::numeric_limits<size_t>::max();

  // We want to leave as many edges as possible,
  // so we remove the minimum number.
  size_t minimum_edges_removed = invalid_number_of_edges;
  size_t best_vertex = std::numeric_limits<size_t>::max();
  for (const auto& existing_vertex_pair : mapping) {
    if (existing_vertex_pair.first != existing_vertex_pair.second) {
      // The vertex is not fixed, so we cannot remove it.
      continue;
    }
    const auto edge_count = get_edge_count(mapping, existing_vertex_pair.first);
    if (edge_count < minimum_edges_removed) {
      best_vertex = existing_vertex_pair.first;
      minimum_edges_removed = edge_count;
    }
  }
  if (minimum_edges_removed < invalid_number_of_edges) {
    TKET_ASSERT(mapping.at(best_vertex) == best_vertex);
    TKET_ASSERT(mapping.erase(best_vertex) == 1);
  }
}

void VertexMapResizing::fill_result_edges(const VertexMapping& mapping) {
  m_result.edges.clear();
  for (auto citer1 = mapping.cbegin(); citer1 != mapping.cend(); ++citer1) {
    auto citer2 = citer1;
    for (++citer2; citer2 != mapping.cend(); ++citer2) {
      const auto edge = get_swap(citer1->first, citer2->first);
      if (m_cached_full_edges.count(edge) != 0) {
        m_result.edges.push_back(edge);
      }
    }
  }
}

}  // namespace tsa_internal
}  // namespace tket
