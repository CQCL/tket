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

#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NeighboursData::NeighboursData() : m_number_of_edges(0) {}

NeighboursData::NeighboursData(const GraphEdgeWeights& edges_and_weights) {
  initialise(edges_and_weights);
}

std::size_t NeighboursData::get_number_of_edges() const {
  return m_number_of_edges;
}

std::vector<VertexWSM> NeighboursData::get_neighbours_expensive(
    VertexWSM v) const {
  std::vector<VertexWSM> result;
  const auto& neighbours_and_weights = get_neighbours_and_weights(v);
  result.reserve(neighbours_and_weights.size());
  for (const auto& entry : neighbours_and_weights) {
    result.push_back(entry.first);
  }
  return result;
}

std::size_t NeighboursData::get_number_of_nonisolated_vertices() const {
  return m_neighbours_and_weights_map.size();
}

std::vector<VertexWSM> NeighboursData::get_nonisolated_vertices_expensive()
    const {
  std::vector<VertexWSM> nonisolated_vertices;
  nonisolated_vertices.reserve(m_neighbours_and_weights_map.size());
  for (auto& entry : m_neighbours_and_weights_map) {
    nonisolated_vertices.push_back(entry.first);
  }
  return nonisolated_vertices;
}

void NeighboursData::initialise(const GraphEdgeWeights& edges_and_weights) {
  m_number_of_edges = edges_and_weights.size();
  m_neighbours_and_weights_map.clear();

  for (const auto& entry : edges_and_weights) {
    const auto& v1 = entry.first.first;
    const auto& v2 = entry.first.second;
    if (v1 == v2) {
      throw std::runtime_error("Loop found in graph; not allowed");
    }
    m_neighbours_and_weights_map[v1].emplace_back(v2, entry.second);
    m_neighbours_and_weights_map[v2].emplace_back(v1, entry.second);
  }
  for (auto& entry : m_neighbours_and_weights_map) {
    auto& neigh_data = entry.second;
    // Automatically sort by vertex first; lexicographic.
    std::sort(neigh_data.begin(), neigh_data.end());
    TKET_ASSERT(is_sorted_and_unique(neigh_data));
  }
}

std::optional<WeightWSM> NeighboursData::get_edge_weight_opt(
    VertexWSM v1, VertexWSM v2) const {
  const auto v1_citer = m_neighbours_and_weights_map.find(v1);
  if (v1_citer == m_neighbours_and_weights_map.cend()) {
    return {};
  }
  const auto& v1_data = v1_citer->second;
  // (x,0) = (x,min) <= (x,w) <= (x+1, y) in lexicographic order
  std::pair<VertexWSM, WeightWSM> key;
  key.first = v2;
  key.second = 0;
  const auto v2_citer = std::lower_bound(v1_data.cbegin(), v1_data.cend(), key);
  if (v2_citer != v1_data.cend() && v2_citer->first == v2) {
    return v2_citer->second;
  }
  return {};
}

std::size_t NeighboursData::get_degree(VertexWSM v) const {
  const auto v_citer = m_neighbours_and_weights_map.find(v);
  if (v_citer == m_neighbours_and_weights_map.cend()) {
    return 0;
  }
  return v_citer->second.size();
}

std::vector<std::size_t> NeighboursData::get_sorted_degree_sequence_expensive(
    VertexWSM v) const {
  std::vector<std::size_t> result;
  const auto v_citer = m_neighbours_and_weights_map.find(v);
  if (v_citer == m_neighbours_and_weights_map.cend()) {
    return result;
  }
  for (const auto& entry : v_citer->second) {
    result.push_back(get_degree(entry.first));
  }
  std::sort(result.begin(), result.end());
  return result;
}

const std::vector<std::pair<VertexWSM, WeightWSM>>&
NeighboursData::get_neighbours_and_weights(VertexWSM v) const {
  const auto v_citer = m_neighbours_and_weights_map.find(v);
  if (v_citer == m_neighbours_and_weights_map.cend()) {
    return m_empty_data;
  }
  return v_citer->second;
}

const NeighboursData::NeighboursMap& NeighboursData::get_map() const {
  return m_neighbours_and_weights_map;
}

bool NeighboursData::is_adjacent_to_assigned_pv(
    VertexWSM pv, const Assignments& assignments) const {
  // Crude: could do fancy back-and-forth iterator trick.
  if (assignments.empty()) {
    return false;
  }
  const auto& neighbours_and_weights = get_neighbours_and_weights(pv);
  for (const auto& entry : neighbours_and_weights) {
    if (assignments.count(entry.first) != 0) {
      return true;
    }
  }
  return false;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
