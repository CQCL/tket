// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/NeighboursData.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

std::size_t NeighboursData::get_number_of_edges() const {
  return m_number_of_edges;
}

std::vector<VertexWSM> NeighboursData::get_neighbours_expensive(
    VertexWSM v) const {
  std::vector<VertexWSM> result;
  const std::vector<std::pair<VertexWSM, WeightWSM>>& neighbours_and_weights =
      get_neighbours_and_weights(v);
  result.reserve(neighbours_and_weights.size());
  for (const std::pair<VertexWSM, WeightWSM>& entry : neighbours_and_weights) {
    result.push_back(entry.first);
  }
  return result;
}

std::size_t NeighboursData::get_number_of_nonisolated_vertices() const {
  return m_neighbours_and_weights.size();
}

NeighboursData::NeighboursData(const GraphEdgeWeights& edges_and_weights) {
  // Only the edges (v1,v2) with v1<v2.
  std::set<std::pair<VertexWSM, VertexWSM>> ordered_edges_seen;
  std::set<VertexWSM> vertices_seen;
  for (const std::pair<const EdgeWSM, WeightWSM>& entry : edges_and_weights) {
    const VertexWSM& v1 = entry.first.first;
    const VertexWSM& v2 = entry.first.second;
    vertices_seen.insert(v1);
    vertices_seen.insert(v2);
    if (v1 == v2) {
      throw std::runtime_error("Loop found in graph; not allowed");
    }
    ordered_edges_seen.insert(get_edge(v1, v2));
  }
  if (vertices_seen.empty()) {
    throw std::runtime_error("No edges passed to NeighboursData");
  }
  if (*vertices_seen.cbegin() != 0 ||
      *vertices_seen.crbegin() != vertices_seen.size() - 1) {
    throw std::runtime_error("Vertices should be [0,1,2,...,N].");
  }
  for (const EdgeWSM& edge : ordered_edges_seen) {
    const std::optional<WeightWSM> weight_opt =
        get_optional_value(edges_and_weights, edge);
    if (weight_opt) {
      const std::optional<WeightWSM> other_weight_opt = get_optional_value(
          edges_and_weights, std::make_pair(edge.second, edge.first));
      if (other_weight_opt && weight_opt.value() != other_weight_opt.value()) {
        throw std::runtime_error("Edge weights mismatch");
      }
    }
  }
  m_number_of_edges = ordered_edges_seen.size();
  m_neighbours_and_weights.resize(vertices_seen.size());

  for (const EdgeWSM& edge : ordered_edges_seen) {
    const VertexWSM& v1 = edge.first;
    const VertexWSM& v2 = edge.second;
    TKET_ASSERT(v1 < v2);
    TKET_ASSERT(v2 < vertices_seen.size());
    WeightWSM weight;
    const std::optional<WeightWSM> weight_opt =
        get_optional_value(edges_and_weights, edge);
    if (weight_opt) {
      weight = weight_opt.value();
    } else {
      weight = edges_and_weights.at(std::make_pair(v2, v1));
    }
    m_neighbours_and_weights[v1].emplace_back(v2, weight);
    m_neighbours_and_weights[v2].emplace_back(v1, weight);
  }
  for (std::vector<std::pair<VertexWSM, WeightWSM>>& neigh_data :
       m_neighbours_and_weights) {
    // Automatically sort by vertex first; lexicographic.
    std::sort(neigh_data.begin(), neigh_data.end());
    TKET_ASSERT(is_sorted_and_unique(neigh_data));
  }
}

std::optional<WeightWSM> NeighboursData::get_edge_weight_opt(
    VertexWSM v1, VertexWSM v2) const {
  if (v1 >= m_neighbours_and_weights.size()) {
    return {};
  }
  const std::vector<std::pair<VertexWSM, WeightWSM>>& v1_data =
      m_neighbours_and_weights[v1];
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
  if (v >= m_neighbours_and_weights.size()) {
    return 0;
  }
  return m_neighbours_and_weights[v].size();
}

std::vector<std::size_t> NeighboursData::get_sorted_degree_sequence_expensive(
    VertexWSM v) const {
  std::vector<std::size_t> result;
  if (v >= m_neighbours_and_weights.size()) {
    return result;
  }
  for (const std::pair<VertexWSM, WeightWSM>& entry :
       m_neighbours_and_weights[v]) {
    result.push_back(get_degree(entry.first));
  }
  std::sort(result.begin(), result.end());
  return result;
}

const std::vector<std::pair<VertexWSM, WeightWSM>>&
NeighboursData::get_neighbours_and_weights(VertexWSM v) const {
  if (v >= m_neighbours_and_weights.size()) {
    return m_empty_data;
  }
  return m_neighbours_and_weights[v];
}

std::vector<WeightWSM> NeighboursData::get_weights_expensive() const {
  std::vector<WeightWSM> weights;
  weights.reserve(m_number_of_edges);
  for (unsigned v1 = 0; v1 < m_neighbours_and_weights.size(); ++v1) {
    // Every edge is implicitly stored twice, for (v1, v2) and (v2, v1).
    // To avoid duplicates, only write the weight when v1>v2.
    // The neighbour edges are stored with increasing v, as always.
    for (const std::pair<VertexWSM, WeightWSM>& inner_entry :
         m_neighbours_and_weights[v1]) {
      const VertexWSM& v2 = inner_entry.first;
      if (v2 > v1) {
        break;
      }
      weights.emplace_back(inner_entry.second);
    }
  }
  TKET_ASSERT(weights.size() == m_number_of_edges);
  return weights;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
