// Copyright Quantinuum
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

#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <stdexcept>
#include <tkassert/Assert.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Intended for use with Architecture and QubitGraph, which are similar
 * but different types. Calculate new VertexWSM vertex labels.
 */
template <class VertexType, class GraphType>
class RelabelledGraphWSM {
 public:
  explicit RelabelledGraphWSM(const GraphType& graph) {
    // Get the new vertex integer labels.
    m_original_vertices.reserve(boost::num_vertices(graph));
    {
      const auto vertex_iter_pair = boost::vertices(graph);
      for (auto citer = vertex_iter_pair.first;
           citer != vertex_iter_pair.second; ++citer) {
        // Add the vertex to the map keys.
        m_old_to_new_vertex_map[graph[*citer]];
      }
      // Now add the new vertex labels.
      for (auto& entry : m_old_to_new_vertex_map) {
        entry.second = m_original_vertices.size();
        m_original_vertices.emplace_back(entry.first);
      }
    }
    // Get the newly labelled edges.
    const auto edge_iter_pair = boost::edges(graph);
    for (auto citer = edge_iter_pair.first; citer != edge_iter_pair.second;
         ++citer) {
      auto edge = get_edge(
          get_relabelled_vertex(graph[boost::source(*citer, graph)]),
          get_relabelled_vertex(graph[boost::target(*citer, graph)]));
      m_relabelled_edges_and_weights[edge] = unsigned(graph[*citer].weight);
    }
    // Now, classify the vertices into isolated and nonisolated categories
    for (const auto& relabelled_edge_entry : m_relabelled_edges_and_weights) {
      const auto& relabelled_edge = relabelled_edge_entry.first;
      m_relabelled_nonisolated_vertices.insert(relabelled_edge.first);
      m_relabelled_nonisolated_vertices.insert(relabelled_edge.second);
    }
    for (const auto& entry : m_old_to_new_vertex_map) {
      if (m_relabelled_nonisolated_vertices.count(entry.second) == 0) {
        m_relabelled_isolated_vertices.insert(entry.second);
      }
    }
    TKET_ASSERT(
        m_old_to_new_vertex_map.size() ==
        m_relabelled_isolated_vertices.size() +
            m_relabelled_nonisolated_vertices.size());
  }

  const GraphEdgeWeights& get_relabelled_edges_and_weights() const {
    return m_relabelled_edges_and_weights;
  }

  const std::set<VertexWSM>& get_relabelled_isolated_vertices() const {
    return m_relabelled_isolated_vertices;
  }

  const std::set<VertexWSM>& get_relabelled_nonisolated_vertices() const {
    return m_relabelled_nonisolated_vertices;
  }

  // element [i] is the vertex which has been relabelled i.
  const std::vector<VertexType>& get_original_vertices() const {
    return m_original_vertices;
  }

  VertexWSM get_relabelled_vertex(const VertexType& original_vertex) const {
    const auto v_opt =
        get_optional_value(m_old_to_new_vertex_map, original_vertex);
    if (!v_opt) {
      throw std::runtime_error("Original vertex has no new label");
    }
    return v_opt.value();
  }

 private:
  std::vector<VertexType> m_original_vertices;
  std::map<VertexType, VertexWSM> m_old_to_new_vertex_map;
  std::set<VertexWSM> m_relabelled_isolated_vertices;
  std::set<VertexWSM> m_relabelled_nonisolated_vertices;

  // All edge weights will be 1, since we're only considering
  // unweighted problems.
  GraphEdgeWeights m_relabelled_edges_and_weights;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
