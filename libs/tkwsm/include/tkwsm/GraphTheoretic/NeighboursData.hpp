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

#pragma once
#include <map>
#include <optional>
#include <utility>

#include "GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** The main object used to search for neighbours of a vertex,
 * check for existing edges, and get edge weights.
 * The graph is undirected. No loops or multiple edges are allowed.
 * The input vertices MUST be {0,1,2,...,N} for some N. Will throw if not.
 * There cannot be any isolated vertices (more precisely, since edges
 * are passed in and vertices are deduced from them, it's simply not possible
 * for this class to know anything about isolated vertices).
 */
class NeighboursData {
 public:
  /** Initialise with a given graph with edge weights.
   * The graph is undirected, therefore edges (v1,v2) and (v2,v1) are
   * counted as the same. It's not necessary to pass in both pairs,
   * but if both ARE passed in, then the weights must be equal,
   * and this is checked.
   * Note that edge weights are nonnegative, but also allowed to be zero.
   * An edge weight 0 has no special meaning.
   * @param edges_and_weights An edge->weight map ("edge" just being a pair of
   * VertexWSM).
   */
  explicit NeighboursData(const GraphEdgeWeights& edges_and_weights);

  /** Return the total number of edges in the graph (remembering that the
   * graph is UNDIRECTED, so (v1,v2) and (v2,v1) are the same edge).
   * @return the number of edges in the graph.
   */
  std::size_t get_number_of_edges() const;

  /** Isolated vertices are ignored (actually, not even ignored; they are
   * simply invisible to this class, since the class only sees edges).
   * @return the number of vertices in the graph (which occurred in at least one
   * input edge).
   */
  std::size_t get_number_of_nonisolated_vertices() const;

  /** The number of neighbours of the given vertex.
   * (Returns 0 if a vertex doesn't exist. This class has no way to know
   * if a vertex is a valid isolated vertex, since isolated vertices
   * won't be listed; loop edges are not allowed).
   * @param v A vertex.
   * @return The number of neighbours of v.
   */
  std::size_t get_degree(VertexWSM v) const;

  /** Returns the edge weight for edge (v1,v2), or null if there is no edge
   * between v1, v2. (Which might mean that the vertices don't exist at all;
   * but it cannot detect this, since it does not know isolated vertices).
   * @param v1 First vertex.
   * @param v2 Second vertex.
   * @return The edge weight if it exists, or null if it does not exist.
   */
  std::optional<WeightWSM> get_edge_weight_opt(
      VertexWSM v1, VertexWSM v2) const;

  /** Return all the neighbouring vertices of v, together with
   * the edge weights, sorted by neighbouring vertex.
   * The data is stored within this class, so time O(log (number of V)).
   * @param v A vertex.
   * @return The neighbours of v and edge weights, sorted by vertex number.
   */
  const std::vector<std::pair<VertexWSM, WeightWSM>>&
  get_neighbours_and_weights(VertexWSM v) const;

  /** Get the list of vertex degrees of all neighbours of v, sorted in
   * increasing order.
   * Constructed afresh each time, so should not be called repeatedly.
   * @param v A vertex.
   * @return A list of the degrees of all neighbours of v, in increasing order.
   */
  std::vector<std::size_t> get_sorted_degree_sequence_expensive(
      VertexWSM v) const;

  /** Reconstructs the list of neighbours of v, each time.
   * Thus time O(log V + |Neighbours|).
   * One should prefer get_neighbours_and_weights whenever possible,
   * which is cheap.
   * @param v A vertex.
   * @return The recalculated neighbours of v, sorted by vertex number.
   */
  std::vector<VertexWSM> get_neighbours_expensive(VertexWSM v) const;

  /** This is "expensive" because it must be constructed afresh each time,
   * taking O(E) time. The weights are not in sorted order.
   * @return A list of all the edge weights, NOT in any particular order.
   * Duplicate values are allowed, of course; the number of weights equals the
   * number of edges.
   */
  std::vector<WeightWSM> get_weights_expensive() const;

 private:
  // Element[i] gives all neighbouring vertices and edge weights
  // for vertex i, sorted by the neighbouring vertex numerical value.
  std::vector<std::vector<std::pair<VertexWSM, WeightWSM>>>
      m_neighbours_and_weights;

  std::vector<std::pair<VertexWSM, WeightWSM>> m_empty_data;

  std::size_t m_number_of_edges;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
