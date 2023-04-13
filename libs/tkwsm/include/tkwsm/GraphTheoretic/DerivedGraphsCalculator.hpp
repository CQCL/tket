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
#include "DerivedGraphStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

/** If    f : V(P) -> V(T)   is a valid subgraph monomorphism,
 * it is easy to show that f is also a monomorphism between P', T'
 * for various derived graphs (e.g., called "supplemental graphs" in
 * the paper "Parallel, Backjumping Subgraph Isomorphism Algorithm using
 * Supplemental Graphs" by Ciaran McCreesh and Patrick Prosser).
 * In particular, for k=2,3,4,... let D(k)(G) be a weighted graph derived
 * from G by setting V(D(k)(G)) = V(G), and (u,v) is an edge in D(k)(G)
 * if and only if there is a path of length k joining u to v.
 * Furthermore, the edge (u,v) has weight n, where n is the number of
 * distinct paths of length k.
 * Thus, if f is a valid mapping, then it remains valid mapping D(k)P
 * to D(k)T, with also Weight(f(u,v)) >= Weight(u,v).
 * We just say D2, D3 more briefly, instead of D(2)(P), D(2)(T), etc.
 * The edges of D2, D3 are calculated lazily since, particularly for T,
 * there may be many unused vertices and so we don't calculate the edges
 * until they're needed.
 */
class DerivedGraphsCalculator {
 public:
  /** Upon demand, calculates the neighbours in D2, D3 of the given vertex.
   * @param ndata The NeighboursData object for the original graph.
   * @param v The vertex in the original graph.
   * @param triangle_count The triangle count in the original graph.
   * @param depth_2_neighbours_and_counts Data for the neighbours of v in the
   * derived graph D2. Sorted by vertex.
   * @param depth_3_neighbours_and_counts Data for the neighbours of v in the
   * derived graph D3. Sorted by vertex.
   */
  void fill(
      const NeighboursData& ndata, VertexWSM v,
      DerivedGraphStructs::Count& triangle_count,
      DerivedGraphStructs::NeighboursAndCounts& depth_2_neighbours_and_counts,
      DerivedGraphStructs::NeighboursAndCounts& depth_3_neighbours_and_counts);

 private:
  // KEY: a vertex v2
  // VALUE: sorted vector of all v1 such that v0--v1--v2 is a path.
  std::map<VertexWSM, std::vector<VertexWSM>>
      m_mid_vertices_for_length_two_paths;

  // KEY: vertex v3
  // VALUE: the number of distinct paths v--v1--v2--v3.
  std::map<VertexWSM, DerivedGraphStructs::Count>
      m_depth_3_neighbours_and_counts_map;

  // Call the following functions, in order, to build up the data.

  void fill_mid_vertices_for_length_two_paths(
      const NeighboursData& ndata, VertexWSM v);

  void fill_d2_neighbours_and_counts(
      DerivedGraphStructs::NeighboursAndCounts& depth_2_neighbours_and_counts);

  void fill_d3_neighbours_and_counts_map(const NeighboursData& ndata);

  void fill_remaining_d3_data(
      VertexWSM v, DerivedGraphStructs::Count& triangle_count,
      DerivedGraphStructs::NeighboursAndCounts& depth_3_neighbours_and_counts);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
