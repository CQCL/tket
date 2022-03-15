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

#pragma once
#include <limits>
#include <set>
#include <string>
#include <utility>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

/** This is for one-off calculation at the beginning;
 * get information I(v) about a vertex v in a graph and nearby vertices,
 * and use it to decide the initial domains of the pattern vertices.
 * (E.g., if F : V(P) -> V(T) is valid, then deg(v) <= deg(F(v)),
 * together with more sophisticated checks.
 *
 * Note that it's all a balancinng act: the more we cut down the domains
 * Dom(v) for pattern vertices v, the smaller the search space,
 * and so (hopefully) the quicker the search. However, this initial
 * calculation takes time, so we must be careful not to spend too long on it.
 *
 * TODO: think about delayed initialisers (lazy evaluation):
 * only at the point where the assignment x --> y is considered for
 * the first time, apply some more expensive filters
 * to decide whether to allow it.
 *
 * TODO: does the order matter? This is based upon higher-order interactions,
 * so it's conceivable that multiple passes or reordering of calculations
 * might lead to better final results.
 *
 * TODO: think about supplemental graphs...:
 * DEFINITION: if H is any fixed graph with vertices h0, h1 in H,
 * and G is any graph, define s(G,H) ["the supplemental graph"]
 * to be a new graph with V(s(G,H)) = V(G). Define edges as follows:
 * if v0,v1 in V(G), then let (v0,v1) in E(s(G,H))
 * (i.e. there is an edge between them in the supplemental graph)
 * if there exists a non-induced subgraph isomorphism
 * f : H -> G with f(h(j)) = v(j), for j=0,1.
 *
 * THEOREM: if   F : P -> T   is a non-induced subgraph isomorphism,
 * (so that F : V(P) -> V(T)) and H is any graph with |V(H)| >= 2,
 * then F also defines a non-induced subgraph isomorphism from s(P,H)
 * to s(T,H).
 *
 * Thus, we can immediately get many NECESSARY conditions
 * for any y in Dom(v).
 */
class DomainInitialiser {
 public:
  /** Adjustable parameters to configure the calculation. */
  struct Parameters {
    /** The larger this is, the more calculation,
     * but the better the results (i.e. cutting down of possibilities,
     * for future searching).
     */
    std::size_t max_path_length;

    Parameters();
  };

  /** Apply all the different filters once at the start to reduce
   * possible_assignments as much as possible.
   * @param possible_assignments All possible pv->tv mappings. We will erase
   * some, if we can show that they are impossible.
   * @param pattern_vertices The list of all P-vertices.
   * @param pattern_neighbours_data The object to calculate information about
   * the pattern graph.
   * @param target_vertices The list of all T-vertices.
   * @param target_neighbours_data The object to calculate information about the
   * target graph.
   * @param params Extra parameters for the calculation.
   * @return False if the filters reduced any domain to empty (so that no
   * monomorphism f exists).
   */
  bool full_initialisation(
      PossibleAssignments& possible_assignments,
      const std::vector<VertexWSM>& pattern_vertices,
      const NeighboursData& pattern_neighbours_data,
      const std::vector<VertexWSM>& target_vertices,
      const NeighboursData& target_neighbours_data,
      const Parameters& params = {});

  /** After calling full_initialisation, if it returned true, then this
   * lists all pattern vertices with domain size 1, i.e. their targets
   * are uniquely determined before we've done any searching.
   * @return A sorted vector of all pattern vertices (i.e., keys in
   * "possible_assignments") whose domain (i.e., value) has size 1.
   */
  const std::vector<VertexWSM>& get_assigned_vertices() const;

 private:
  /** The most obvious initialisation: just look at degree sequences,
   * i.e. the vertex degrees of all neighbouring vertices.
   */
  bool degree_sequence_initialisation(
      PossibleAssignments& possible_assignments,
      const std::vector<VertexWSM>& pattern_vertices,
      const NeighboursData& pattern_neighbours_data,
      const std::vector<VertexWSM>& target_vertices,
      const NeighboursData& target_neighbours_data);

  /** Use the class DistanceCounts.
   * If  F : V(P) -> V(T) is valid, v in V(P),
   * and there are N vertices x with Dist(v,x) = d,
   * then there are also >= N vertices y with
   * Dist(F(v),y) <= d.
   */
  bool distance_counts_reduction(
      PossibleAssignments& possible_assignments,
      const NeighboursData& pattern_neighbours_data,
      const NeighboursData& target_neighbours_data, const Parameters& params);

  /** If any PV have domains reduced to size 1,
   * remove their targets from the other domains.
   */
  bool alldiff_reduction(PossibleAssignments& possible_assignments);

  std::vector<VertexWSM> m_work_vector;
  std::vector<VertexWSM> m_assigned_vertices;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
