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

#pragma once
#include <set>

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NearNeighboursData;
class NeighboursData;

/** For constructing the initial domains Dom(pv) for each pattern vertex pv.
 * There is always a balancing act involved: to make the searching quicker,
 * we want smaller domains (both for pruning, and to reduce the amount
 * of data copying). However, making smaller domains takes more computation.
 */
class DomainInitialiser {
 public:
  /** Element[pv] is Domain(pv), i.e. the set of all target vertices
   * which pv could map to.
   */
  typedef std::vector<boost::dynamic_bitset<>> InitialDomains;

  /** Fill initial_domains, using the different filters to try to make domains
   * small; the whole purpose of this class. The necessary data about the
   * pattern and target graphs is already contained within the NeighboursData
   * objects.
   * @param initial_domains All possible pv->tv mappings, returned by reference.
   * We will erase some, if we can show that they are impossible.
   * @param pattern_neighbours_data The object to calculate information about
   * the pattern graph (immediate neighbours of a vertex).
   * @param pattern_near_neighbours_data An object to calculate nearby vertices
   * in the pattern graph (and cache calculated results).
   * @param target_neighbours_data The object to calculate information about the
   * target graph.
   * @param target_near_neighbours_data The object for calculating nearby
   * vertices in the target graph (not just immediate neighbours).
   * @param max_path_length Extra parameter for the calculation. Larger values
   * need more computation, but give smaller domains.
   * @return False if the filters reduced any domain to empty (so that no
   * monomorphism f exists).
   */
  static bool full_initialisation(
      InitialDomains& initial_domains,
      const NeighboursData& pattern_neighbours_data,
      NearNeighboursData& pattern_near_neighbours_data,
      const NeighboursData& target_neighbours_data,
      NearNeighboursData& target_near_neighbours_data,
      unsigned max_path_length);

 private:
  /** The most obvious initialisation: just look at degree sequences,
   * i.e. the vertex degrees of all neighbouring vertices.
   * Begin here by filling initial_domains.
   */
  static bool degree_sequence_initialisation(
      InitialDomains& initial_domains,
      const NeighboursData& pattern_neighbours_data,
      const NeighboursData& target_neighbours_data);

  /** If  F : V(P) -> V(T) is valid, v in V(P),
   * and there are N vertices x with Dist(v,x) = d,
   * then there are also >= N vertices y with
   * Dist(F(v),y) <= d.
   */
  static bool distance_counts_reduction(
      InitialDomains& initial_domains,
      NearNeighboursData& pattern_near_neighbours_data,
      NearNeighboursData& target_near_neighbours_data,
      unsigned max_path_length);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
