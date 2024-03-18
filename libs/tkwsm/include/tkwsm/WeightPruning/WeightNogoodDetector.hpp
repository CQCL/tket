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
#include <optional>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;
class NeighboursData;

/** Can we get a good guaranteed lower bound on the extra weight
 * to be added from all unassigned pattern vertices,
 * HOWEVER it's done?
 * If so, and if this extra weight would take the total weight
 * above the constraint value (i.e., the imposed upper bound
 * on the weight), then we can prune immediately (backtrack).
 *
 * We want something fast enough not to slow the
 * calculation down too much, but accurate enough (i.e. giving
 * large lower bounds on the extra weight) to be useful,
 * i.e. actually detect weight nogoods sometimes.
 * A delicate balancing act!
 *
 * The ACTUAL optimal lower bound can only be found, presumably, by solving
 * a complete WSM problem. More accurate estimates can be obtained
 * by weighted bipartite matching, etc. etc., but take longer.
 */
class WeightNogoodDetector {
 public:
  /** The detector needs to know all target vertices which are possible.
   * But note that this
   * need not be done at the START of the search; it can be delayed
   * until needed, when search progress may have cut down some
   * possibilities, and hence reduced the domains slightly,
   * giving better estimates.
   */
  WeightNogoodDetector(
      const NeighboursData& pattern_neighbours_data,
      const NeighboursData& target_neighbours_data,
      std::set<VertexWSM> initial_used_target_vertices,
      std::set<VertexWSM>& invalid_target_vertices);

  /** Return a guaranteed lower bound for the scalar product increase
   * if the current partial solution could be extended to a complete
   * valid solution, or null if in fact we can prove that it's impossible
   * (EITHER for graph theoretic reasons - maybe the graphs simply aren't
   * compatible, regardless of weights - OR because it would exceed
   * the allowed maximum).
   * @param accessor An object to retrieve information about domains.
   * @param max_extra_scalar_product The maximum extra weight we allow (the
   * whole point of this nogood detection attempt).
   * @return Information about the min possible scalar product increase (null if
   * it's a nogood)
   */
  std::optional<WeightWSM> get_extra_scalar_product_lower_bound(
      const DomainsAccessor& accessor, WeightWSM max_extra_scalar_product);

  /** How many target vertices are still valid (i.e., could possibly
   * be mapped to by some PV)?
   * @return the size of m_valid_target_vertices
   */
  std::size_t get_number_of_possible_tv() const;

 private:
  const NeighboursData& m_pattern_neighbours_data;
  const NeighboursData& m_target_neighbours_data;

  mutable std::set<VertexWSM> m_valid_target_vertices;

  // As soon as a TV is newly discovered to be invalid,
  // meaning that no p-vertex could ever be assigned to it
  // (quite a rare event), store it here.
  //
  // This of course is even stronger than a normal impossible
  // assignment PV->TV.
  //
  // Moreover (in case the caller ever tries parallel searches, etc.)
  // the TV was ALWAYS invalid, so can be erased from ALL data.
  std::set<VertexWSM>& m_invalid_target_vertices;

  // There are many mutable data members.
  // Even though this class is "logically" const,
  // it's not "physically" const because of lazy evaluation,
  // as well as work data to avoid memory reallocation.

  // KEY: a target vertex TV.
  // VALUE: the smallest t-edge weight which can arise
  //    from an edge containing TV.
  // Excludes target edges, if any, which can never arise
  // (i.e. containing at least one vertex which no pattern vertex
  // can ever map to).
  //
  // Mutable because it's lazily initialised.
  mutable std::map<VertexWSM, WeightWSM> m_minimum_t_weights_from_tv;

  // Calculated only on first use, and cached;
  // if non-null, the minimum edge weight of any target edge
  // containing the target vertex tv.
  std::optional<WeightWSM> get_min_weight_for_tv(VertexWSM tv) const;

  // This really is mutable, to avoid reallocation.
  // In each pair: the FIRST is an UNASSIGNED pattern vertex pv.
  // The second is a guaranteed lower bound for the t-weight which will be
  // matched to any p-edge containing pv. (The p-edge therefore must necessarily
  // be unassigned; but BEWARE: the other vertex could be assigned already).
  //
  // It just so happens, from the way it's created, that it's ALWAYS sorted!
  // Thus, we can use a sorted vector instead of a map.
  mutable std::vector<std::pair<VertexWSM, WeightWSM>>
      m_t_weight_lower_bounds_for_p_edges_containing_pv;

  // Fills m_t_weight_lower_bounds_for_p_edges_containing_pv,
  // relevant for this current search node only.
  // Returns false if it found we're already at a nogood.
  bool fill_t_weight_lower_bounds_for_p_edges_containing_pv(
      const DomainsAccessor& accessor) const;

  // Simply looks up and returns the weight in
  // m_t_weight_lower_bounds_for_p_edges_containing_pv.
  // So, out of all t-edges which could possibly contain
  // F(pv), for any valid mapping F from this current node,
  // this is the lowest weight.
  // We're GUARANTEED that it exists, because
  // m_t_weight_lower_bounds_for_p_edges_containing_pv was already filled
  // with EVERY pv, except when there were no valid t-edges,
  // in which case it was a nogood and we should already have returned.
  WeightWSM get_t_weight_lower_bound(VertexWSM pv) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
