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
#include <map>
#include <optional>
#include <set>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct SearchNode;
struct FixedData;

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
  /** The result of our attempted detection. */
  struct Result {
    /** The extra weight from adding all unassigned variables
     * will DEFINITELY be at least as big as this.
     *
     * (Or, if null, we're ALREADY at a nogood).
     *
     * We will break off calculation early of course
     * as soon as we exceed the bound.
     */
    std::optional<WeightWSM> extra_weight_lower_bound;

    /** If not null, gives a (PV,TV) pair where a PV->TV assignment
     * has already taken place,
     * but TV has been newly discovered to be invalid
     * (and in fact was ALWAYS invalid),
     * meaning that no p-vertex could ever be assigned to it.
     *
     * Thus, we have a very "strong" nogood,
     * giving the caller the option to backtrack
     * (possibly many levels up) to where the assignment PV->TV
     * was first made,
     * and ALSO erase TV from ALL domains.
     *
     * (As opposed to a "weak nogood", i.e. an ordinary nogood:
     * some assignment PV->TV has been made, which we discover
     * is invalid at the current position only, but might be possible
     * after backtracking.
     *
     * A "middle" nogood might be, say, the fresh discovery that PV->TV is
     * always invalid, even though TV is valid).
     */
    std::optional<std::pair<VertexWSM, VertexWSM>>
        assignment_with_invalid_t_vertex;
  };

  /** Try to find a weight nogood (i.e., an impossible current position).
   * @param fixed_data Contains edge/weight data, etc. needed for the
   * calculation.
   * @param possible_assignments Data for all unassigned vertices.
   * @param assignments Data for all assignments.
   * @param max_extra_weight The maximum extra weight we allow (the whole point
   * of this nogood detection attempt).
   * @return Information about a possible weight nogood.
   */
  Result operator()(
      const FixedData& fixed_data,
      const PossibleAssignments& possible_assignments,
      const Assignments& assignments, WeightWSM max_extra_weight) const;

 private:
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

  // Lazy initialised, with those target vertices which actually
  // are used somewhere (i.e., lie in the domain of some pattern v).
  mutable std::set<VertexWSM> m_valid_target_vertices;

  bool target_vertex_is_valid(VertexWSM tv, const FixedData& fixed_data) const;

  // Calculated only on first use, and cached;
  // if non-null, the minimum edge weight of any target edge
  // containing the target vertex tv.
  std::optional<WeightWSM> get_min_weight_for_tv(
      VertexWSM tv, const FixedData& fixed_data) const;

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
      const FixedData& fixed_data,
      const PossibleAssignments& possible_assignments) const;

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
