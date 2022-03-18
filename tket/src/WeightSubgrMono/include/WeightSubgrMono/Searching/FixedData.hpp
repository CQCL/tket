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

#include "../DomainInitialising/DomainInitialiser.hpp"
#include "../GraphTheoretic/NeighboursData.hpp"
#include "../Reducing/AllDiffPropagator.hpp"
#include "../Reducing/HallSetReducer.hpp"
#include "../WeightPruning/WeightNogoodDetector.hpp"
#include "SearchNode.hpp"
#include "WeightUpdater.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {


/** Store all data which is unchanged throughout the backtracking
 * conveniently in one place.
 * TODO: however, we may introduce some lazy initialisation; i.e.,
 * allowing data to change over time, but ONLY with results which were
 * always true from the start, just maybe unknown/uncomputed at that time.
 * (E.g., maybe a pattern vertex PV has a nearby vertex at a certain distance
 * with high degree, but a target vertex TV initially in Dom(PV) has all
 * nearby vertices of too low a degree; but this computation was not carried
 * out until later. It means that we could have removed TV from Dom(PV)
 * right at the start, if we had known.
 *
 * A more extreme possibility is: we may find that a target vertex TV
 * is removed from ALL domains; in which case, we can ERASE TV completely
 * from the target graph. This would also change the neighbours data,
 * possibly causing new reductions, etc.
 * As always, the domains would be pruned and reduced as much as possible
 * if we updated all this information; but we'd have to be careful that
 * the extra computation time from doing all this is not
 * too much to be worthwhile).
 */
struct FixedData {
  /** Data about the pattern graph. */
  NeighboursData pattern_neighbours_data;

  /** Data about the target graph. */
  NeighboursData target_neighbours_data;

  /** Does the target graph contain every possible edge?
   * Obviously, if the target graph is complete
   * (at least, for all the target vertices mentioned in the edge weights)
   * then no graph theoretic vertex filtering is possible;
   * EVERY tv lies in every PV domain.
   */
  bool target_is_complete;

  /** If all pattern weights are equal, AND all target weights are equal,
   * it's effectively unweighted (a pure subgraph isomorphism problem).
   */
  bool problem_is_unweighted;

  /** Contains the initial lists of possible pv->tv mappings, etc.
   * As each new search starts - if we are interleaving searches -
   * this should be copied to the first node, but also reduced
   * (e.g. if other search branches have uncovered nogoods, etc.
   * then the domains will reduce).
   */
  SearchNode initial_node;

  /** Simply store the sum of all p edge weights. */
  WeightWSM total_p_edge_weights;

  // Store reduction/updating/detector objects which don't need
  // extra state data attached to a specific node.
  // Really just for convenience, as they have to be stored SOMEWHERE
  // (they contain internal data structures etc. which it would be
  // wasteful to deallocate and reallocate each time).

  /** Detects subsets of p-vertices whose combined domain is small.
   * This may be viewed as a generalisation of an AllDiffPropagator.
   */
  HallSetReducer hall_set_reducer;

  /** Removes newly assigned TV from other domains. */
  AllDiffPropagator alldiff_propagator;

  /** Calculates the total weight (scalar product) so far
   * of all assignments; this is done as each new pattern edge
   * becomes assigned.
   */
  WeightUpdater weight_updater;

  /** Unlike weight_updater, which finds the actual weight so far,
   * this tries to predict how much FUTURE weight might be added,
   * at a minimum; if a guaranteed lower bound is found which would
   * take the total weight over a given maximum, then we can immediately
   * backtrack.
   */
  WeightNogoodDetector weight_nogood_detector;

  /** Returns FALSE if impossible; fills all the above data. */
  bool initialise(
      const GraphEdgeWeights& p_data, const GraphEdgeWeights& t_data,
      DomainInitialiser::Parameters parameters = {});
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
