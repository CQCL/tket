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

#include "TokenSwapping/DistancesInterface.hpp"
#include "TokenSwapping/NeighboursInterface.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** Contains information about a cyclic shift. Note that "moves"
 *  are not swaps; they are "half swaps". I.e., a move v1->v2
 *  means that we pretend that v2 has no token on it, and see
 *  what would happen if we moved the token on v1 to v2, ignoring
 *  whatever token is on v2.
 *  It's important to realise that moves are impossible to do by themselves,
 *  if both vertices contain tokens; it is only SEQUENCES of moves
 *  which may sometimes be converted into swaps.
 *  For example, assuming that edges v0-v1 and v1-v2 exist,
 *  the length 3 move sequence v0->v1->v2->v0
 *  may be enacted by 2 swaps (v0, v1) . (v1, v2).
 *  Notice that the edge v0-v2 does NOT have to exist. Also, this cyclic shift
 *  is still possible in 2 swaps if any 2 of the 3 edges v0-v1, v1-v2, v0-v2
 * exist.
 */
struct Cycle {
  /** By how much would L (the sum of distances from current vertex
   *  to target vertex) decrease? Can be positive or negative.
   *  It has two different interpretations:
   *  for the first, for OPEN cycles,
   *  we simply IGNORE the token on v(N), the last vertex,
   *  and store the decrease for the partial cyclic shift
   *  v0->v1->v2->v3-> ... -> v(N), AS IF there were no token on v(N).
   *
   *  For the second interpretation, once "attempt_to_close_cycles"
   *  has returned true, this switches meaning to the L-decrease for
   *  the FULL cycle, i.e. including the v(N)->v(0) decrease.
   */
  int decrease;

  /** The abstract move sequence moves each vertex to the next in the list.
   *  When the cycle is closed, the final vertex moves back to the start.
   *  [v0,v1,v2,v3,...,vN] must be a genuine path (the edges must exist),
   *  BUT the edge vN -> v0 to close the cycle does NOT have to exist.
   */
  std::vector<size_t> vertices;

  /** We need this to maintain paths without duplicate vertices.
   *  Maintaining a std::set of vertices for quick lookup would work,
   *  BUT actually "vertices" is always quite small,
   *  so just do a linear search.
   *  @param vertex A vertex
   *  @return whether that vertex already exists in "vertices".
   */
  bool contains(size_t vertex) const;
};

typedef VectorListHybrid<Cycle> Cycles;

/** Concerned only with growing cycles and closing them.
 *  For use in CyclesPartialTsa. We build up collections of cycles
 *  with information about what would happen if it were closed,
 *  i.e. the complete cycle were performed somehow,
 *  and also ensure when we grow cycles that the newly added vertex
 *  is not already present.
 *
 *  Note that longer cycles need more swaps, so our heuristic is to prefer
 *  shorter cycles, if all else is equal.
 *  (In the best possible case, if every
 *  abstract token move v(i)->v(i+1) moved one token closer to home,
 *  then the total L-decrease would be V, for V vertices, but would need
 *  V+1 swaps to perform, for a "power" (L-decrease per swap) of (V+1)/V
 *  which is actually decreasing in V).
 *  [Of course it's only a heuristic, not necessarily optimal, because
 *  doing short-term worse moves now might allow better moves
 *  in the long term - always the problem with optimisation].
 */
class CyclesGrowthManager {
 public:
  /** These control the behaviour; experimentation is needed
   *  to find the best values.
   */
  struct Options {
    size_t max_cycle_size = 6;

    /** The worst-case total number of cycles grows exponentially,
     *  e.g. the complete graph with n vertices has ~ 0.5 n^2 edges,
     *  but >> 2^n cycles.
     *
     *  We avoid exponential time/space blowup by limiting the number
     *  of cycles under consideration; any more are just discarded.
     */
    size_t max_number_of_cycles = 1000;

    /** Discard a partially built up cycle as soon as the L-decrease
     *  (decrease of total distances of vertices from their targets)
     * drops below this.
     *
     * Larger values should lead to a more "aggressive", "greedy-like"
     * algorithm, which MAY be better - more experimentation needed.
     *
     * This can even be negative, giving cycles the chance to be initially bad,
     * but later turn good.
     */
    int min_decrease_for_partial_path = 0;

    /** Similar to "min_decrease_for_partial_path", but expressed
     *  in terms of "power". Power is defined as (L-decrease)/(number of moves),
     *  which is always between -1 and +1 since each move (NOT a swap!)
     *  changes L by one of -1,0,+1.
     *  Express as a percentage to handle fractions.
     *  The partial cycle will be discarded unless BOTH criteria using
     *  min_decrease_for_partial_path AND this are satisfied.
     */
    int min_power_percentage_for_partial_path = 0;
  };

  /** Access the options, to change if desired. */
  Options& get_options();

  /** Simply returns the stored cycles. For an extra security check:
   *  unless you request otherwise (e.g., for debugging purposes),
   *  you can ONLY extract the cycles once they are
   *  converted into good candidates by "attempt_to_close_cycles".
   *  Note that some cycles may be repeated, e.g. [v0, v1, v2] and [v1, v2, v0]
   *  might both occur; further filtering is necessary.
   *
   *  Of course [v2, v1, v0] would be a totally different cycle,
   *  as the direction is reversed.
   *
   *  @param throw_if_cycles_are_not_candidates The intended use is to call this
   *    function only once candidates arise. If you call this function
   *    without having candidates, then it throws if this is set to true (the
   * default value). But for testing/debugging, it is helpful to call this
   * function just to inspect the cycles, and so this parameter should be set to
   * false.
   * @return The stored cycles.
   */
  const Cycles& get_cycles(
      bool throw_if_cycles_are_not_candidates = true) const;

  /** Start a new problem. The next function to call is
   * "attempt_to_close_cycles". Of course, swaps are just cycles with 2
   * vertices.
   *  @param vertex_mapping Where does each vertex want to move? (I.e., it's
   *    a current vertex -> target vertex map). Can be partial, i.e. not every
   * vertex has to have a token.
   *  @param distances Object to calculate distances between vertices
   *  @param neighbours Object to calculate vertices adjacent to a given vertex
   *  @return True if it found at least some good moves, false if it couldn't
   *     find ANY good moves (which must mean that all tokens are home).
   *     Recall that a move is only a "half swap".
   */
  bool reset(
      const VertexMapping& vertex_mapping, DistancesInterface& distances,
      NeighboursInterface& neighbours);

  /** For each cycle, see what would happen if we performed the full cycle
   *  (i.e., "closed the cycle").
   *  The current cycles are stored as paths [v0, v1, ..., vn], where the edges
   *  v(i) <-> v(i+1) exist, for 0 <= i < n.
   *  Even if the edge v(n)->v(0) does not exist, the cycle is POSSIBLE
   *  by "swapping along" the path [v0, v1, ..., vn]. The end result is a cyclic
   * shift. If at least one cycle could be closed to create a viable candidate
   *  (giving a net decrease in L), return true and delete all cycles which are
   *  NOT candidates, and also fill in the L-decrease values for the CLOSED
   * cycle. If NO cycle closures give a good result, do nothing and return
   * false.
   *  @param vertex_mapping The desired (source vertex->target vertex) mapping,
   *    for the current locations of tokens on vertices.
   *  @param distances Object to calculate distances, used to calculate
   * L-decreases.
   *  @return True if at least one good closed cycle exists (i.e., giving net
   * strict decrease of L). If so, all non-good cycles are deleted. If no good
   * closed cycle exists, do nothing and return false.
   */
  bool attempt_to_close_cycles(
      const VertexMapping& vertex_mapping, DistancesInterface& distances);

  /** Record what happens when we try to GROW cycles (i.e., increase the length
   * of each stored cycle by one, discarding all those which could not grow).
   */
  struct GrowthResult {
    /** If TRUE, there are no more cycles to consider; finish. */
    bool empty = false;

    /** If we're already at the length limit, delete all cycles.
     *  (There is no further use for them, so this is safest).
     *  However, this is not the only possible way for all cycles to be deleted.
     *  There might not be any other vertices in the graph to add;
     *  or they might all be bad cycles (i.e., not decreasing L by enough to
     * keep them).
     */
    bool hit_cycle_length_limit = false;
  };

  /** For each existing cycle, try all possible ways to extend it
   *  by one step from the last vertex.
   *  Keep all new cycles generated in this way with a good L decrease,
   *  and discard all others (including the original cycle).
   *  Thus, all cycles should have the same number of vertices, increasing
   *  by one each time this function is called (unless they are all deleted).
   *  @param vertex_mapping The current desired (source vertex -> target vertex)
   * mapping.
   *  @param distances Object to calculate distances, used to calculate
   * L-decreases.
   *  @param neighbours Object to calculate adjacent vertices to a given vertex.
   *  @return What happened when we tried to grow the cycles.
   */
  GrowthResult attempt_to_grow(
      const VertexMapping& vertex_mapping, DistancesInterface& distances,
      NeighboursInterface& neighbours);

 private:
  Cycles m_cycles;
  Options m_options;
  bool m_cycles_are_candidates = false;
};

}  // namespace tsa_internal
}  // namespace tket
