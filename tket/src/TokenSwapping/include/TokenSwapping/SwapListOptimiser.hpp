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

#include "DynamicTokenTracker.hpp"

namespace tket {
namespace tsa_internal {

/** Can be reused, faster than constructing a new object.
 *  (This is intended: the final full algorithm may involve much
 *  random chopping and changing, simulated annealing, etc. and so
 *  call this many times).
 *  This is about directly optimising a list of swaps,
 *  knowing nothing about target vertices and tokens.
 *  Each optimisation pass may reorder the swaps and erase some,
 *  but always such that the resultant start-to-end vertex permutation
 *  is unchanged.
 *  (Thus, unmentioned swaps can never be added, because this class
 *  has no way of knowing if such swaps are possible).
 *  Thus this can convert bad solutions into better ones, but if our problem
 *  has some empty tokens, i.e. tokens with no target, which can end up
 *  anywhere, then algorithms will need that data to get good solutions.
 */
class SwapListOptimiser {
 public:
  using ID = SwapList::ID;

  /** The most trivial O(1) optimisation: the new swap is added to the back,
   *  unless it equals the current last swap, in which case they cancel
   *  each other and are simply removed. However, all other
   *  optimisation passes also do this as a byproduct.
   *  @param list The swaps.
   *  @param swap The single swap to be pushed back to the list (but
   *    adjacent equal swaps will be erased).
   */
  static void push_back(SwapList& list, const Swap& swap);

  /** The slowest but most accurate end-to-end optimisation pass.
   *  It will include other optimisation passes, so don't bother
   *  calling them also. This is possibly O(N^3.log N)
   *  in the worst case (but a clever proof might reduce this),
   *  but hopefully much faster in practice.
   *  @param list The swaps to be optimised.
   */
  void full_optimise(SwapList& list);

  /** Call when you also know the tokens. The slowest but hopefully best
   *  optimisation pass, which also removes empty swaps (swaps in which
   *  neither vertex has a token, so that they have no effect).
   *  @param list The swaps to be optimised.
   *  @param vertex_mapping The desired source->token mapping.
   */
  void full_optimise(SwapList& list, const VertexMapping& vertex_mapping);

  // Most optimisation passes below are O(N^2.log N) in the worst case,
  // but in practice will hopefully be a lot faster.
  // It's hard to compare passes;
  // for any two passes A, B there are probably examples where
  // pass A is better than B, but others where B is better than A.
  // Also, passes are not commutative; reordering the passes
  // can give different results! Experimentation needed.

  /** Do not move any swaps, unless it cancels with a previous copy of itself
   *  (in which case, delete both). The fastest pass.
   *  @param list The swaps to be optimised.
   */
  void optimise_pass_with_zero_travel(SwapList& list);

  /** Starting from the front and working along, every swap is moved
   *  as far towards the front as possible, until it hits a non-disjoint
   *  swap (so that it cannot pass through it; it doesn't commute),
   *  or an identical copy of itself (so that they cancel each other).
   *  The overall reduction should be the same as
   *  optimise_pass_with_zero_travel, which is cheaper
   *  (because swaps do not move), but interacting swaps should
   *  cluster together, which may be useful for certain algorithms.
   *  @param list The swaps to be optimised.
   */
  void optimise_pass_with_frontward_travel(SwapList& list);

  /** Erase two swaps if they do the same TOKEN swap (which means that
   *  they can be removed). Knows nothing about the problem-specific tokens,
   *  instead this creates artificial tokens.
   *  This is slower than optimise_pass_with_zero_travel and
   *  optimise_pass_with_frontward_travel, but is strictly more powerful
   *  (any reduction by those passes will also occur with this pass,
   *  but some additional reductions are possible with this pass. E.g.,
   *  this pass reduces (01)(12)(01)(12)(01)(12), the cube of a 3-cycle,
   *  to zero swaps, which the other passes cannot, since (01) and (12)
   *  are not disjoint and hence cannot pass through each other).
   *  However, NOTE that this pass can introduce EMPTY swaps, w.r.t.
   *  the problem-specific tokens, so further passes to remove
   *  problem-specific empty token swaps are necessary
   *  to get the full reduction.
   *  @param list The swaps to be optimised.
   */
  void optimise_pass_with_token_tracking(SwapList& list);

  /** O(N log N): simply discard any swap between two empty tokens.
   *  (Recall that optimise_pass_with_token_tracking does NOT know
   *  about these specific tokens, it creates internal artificial ones
   *  just for the pass. That pass can be much slower than this pass,
   *  but also can make some reductions which this pass cannot).
   *  @param list The swaps to be optimised.
   *  @param vertex_mapping The desired source->target mapping (so that
   *    we can determine which vertices have tokens on them).
   */
  void optimise_pass_remove_empty_swaps(
      SwapList& list, VertexMapping vertex_mapping);

 private:
  std::map<Swap, size_t> m_data;

  DynamicTokenTracker m_token_tracker;

  /** What would happen if you tried to move the swap towards the front?
   *  Doesn't actually move the swap, just returns the ID of the first
   *  blocking swap it hits (or null if there is none and it could move
   *  all the way to the front), UNLESS it actually hits another copy of itself,
   *  in which case it DOES erase (and the caller can tell by checking the
   * size).
   *  @param list The swaps to be optimised.
   *  @param id The ID of the current swap which we might move frontwards.
   *  @return The ID of the previous reachable distinct swap which is
   *    non-disjoint (has a vertex in common, so doesn't commute),
   *    or empty if none exists.
   */
  static std::optional<ID> get_id_of_previous_blocker(
      SwapList& list, SwapID id);

  /** Actually move the swap as far towards the front as possible until
   *  blocked, erasing it if it cancelled with another copy of itself.
   *  @param list The swaps to be optimised.
   *  @param id The ID of the current swap to move frontwards.
   *  @return true if the swap cancelled with a copy of itself.
   */
  static bool move_swap_towards_front(SwapList& list, SwapID id);

  /** The same as optimise_pass_with_token_tracking,
   *  but without calling "clear" OR "reset" on m_token_tracker first).
   *  @param list The swaps to be optimised.
   */
  void optimise_pass_with_token_tracking_without_clearing_tracker(
      SwapList& list);
};

}  // namespace tsa_internal
}  // namespace tket
