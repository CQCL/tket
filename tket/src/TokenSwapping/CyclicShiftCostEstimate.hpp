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

#include "DistancesInterface.hpp"

namespace tket {
namespace tsa_internal {

/** Used in the TrivialTSA class (NOT in CyclesPartialTsa!)
 *  Given a desired abstract cyclic shift on [v0, v1, v2, ..., vn],
 *  i.e. abstract moves v(0)->v(1)->v(2)-> ... ->v(n)->v(0),
 *  [meaning that v(i), v(i+1) need not actually be adjacent in the graph,
 *  so we must decide how to represent the desired moves as actual swaps],
 *  there are n+1 possible obvious ways to enact it
 *  (and of course, maybe some "nonobvious" ways.
 *  Finding a good way is, of course, a special case of the Token Swapping
 *  problem which we're trying to solve!)
 *  (Of course, also maybe more than n+1 "obvious" ways because paths
 *  from v[i] to v[i+1] might not be unique).
 *
 *  It's important that the overall effect of the complete cycle
 *  doesn't move any OTHER tokens, so that we can GUARANTEE that
 *  the final TrivialTSA solution really does terminate in all cases.
 *
 *  We can "swap along" the path v(i), v(i+1), ..., v(i+n) for any 0 <= i <= n
 *  (regarding the v indices as wrapping around cyclicly,
 *  i.e. reducing (i+n) mod (n+1).)
 *
 *  This finds a choice giving the smallest number of concrete swaps,
 *  assuming no additional swap optimisation,
 *  and disregarding the tokens on the vertices.
 *
 *  It may not be the genuinely best solution
 *  because (1) swap sequences can often be optimised;
 *  (2) some of the swaps may be empty, and hence removable.
 *  But finding a truly optimal solution, taking these into account,
 *  is probably about as hard as the general token swapping problem.
 */
struct CyclicShiftCostEstimate {
  /** A simple estimate of how many swaps will be needed. */
  size_t estimated_concrete_swaps = 0;

  /** If the stored vertices are v[0], v[1], ..., v[n],
   *  this is the value of i such that swapping along the abstract path
   *  v[i], v[i+1], ..., v[i+n] gives the smallest number of swaps.
   *  (Remembering that each abstract move is v[j] -> v[j+1]).
   */
  size_t start_v_index = std::numeric_limits<size_t>::max();

  /** Calculate the data upon construction.
   *  @param vertices The list of vertices, in order, for a cyclic shift.
   *    Must have size >= 2.
   *  @param distances An object to calculate distances (we don't need to know
   *    WHICH path between vertices will be used, at this stage).
   */
  CyclicShiftCostEstimate(
      const std::vector<size_t>& vertices, DistancesInterface& distances);
};

}  // namespace tsa_internal
}  // namespace tket
