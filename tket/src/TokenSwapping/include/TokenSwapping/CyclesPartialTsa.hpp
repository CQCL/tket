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

#include "CyclesCandidateManager.hpp"
#include "PartialTsaInterface.hpp"

namespace tket {
namespace tsa_internal {

/** A partial TSA (token swapping algorithm), similar to the cycle-finding
 *  algorithm as described in the 2016 paper "Approximation and Hardness of
 * Token Swapping" by T.Miltzow and others:
 *
 *  https://arxiv.org/abs/1602.05150
 *
 *  However, our algorithm differs from the paper in several important ways:
 *   (1) We also accept partial mappings, i.e. there might not be a token on
 * every vertex. (2) It is only a partial TSA, not a full TSA (it may give up
 * early). Thus, a full end-to-end solution must combine this with another TSA).
 *  (3) It does not detect long cycles. (4) It never returns "unhappy
 * swaps", and is strictly monotonic: L, the sum of distances of a vertex to its
 * target, either strictly decreases, or stays the same and no swaps are
 * performed. (However, within each cycle, it is possible to have bad swaps
 *          which don't decrease L much, or even increase L, as long as the
 * overall result is a decrease in L). (5) The closing edge of a cycle is not
 * required to exist in the graph.
 *
 * Thus, neither this nor the algorithm in the paper is a generalisation of or
 * necessarily better/worse than the other.
 *
 * One of the ideas in the 2016 paper is to detect good cycles (cyclic shifts)
 * v0->v1-> ... ->vn->v0, by searching for cycles in a directed graph.
 * It is guaranteed to find cycles if they exist, no matter the length. So, by
 * (3), it is better than ours in this sense. However, we don't need the full
 * cycle to exist, by (5), since we swap along the path [v0,v1,v2,...,vn].
 * Hence, the paper algorithm is worse than ours in this respect. Regarding (2)
 * and (4), the paper is better than ours because it always completes.
 */
class CyclesPartialTsa : public PartialTsaInterface {
 public:
  CyclesPartialTsa();

  /** Calculate a solution to improve the current token configuarion,
   *  add the swaps to the list, and carry out the swaps on "vertex_mapping".
   *  We don't need a path finder because the cycles are built up one vertex
   *  at a time, so we only need distances and neighbours.
   *  There is no point in calling this multiple times;
   *  it will continue until EITHER all tokens are home, OR it gives up.
   *  @param swaps The list of swaps to add to.
   *  @param vertex_mapping The current state, giving vertex->target mappings.
   *    Will be updated if any new swaps are performed.
   *  @param distances An object to calculate distances between vertices.
   *  @param neighbours An object to calculate the neighbours of a vertex.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices. (Of course, paths might not be unique if the graph
   *    is not a tree, so it is an important part of the heuristics that
   *    the returned paths are fairly "consistent", i.e. "nearby" vertex pairs
   *    should return "nearby" paths).
   */
  virtual void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& neighbours,
      RiverFlowPathFinder& path_finder) override;

 private:
  /** Stores cycles, and controls the growth and discarding of cycles.
   *  We grow the cycles one vertex at a time until we reach a good cycle
   *  which is worth turning into swaps.
   *  If we never find a good cycle then we give up without returning a
   * solution.
   */
  CyclesGrowthManager m_growth_manager;

  /** Controls the final selection of cycles to perform. Once we've found
   *  some good cycles, we may not be able to perform all of them
   *  (because they might not be disjoint, so interfere with each other).
   *  We may not even want to perform them all, depending upon the options.
   */
  CyclesCandidateManager m_candidate_manager;

  /** "append_partial_solution" simply loops, calling this repeatedly until
   *  it gives up, or all tokens are home.
   */
  void single_iteration_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& neighbours);
};

}  // namespace tsa_internal
}  // namespace tket
