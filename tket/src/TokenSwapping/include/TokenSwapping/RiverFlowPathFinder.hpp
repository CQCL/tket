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
#include <memory>
#include <optional>

#include "DistancesInterface.hpp"
#include "NeighboursInterface.hpp"
#include "Utils/RNG.hpp"

namespace tket {
namespace tsa_internal {

/** Given two vertices in a graph, find a shortest path between them;
 * of course paths might not be unique.
 * The aim is to make paths overlap;
 * if we move tokens along paths with many edges in common, it is more likely
 * that some basic swap optimisation will reduce the number of swaps.
 * (Disjoint swaps are the worst kind to optimise, of course;
 * no reduction is possible).
 *
 * We think of flowing water: if water has already flowed through,
 * it creates channels along which it is more likely to flow next time.
 * We do a similar thing: by remembering which edges have already been used,
 * whenever we have a choice of edge to continue a path, choose one which
 * has already been used frequently.
 *
 * Repeated calls to operator()(v1,v2)
 * are likely to return the same path, but may change slightly over time.
 */
class RiverFlowPathFinder {
 public:
  /** All the objects should remain valid throughout
   *  the lifetime of this object.
   *  @param distances An object to calculate distances between vertices.
   *  @param neighbours An object to calculate adjacent vertices to any given
   * vertex.
   *  @param rng A source of (pseudo) randomness.
   */
  RiverFlowPathFinder(
      DistancesInterface& distances, NeighboursInterface& neighbours, RNG& rng);

  /** For reuse in different problems (but still the same architecture;
   *  the same "distances" and "neighbours" objects are used),
   *  but constructing paths anew
   *  (which is appropriate because completely different problems will
   *  probably need different paths). This also resets the RNG with its
   *  default seed, for better reproducibility.
   *
   *  (This may be suitable for simulated annealing-type algorithms
   *  which involve solving with many different token positions, i.e.
   *  partially finished problems, even though the end-to-end problem
   *  is the same).
   */
  void reset();

  /** Get the path from v1 to v2. May change over time, and
   *  path(v1, v2) is NOT necessarily the reverse of path(v2, v1).
   *  @param vertex1 First vertex v1.
   *  @param vertex2 Second vertex v2.
   *  @return A list of vertices, starting with v1 and ending with v2,
   *    giving a shortest path from v1 to v2.
   */
  const std::vector<size_t>& operator()(size_t vertex1, size_t vertex2);

  ~RiverFlowPathFinder();

  /** Whenever an edge is used, i.e. we swap tokens along it, tell this
   * object; the proper functioning of this class depends on
   * knowing which edges have been used in the solution so far.
   * @param vertex1 First vertex v1 of an edge v1-v2 that was used in the
   * solution.
   * @param vertex2 Second vertex v2 of the edge.
   */
  void register_edge(size_t vertex1, size_t vertex2);

 private:
  struct Impl;
  /** Pimpl idiom. */
  std::unique_ptr<Impl> m_pimpl;
};

}  // namespace tsa_internal
}  // namespace tket
