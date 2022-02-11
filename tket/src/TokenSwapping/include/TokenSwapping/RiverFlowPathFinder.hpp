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
#include "PathFinderInterface.hpp"
#include "RNG.hpp"

namespace tket {
namespace tsa_internal {

/** Think of flowing water: if it has already flowed through, it creates
 *  channels along which it is more likely to flow next time.
 *  We do a similar idea: the PURPOSE is to try to make paths overlap;
 *  if we move tokens along paths with many edges in common, it is more likely
 *  that some basic swap optimisation will reduce the number of swaps.
 *  (Disjoint swaps are the worst kind to optimise, of course;
 *  no reduction is possible).
 *
 *  This is supposed to be reasonably fast. Repeated calls to operator()(v1,v2)
 *  are likely to return the same path, but may change slightly over time.
 */
class RiverFlowPathFinder : public PathFinderInterface {
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
  virtual void reset() override;

  /** Get the path from v1 to v2. As always, may change over time;
   *  path(v1, v2) is NOT necessarily the reverse of path(v2, v1).
   *  @param vertex1 First vertex v1.
   *  @param vertex2 Second vertex v2.
   *  @return A list of vertices, starting with v1 and ending with v2,
   *    giving a shortest path from v1 to v2.
   */
  virtual const std::vector<size_t>& operator()(
      size_t vertex1, size_t vertex2) override;

  virtual ~RiverFlowPathFinder();

  /** We really do want to know which edges have been used in the solution so
   * far, that's the whole point of this class.
   * @param vertex1 First vertex v1 of an edge v1-v2 that was used in the
   * solution.
   * @param vertex2 Second vertex v2 of the edge.
   */
  virtual void register_edge(size_t vertex1, size_t vertex2) override;

  /** Returns true for this object, since we definitely do want to remember
   * previous edges.
   * @return True, always, for this class.
   */
  virtual bool edge_registration_has_effect() const override;

 private:
  struct Impl;
  std::unique_ptr<Impl> m_pimpl;
};

}  // namespace tsa_internal
}  // namespace tket
