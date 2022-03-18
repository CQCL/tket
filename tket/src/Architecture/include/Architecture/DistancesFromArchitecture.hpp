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

#include "ArchitectureMapping.hpp"
#include "TokenSwapping/DistancesInterface.hpp"
#include "TokenSwapping/SwapFunctions.hpp"

namespace tket {

/** Directly get distances from an architecture object,
 *  but evaluated lazily.
 */
class DistancesFromArchitecture : public DistancesInterface {
 public:
  /** The ArchitectureMapping object already handles the Node <-> vertex size_t
   * conversion.
   *  @param arch_mapping Object containing a reference to an Architecture,
   *    which has decided upon Node <-> vertex size_t conversions.
   */
  explicit DistancesFromArchitecture(const ArchitectureMapping& arch_mapping);

  /** Get the distance from v1 to v2. Throws if distinct vertices return
   *  distance 0, which probably means a disconnected graph.
   *  @param vertex1 First vertex
   *  @param vertex2 Second vertex
   *  @return distance from v1 to v2 within the Architecture graph, throwing if
   * they are disconnected (so the distance is +infinity).
   */
  virtual size_t operator()(size_t vertex1, size_t vertex2) override;

  /** May save computation time later; by some method, the caller
   *  has determined a path from v1 to v2, and hence all along the path
   *  we know the distance between any two points.
   *  However, avoids quadratic time blowup by discarding some information
   *  for long paths.
   *  @param path A sequence [v0,v1, v2, ..., vn] of vertices, KNOWN to be a
   * shortest path from v0 to vn. The caller must not call this without being
   * SURE that it really is a shortest path, or incorrect results may occur.
   */
  virtual void register_shortest_path(const std::vector<size_t>& path) override;

  /** The caller has determined that v1, v2 are adjacent, and therefore
   *  the distance from v1 to v2 equals one. Store this.
   *  @param vertex1 First vertex
   *  @param vertex2 Second vertex
   */
  virtual void register_edge(size_t vertex1, size_t vertex2) override;

 private:
  /** Reference to the original object passed into the constructor;
   *  the caller must ensure that it remains valid and unchanged.
   */
  const ArchitectureMapping& m_arch_mapping;

  /** The key is the vertex pair (v1, v2), but always sorted with v1<v2
   *  to use half the space.
   */
  std::map<Swap, size_t> m_cached_distances;

  /** The main register_shortest_path wraps around this; we want to avoid
   *  quadratic timings growth by cutting off long paths.
   *  This stores the quadratic number of distances between all vertex pairs
   *  within the given subpath.
   *  @param path A sequence [v0,v1, v2, ..., vn] of vertices,
   *    KNOWN to be a shortest path from v0 to vn.
   *  @param begin The first index in path to use.
   *  @param end Like end(), an index one past the last index in path to use.
   */
  void register_shortest_path_with_limits(
      const std::vector<size_t>& path, size_t begin, size_t end);
};

}  // namespace tket
