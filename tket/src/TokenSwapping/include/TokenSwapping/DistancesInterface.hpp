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

#include <cstdint>
#include <vector>

namespace tket {

/** What is the distance between any two vertices on a graph?
 *  To save time and cope with larger, sparse graphs, it may
 *  calculate distances only when required.
 */
class DistancesInterface {
 public:
  /** Not const because there might be caching, dynamic stuff going on.
   *  Find the distance between v1,v2.
   *  @param vertex1 First vertex
   *  @param vertex2 Second vertex
   *  @return distance from v1 to v2 within the graph.
   */
  virtual size_t operator()(size_t vertex1, size_t vertex2) = 0;

  /** If you KNOW a path from v1 to v2 which is shortest, then
   *  extra information about distances can be deduced from subpaths
   *  (each subpath must also be a shortest path: otherwise, the whole path
   *  would not be of minimum length).
   *  Does nothing unless overridden.
   *  @param path A sequence [v0,v1, v2, ..., vn] of vertices, KNOWN to be a
   * shortest path from v0 to vn. The caller must not call this without being
   * SURE that it really is a shortest path, or incorrect results may occur.
   */
  virtual void register_shortest_path(const std::vector<size_t>& path);

  /** If you know the neighbours of a vertex, you can tell this class
   *  and it MIGHT choose to cache the distances.
   *  Simply calls register_neighbours(v1, v2) repeatedly, unless overridden.
   *  @param vertex A vertex.
   *  @param neighbours A list of vertices adjacent to the given vertex.
   */
  virtual void register_neighbours(
      size_t vertex, const std::vector<size_t>& neighbours);

  /** Does nothing unless overridden. Stores the fact that v1,v2 are adjacent,
   *  to save later recalculation.
   *  @param vertex1 First vertex
   *  @param vertex2 Second vertex
   */
  virtual void register_edge(size_t vertex1, size_t vertex2);

  virtual ~DistancesInterface();
};

}  // namespace tket
