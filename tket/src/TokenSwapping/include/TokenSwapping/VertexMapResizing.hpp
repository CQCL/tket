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
#include <map>
#include <set>
#include <vector>

#include "TokenSwapping/NeighboursInterface.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** If a vertex mapping { u -> v } has too few vertices, try to add extra
 * vertices, fixed by the new mapping, to get to the desired size. This may
 * allow extra optimisations to be found in the table. E.g., imagine a vertex in
 * a graph which is not moved by the mapping. Imagine that removing it makes the
 * graph disconnected. If the desired mapping moves a token
 * between different components, it is then impossible for any swap
 * sequence within the subgraph to perform that mapping.
 * However, adding the vertex back makes it possible.
 *
 * If instead there are too many vertices to look up in the table, it tries
 * to remove vertices which are fixed by the mapping to get it down to size.
 */
class VertexMapResizing : public NeighboursInterface {
 public:
  /** Store a Neighbours object, to be used throughout when required to find
   * all neighbours of a given vertex. The caller must ensure that the
   * object remains valid.
   * @param neighbours The object to calculate neighbours of a vertex.
   */
  explicit VertexMapResizing(NeighboursInterface& neighbours);

  /** Gets the data by calling the NeighboursInterface object which was passed
   * into the constructor. HOWEVER, it does internal caching, so doesn't call it
   * multiple times.
   * @param vertex A vertex in the graph.
   * @return A cached list of neighbours of that vertex, stored internally.
   */
  virtual const std::vector<size_t>& operator()(size_t vertex) override;

  /** The result of resizing a mapping by deleting fixed vertices if too big,
   * or adding new vertices if too small.
   */
  struct Result {
    /** It is still a success if we have fewer vertices than the desired number
     * (as this can still be looked up in the table). However, it's a failure if
     * there are too many vertices (which than cannot be looked up).
     */
    bool success;

    /** If successful, the edges of the subgraph containing only the vertices in
     * the new mapping. */
    std::vector<Swap> edges;
  };

  /** The mapping may be altered, even upon failure, so obviously the caller
   * should make a copy if it needs to be preserved. Increase the map size as
   * much as possible if too small (still a success even if it cannot reach the
   * size). Decrease the size if too large (and not reaching the szie is then a
   * failure). Newly added or removed vertices are all fixed, i.e. map[v]=v.
   * @param mapping The mapping which will be altered and returned by reference.
   * @param desired_size The size we wish to reach, or as close as possible if
   * the mapping is currently too small.
   */
  const Result& resize_mapping(
      VertexMapping& mapping, unsigned desired_size = 6);

 private:
  NeighboursInterface& m_neighbours;
  Result m_result;

  // KEY: a vertex. VALUE: all its neighbours.
  std::map<size_t, std::vector<size_t>> m_cached_neighbours;
  std::set<Swap> m_cached_full_edges;

  /** How many edges join the given vertex to other existing vertices?
   * @param mapping The current vertex permutation which we may expand or
   * contract.
   * @param vertex A vertex which may or may not be already within the mapping.
   * @return The total number of edges within the LARGER graph joining the
   * vertex to other vertices within the mapping.
   */
  size_t get_edge_count(const VertexMapping& mapping, size_t vertex);

  /** Try to add a single new fixed vertex to the mapping, i.e. a new v with
   * map[v]=v.
   * @param mapping The current vertex permutation which we wish to expand by
   * one vertex.
   */
  void add_vertex(VertexMapping& mapping);

  /** Try to remove a single vertex within the mapping, but only if it is fixed,
   * i.e. map[v]==v.
   * @param mapping The current vertex permutation which we wish to shrink by
   * one vertex.
   */
  void remove_vertex(VertexMapping& mapping);

  /** Within the m_result object, fill "edges" for the new mapping. */
  void fill_result_edges(const VertexMapping& mapping);
};

}  // namespace tsa_internal
}  // namespace tket
