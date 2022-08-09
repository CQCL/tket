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
#include <utility>

#include "FilterUtils.hpp"
#include "GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

class NearNeighboursData {
 public:
  enum class Type { PATTERN_GRAPH, TARGET_GRAPH };

  NearNeighboursData(const NeighboursData& ndata, Type type);

  /** Calculated lazily, on demand; returns a sorted list of vertices
   * exactly at distance d from v. We must have d >= 2
   * (neighbours with d=1 come from the NeighboursData object instead).
   * @param v A vertex in the graph.
   * @param distance A distance d >= 2.
   * @return A reference to a sorted vector (stored within this class) of all
   * vertices at distance exactly d from v. The reference may be invalidated, of
   * course, if further non-const calls are made to this class.
   */
  const std::vector<VertexWSM>& get_vertices_at_distance(
      VertexWSM v, unsigned distance);

  /** Returns information about the degrees of vertices up to
   * the given distance from the root vertex.
   * If this is a pattern graph, returns data for vertices EXACTLY at
   * the specified distance; if a target graph, for vertices at
   * or closer than the given distance.
   * Requires the distance to be >= 2.
   */
  const FilterUtils::DegreeCounts& get_degree_counts(
      VertexWSM v, unsigned distance);

  /** Cached. Does NOT include the vertex v itself. */
  std::size_t get_n_vertices_at_max_distance(VertexWSM v, unsigned distance);

 private:
  const NeighboursData& m_ndata;
  const Type m_type;

  struct VertexData {
    /** element[i] is all the vertices at distance i+2, sorted by vertex number.
     * We don't list immediate neighbours, to save space,
     * since they're already stored in the NeighboursData object.
     */
    std::vector<std::vector<VertexWSM>> vertices_at_distance;
    std::vector<std::size_t> n_vertices_at_max_distance;

    /** Element[i] is for distance i+2.
     * If a pattern graph, includes vertices at distance EXACTLY i+2;
     * but if a target graph, includes vertices at distance <= i+2.
     */
    std::vector<FilterUtils::DegreeCounts> degree_counts_for_distance;
  };

  // Element[i] gives data for vertex i.
  // Lazy initialisation, but will be resized correctly upon construction.
  std::vector<VertexData> m_data;

  std::set<VertexWSM> m_vertices_workset;

  // Will be used to fill "degree_counts_at_distance".
  std::map<std::size_t, std::size_t> m_work_map;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
