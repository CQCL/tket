// Copyright 2019-2024 Cambridge Quantum Computing
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
  explicit NearNeighboursData(const NeighboursData& ndata);

  std::size_t get_number_of_vertices() const;

  /** Calculated lazily, on demand; returns a sorted list of vertices
   * exactly at distance d from v. We must have d >= 2
   * (neighbours with d=1 come from the NeighboursData object instead).
   * @param v A vertex in the graph.
   * @param distance A distance d >= 2.
   * @return A reference to a sorted vector (stored within this class) of all
   * vertices at distance exactly d from v. The reference may be invalidated, of
   * course, if further non-const calls are made to this class.
   */
  const boost::dynamic_bitset<>& get_vertices_at_exact_distance(
      VertexWSM v, unsigned distance);

  /** Not including the input vertex itself, i.e. 1 <= dist(u,v) <= d. */
  const boost::dynamic_bitset<>& get_vertices_up_to_distance(
      VertexWSM v, unsigned distance);

  /** Returns information about the degrees of vertices at
   * the given distance from the root vertex.
   * Requires the distance to be >= 1.
   */
  const FilterUtils::DegreeCounts& get_degree_counts_at_exact_distance(
      VertexWSM v, unsigned distance);

  /** If combining all the vertices up to the given distance
   * (not including v itself) into a single set of vertices,
   * return a complete count of all the vertex degrees.
   */
  const FilterUtils::DegreeCounts& get_degree_counts_up_to_distance(
      VertexWSM v, unsigned distance);

  /** The number of vertices V' with 1 <= Dist(V,V') <= distance.
   * Thus, does NOT include the vertex v itself
   * (more convenient for applications).
   */
  std::size_t get_n_vertices_up_to_distance(VertexWSM v, unsigned distance);

  std::size_t get_n_vertices_at_exact_distance(VertexWSM v, unsigned distance);

 private:
  const NeighboursData& m_ndata;

  /** In each vector, element[i] is the data for distance i+1.
   * Some of this data is needed for pattern graphs but not
   * target graphs, or vice versa; but for simplicity don't bother
   * separating them out carefully by intended use.
   */
  struct VertexData {
    // This will be the "primary" data, i.e. everything else is
    // calculated from it.
    std::vector<boost::dynamic_bitset<>> vertices_at_exact_distance;
    std::vector<boost::dynamic_bitset<>> vertices_up_to_distance;
    std::vector<FilterUtils::DegreeCounts> degree_counts_for_exact_distance;
    std::vector<FilterUtils::DegreeCounts> degree_counts_up_to_max_distance;
  };

  // Element[v] gives data for vertex v. Lazy initialisation.
  std::vector<VertexData> m_data;

  // Will be used to fill the degree_counts.
  // The maximum size this can reach is the number of vertices.
  // Seems crude, but actually clearing and refilling this is
  // probably faster than using a std::map, over time.
  // (std::set and std::map are a LOT slower than std::vector!)
  std::vector<std::size_t> m_degree_counts_work_vector;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
