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

#include "GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

class NearNeighboursData {
 public:
  explicit NearNeighboursData(const NeighboursData& ndata);

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

  /** Cached. Does NOT include the vertex v itself. */
  std::size_t get_n_vertices_at_max_distance(VertexWSM v, unsigned distance);

 private:
  const NeighboursData& m_ndata;

  struct VertexData {
    /** element[i] is all the vertices at distance i+2, sorted by vertex number.
     * We don't list immediate neighbours, to save space,
     * since they're already stored in the NeighboursData object.
     */
    std::vector<std::vector<VertexWSM>> vertices_at_distance;
    std::vector<std::size_t> n_vertices_at_max_distance;
  };

  std::set<VertexWSM> m_vertices_workset;

  // KEY: the vertex  VALUE: data about that vertex (including lazy
  // initialisation).
  std::map<VertexWSM, VertexData> m_data;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
