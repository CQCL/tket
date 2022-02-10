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
#include <set>

#include "GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

/** For finding all vertices within a specific distance of a vertex. */
class CloseNeighboursCalculator {
 public:
  /** Element[i] lists all vertices at distance exactly i+1
   * from the given vertex.
   */
  typedef std::vector<std::vector<VertexWSM>> CloseVertices;

  /** Fills close_vertices. Sorts each vector (necessary for some algorithms).
   * Leaves empty vectors in the higher levels if it runs out of vertices.
   * @param num_levels The final desired size of CloseVertices
   * @param close_vertices The data to overwrite.
   * @param ndata An object for computing neighbours in the graph.
   * @param v The seed vertex in the graph, from which distances are calculated.
   */
  void operator()(
      unsigned num_levels, CloseVertices& close_vertices,
      const NeighboursData& ndata, VertexWSM v);

 private:
  std::set<VertexWSM> m_vertices_seen;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
