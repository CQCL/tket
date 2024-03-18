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

#include "PlacementCostModelInterface.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

/** For testing only; for simplicity, same width and height;
 * will be auto-deduced from the number of weights.
 */
class WeightedSquareGrid : public PlacementCostModelInterface {
 public:
  /** The number of weights will be used to deduce the size of the grid
   * (we reuqire it to have same width and height).
   * @param weights weights[i] is the edge weight from vertex i to its parent.
   * (Thus weights[0], weights[1] are "dummy" weights, as 0,1 have no parent).
   * @param number_of_primitive_gates_in_swap How many primitive 2-qubit gates
   * (with cost equal to the edge weight) does it take to make a single SWAP
   * gate, in our model?
   */
  explicit WeightedSquareGrid(
      std::vector<WeightWSM> weights,
      unsigned number_of_primitive_gates_in_swap = 3);

  virtual GraphEdgeWeights get_graph_data() const override;

 private:
  const std::vector<WeightWSM> m_weights;
  const unsigned m_width;
  const unsigned m_number_of_vertices;

  // Extremely crude: to get from (x1,y1) to (x2,y2), we'll just try various
  // "random" V,H strings like "VHHVHVHVVVHVHHHVHVVHVHHHV..." etc.
  // encoded as these bits, 0 for H.
  // If x1<x2 then H will mean (dx,dy) = (1,0), whereas if x1>x2
  // then H will mean (dx,dy) = (-1,0), etc.
  // We then try a fixed number of these paths and choose the path
  // with the lowest total weight.
  const std::vector<std::uint64_t> m_vertical_horizontal_patterns;

  VertexWSM get_vertex(unsigned x, unsigned y) const;

  // Return the index i such that
  //      m_weights[i] = Weight( (x y), (x+1 y) ).
  unsigned get_horizontal_weight_index(unsigned x, unsigned y) const;

  // Return the index i such that
  //      m_weights[i] = Weight( (x y), (x y+1) ).
  unsigned get_vertical_weight_index(unsigned x, unsigned y) const;

  std::pair<unsigned, unsigned> get_xy(VertexWSM v) const;

  // KEY: is the pair (v1,v2), NOT required to have v1<v2.
  //    But we shall require Path(v1, v2) to be the reverse of Path(v2, v1),
  //    and always calculate Path(u, v) with u<v first,
  //    and take Path(v,u) to be the reverse of that.
  //    Thus, we can calculate paths lazily, in any order,
  //    and the result is the same.
  //    There are no relations between subpaths, i.e.
  //      Path(v1, v2),  Path(v2, v3) are not necessarily subpaths of Path(v1,
  //      v3).
  //    This is actually realistic: when ACTUALLY doing token swapping,
  //    there's no "relation" between "similar" vertex mappings and returned
  //    token swapping solutions - although these are not precise notions.
  //    We think of the function
  //          {desired vertex mapping} -> {token swapping solution}
  //    as "discontinuous", although of course that's meaningless
  //    as everything is discrete.
  //
  // VALUE: is the path from v1 to v2.
  mutable std::map<std::pair<VertexWSM, VertexWSM>, Path> m_paths;

  // Fills the given path with the best path found from v1 to v2,
  // using path_work_vector as temporary storage; so it will be messed up.
  void fill_path(
      Path& path, Path& path_work_vector, VertexWSM vertex1,
      VertexWSM vertex2) const;

  virtual const Path& get_path_to_use(
      VertexWSM vertex1, VertexWSM vertex2) const override;
};

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
