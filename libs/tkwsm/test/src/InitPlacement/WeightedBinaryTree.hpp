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

#include <string>

#include "PlacementCostModelInterface.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

/** Just for testing; we want a reaonsably large graph where
 * sensible token swapping is easy to work out.
 * (And because it's a TREE, shortest paths are unique).
 *
 * Thus this gives a reasonable simplified model for the cost of applying
 * a sequence of gates, once we've got an initial qubit placement.
 *
 * However, this is only for very quick sanity tests; a full benchmark
 * should be carried out with many different graphs,
 * token swapping strategies, etc. (Even GIVEN an initial placement,
 * and GIVEN a sequence of gates to apply, there may be many
 * possible reorderings from parallel gates;
 * so there are probably many possible solutions, even for trees.
 * It is presumably impossible to find an OPTIMUM solution
 * in less than exponential time in general. Need to research this!)
 *
 * Represents, implicitly, a binary tree. The vertices are labelled
 * {1,2,3,...}. Let 1 be the root vertex, at level 0.
 * As we descend down, vertex v has left child 2v, right child 2v+1.
 * Thus, the trees look like:
 *
 *            1
 *          /   \
 *         /     \
 *        /       \
 *       /         \
 *      2           3
 *    /   \       /   \
 *   4     5     6     7
 *  / \   / \   / \   / \
 * 8   9 10 11 12 13 14  15
 *
 * ...etc. It's very easy to find the level of a vertex
 * from its binary representation.
 */
class WeightedBinaryTree : public PlacementCostModelInterface {
 public:
  /**
   * @param weights weights[i] is the edge weight from vertex i to its parent.
   * (Thus weights[0], weights[1] are "dummy" weights, as 0,1 have no parent).
   * @param number_of_primitive_gates_in_swap How many primitive 2-qubit gates
   * (with cost equal to the edge weight) does it take to make a single SWAP
   * gate, in our model?
   */
  explicit WeightedBinaryTree(
      std::vector<WeightWSM> weights,
      unsigned number_of_primitive_gates_in_swap = 3);

  // Calculate the edges and weights in standard format.
  virtual GraphEdgeWeights get_graph_data() const override;

  // The vertices are 1,2,...,N. Returns N.
  // This is, of course, equal to the number of vertices.
  unsigned get_max_vertex_number() const;

 private:
  // Element[i] is the edge weight from vertex i to its parent
  // (with weight 0 for i=0,1, which are "dummy" values; vertices 0,1 have no
  // parent.)
  const std::vector<WeightWSM> m_weights;

  // For the path u -> v, e.g. [u, x, y, z, a, b, c, v],
  // it's convenient to build it up from both ends:
  // [u, x, y, ...]  and [v, c, b, a, ...] until they meet in the middle.
  // This will be the tail, i.e. starting with the end vertex;
  // thus the order will later need to be reversed.
  // The weights will also be shifted; thus element[0], corresponding to v,
  // will be (v, weight(v--c)), rather than (v,0). This is more convenient
  // when reversing and concatenating to make the complete path.
  mutable Path m_path_work_vector;
  mutable Path m_tail_path_work_vector;

  // Choose a path between the VERTICES v1, v2, and fill m_path_work_vector.
  void fill_path(VertexWSM vertex1, VertexWSM vertex2) const;

  virtual const Path& get_path_to_use(
      VertexWSM vertex1, VertexWSM vertex2) const override;
};

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
