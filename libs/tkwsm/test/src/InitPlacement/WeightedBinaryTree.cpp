// Copyright Quantinuum
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

#include "WeightedBinaryTree.hpp"

#include <catch2/catch_test_macros.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

WeightedBinaryTree::WeightedBinaryTree(
    std::vector<WeightWSM> weights, unsigned number_of_primitive_gates_in_swap)
    : m_weights(std::move(weights)) {
  m_number_of_primitive_gates_in_swap = number_of_primitive_gates_in_swap;
  // We want a nontrivial size...
  REQUIRE(m_weights.size() >= 6);

  // Start with dummy weights!
  // Vertex 0 doesn't exist, and vertex 1 has no parent.
  REQUIRE(m_weights[0] == 0);
  REQUIRE(m_weights[1] == 0);
  REQUIRE(m_number_of_primitive_gates_in_swap >= 1);
  REQUIRE(m_number_of_primitive_gates_in_swap <= 100);
}

unsigned WeightedBinaryTree::get_max_vertex_number() const {
  return m_weights.size() - 1;
}

GraphEdgeWeights WeightedBinaryTree::get_graph_data() const {
  GraphEdgeWeights result;
  for (unsigned ii = 2; ii < m_weights.size(); ++ii) {
    const VertexWSM child_vertex = ii;
    const VertexWSM parent_vertex = child_vertex / 2;
    result[std::make_pair(parent_vertex, child_vertex)] = m_weights[ii];
  }
  return result;
}

void WeightedBinaryTree::fill_path(VertexWSM vertex1, VertexWSM vertex2) const {
  REQUIRE(vertex1 != vertex2);
  m_path_work_vector.clear();
  m_tail_path_work_vector.clear();
  m_path_work_vector.emplace_back(vertex1, 0);
  m_tail_path_work_vector.emplace_back(vertex2, 0);

  // The path goes up, then goes down; thus the vertex numbers decrease,
  // then increase. Both vectors have decreasing vertex number.
  // Thus when they meet, they've collided,
  // and we then reverse the tail vertices list.
  for (;;) {
    const VertexWSM head_v = m_path_work_vector.back().first;
    const VertexWSM tail_v = m_tail_path_work_vector.back().first;
    if (head_v == tail_v) {
      break;
    }
    if (head_v < tail_v) {
      m_tail_path_work_vector.back().second = m_weights.at(tail_v);
      m_tail_path_work_vector.emplace_back(tail_v / 2, 0);
    } else {
      m_path_work_vector.emplace_back(head_v / 2, m_weights.at(head_v));
    }
  }
  REQUIRE(
      m_path_work_vector.back().first == m_tail_path_work_vector.back().first);

  // Erase the common vertex.
  m_tail_path_work_vector.pop_back();

  // Now treat the tail as a stack.
  while (!m_tail_path_work_vector.empty()) {
    m_path_work_vector.push_back(m_tail_path_work_vector.back());
    m_tail_path_work_vector.pop_back();
  }
}

const WeightedBinaryTree::Path& WeightedBinaryTree::get_path_to_use(
    VertexWSM vertex1, VertexWSM vertex2) const {
  fill_path(vertex1, vertex2);
  return m_path_work_vector;
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
