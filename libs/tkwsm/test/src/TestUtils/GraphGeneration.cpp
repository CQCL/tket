// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "GraphGeneration.hpp"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

GraphGeneration::LimitedSizeGraphGeneral::LimitedSizeGraphGeneral(
    LimitedSizeGraphSeed seed)
    // There are 2 bits per edge (0 to denote no edge, 01,10,11 denote weights
    // w1,w2,w3). 9*8/2 = 36 > 32, so 4 edges are missing. No point in shoving
    // in loads more vertices; there'd just be lots of permanently missing
    // edges, and so lots of isolated vertices.
    : max_number_of_vertices(9) {
  const std::array<WeightWSM, 3> weights{1, 4, 9};
  for (unsigned ii = 0; ii < max_number_of_vertices; ++ii) {
    for (unsigned jj = ii + 1; jj < max_number_of_vertices; ++jj) {
      const auto code = seed & 3;
      if (code > 0) {
        data[get_edge(ii, jj)] = weights.at(code - 1);
      }
      seed >>= 2;
    }
  }
}

GraphEdgeWeights GraphGeneration::get_cycle(
    unsigned vertices, bool mix_weights) {
  GraphEdgeWeights data;
  for (VertexWSM v = 0; v < vertices; ++v) {
    data[get_edge(v, (v + 1) % vertices)] = mix_weights ? v + 1 : 1;
  }
  return data;
}

GraphEdgeWeights GraphGeneration::get_line(
    unsigned vertices, bool mix_weights) {
  GraphEdgeWeights data;
  for (VertexWSM v = 0; v + 1 < vertices; ++v) {
    data[get_edge(v, v + 1)] = mix_weights ? v + 1 : 1;
  }
  return data;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
