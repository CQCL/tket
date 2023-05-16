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

#include "ProblemGeneration.hpp"

#include <catch2/catch_test_macros.hpp>

#include "SquareGridGeneration.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

GraphEdgeWeights ProblemGeneration::get_target_graph_for_encoded_square_grid(
    const EncodedSquareGrid& encoding) {
  REQUIRE(encoding.size() >= 4);
  for (unsigned ii = 1; ii <= 3; ++ii) {
    // force the edge weights to be strictly increasing.
    REQUIRE(encoding[ii] > 1);
    if (ii > 1) {
      REQUIRE(encoding[ii] > encoding[ii - 1]);
    }
  }
  // Make everything comfortably small for 32 bits, etc.
  REQUIRE(encoding[3] <= 100000);
  SquareGrid grid;
  grid.width = 4;
  grid.height = 4;
  std::uint_fast64_t edges_encoding = encoding[0];

  const auto set_weights = [&edges_encoding,
                            &encoding](std::vector<WeightWSM>& weights) {
    weights.resize(20);
    for (unsigned ii = 0; ii < weights.size(); ++ii) {
      const auto code = edges_encoding & 3;
      const WeightWSM weight = (code == 0) ? 1 : encoding[code];
      weights[ii] = weight;
      edges_encoding >>= 2;
    }
  };

  set_weights(grid.horiz_weights);
  set_weights(grid.vert_weights);
  return grid.get_graph_edge_weights();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
