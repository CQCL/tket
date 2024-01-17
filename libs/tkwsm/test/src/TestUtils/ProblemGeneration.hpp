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
#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct ProblemGeneration {
  /* Imagine a square grid of 25 vertices.
  There are 2*(5*4) = 40 edges in total. If each edge could have
  one of 4 weights, we could represent the data in 80 bits.
  For a 64-bit int, let's instead allow 32 edges (with the other edges
  just being set to weight 1), so it fits inside a 64-bit uint.

  We'll also allow the 4 weights to be specified. E.g.,
  weights 1,2,3,4 are more flexible than 1,10,100,1000,
  meaning that if one graph has weights 1,2,3,4 and we try to embed
  another graph into it, it is probably harder than if we take the same
  graph but replace the weights 2->10, 3->100, 4->1000.
  The reason being, that we are far less likely to use weight 1000
  in an optimal solution, so it will be pruned more quickly.
  */
  // Element[0]: the 64-bit uint representing the weights
  // Elements 1,2,3: the edge weights w1,w2,w3
  // (it's assumed that the smallest weight w0 is 1).
  // So, bits 00,01,10,11 will represent weights w0,w1,w2,w3.
  typedef std::vector<std::uint_fast64_t> EncodedSquareGrid;
  static GraphEdgeWeights get_target_graph_for_encoded_square_grid(
      const EncodedSquareGrid& encoding);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
