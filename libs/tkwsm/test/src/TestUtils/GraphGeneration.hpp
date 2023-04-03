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

#pragma once

#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct GraphGeneration {
  // We will generate a graph on vertices {0,1,...,N-1}
  // from a single large int of this type,
  // in a fully portable way, where N is fixed.
  typedef std::uint64_t LimitedSizeGraphSeed;

  // Edge weights can vary. Obviously the number of vertices
  // is not so large, since we've only got 64 bits of information.
  struct LimitedSizeGraphGeneral {
    unsigned max_number_of_vertices;
    GraphEdgeWeights data;

    // Uses the bits of the seed to decide the subset of N(N-1)/2 edges
    // which are to exist; N is fixed, each edge weight is 1.
    explicit LimitedSizeGraphGeneral(LimitedSizeGraphSeed seed);
  };

  static GraphEdgeWeights get_limited_size_graph(LimitedSizeGraphSeed seed);

  static GraphEdgeWeights get_cycle(unsigned vertices, bool mix_weights = true);
  static GraphEdgeWeights get_line(unsigned vertices, bool mix_weights = true);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
