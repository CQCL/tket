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

#include "TestWeightedGraphData.hpp"

#include <catch2/catch_test_macros.hpp>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

GraphEdgeWeights get_graph_data(
    RNG& rng, unsigned num_vertices, unsigned approx_edges,
    WeightWSM approx_min_weight, WeightWSM approx_max_weight) {
  REQUIRE(approx_min_weight > 0);
  REQUIRE(approx_max_weight >= approx_min_weight);
  GraphEdgeWeights graph_data;
  const std::size_t max_vertex = num_vertices - 1;
  for (unsigned ii = 0; ii < num_vertices; ++ii) {
    for (;;) {
      const unsigned jj = rng.get_size_t(max_vertex);
      if (ii != jj) {
        graph_data[get_edge(ii, jj)] = approx_min_weight;
        break;
      }
    }
  }
  approx_edges -= graph_data.size();
  const WeightWSM weight_diff = approx_max_weight - approx_min_weight;
  for (unsigned nn = 1; nn <= approx_edges; ++nn) {
    const unsigned ii = rng.get_size_t(max_vertex);
    const unsigned jj = rng.get_size_t(max_vertex);
    if (ii == jj) {
      continue;
    }
    // Steadily increase the weight from the min to the max
    graph_data[get_edge(ii, jj)] =
        approx_min_weight + (nn * weight_diff) / approx_edges;
  }
  return graph_data;
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
