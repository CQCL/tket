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

#include <catch2/catch.hpp>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("test neighbours_data on cycles") {
  const unsigned cycle_length = 5;
  const std::vector<std::size_t> deg_seqs{2, 2};
  GraphEdgeWeights edge_weights;
  std::vector<VertexWSM> neighbours_recalc(2);

  for (VertexWSM ii = 0; ii < cycle_length; ++ii) {
    const VertexWSM jj = (ii + 1) % cycle_length;
    edge_weights[get_edge(ii, jj)] = ii + jj;
  }
  const NeighboursData ndata(edge_weights);

  // VertexWSM funcs
  for (VertexWSM ii = 0; ii < cycle_length; ++ii) {
    CHECK(ndata.get_degree(ii) == 2);
    CHECK(ndata.get_sorted_degree_sequence_expensive(ii) == deg_seqs);
    const auto neighbours = ndata.get_neighbours_expensive(ii);
    CHECK(neighbours.size() == 2);
    CHECK(is_sorted_and_unique(neighbours));

    neighbours_recalc[0] = ((ii + cycle_length) - 1) % cycle_length;
    neighbours_recalc[1] = (ii + 1) % cycle_length;
    std::sort(neighbours_recalc.begin(), neighbours_recalc.end());
    CHECK(neighbours_recalc == neighbours);
  }
  const auto vertices = ndata.get_nonisolated_vertices_expensive();
  CHECK(vertices.size() == cycle_length);
  CHECK(vertices[0] == 0);
  CHECK(vertices.back() + 1 == cycle_length);
  CHECK(is_sorted_and_unique(vertices));

  // Now, edge functions
  for (VertexWSM ii = 0; ii < cycle_length; ++ii) {
    for (VertexWSM jj = 0; jj < cycle_length; ++jj) {
      const auto edge_weight_opt = ndata.get_edge_weight_opt(ii, jj);
      // Is i-j == +/-1 (mod cycle length)??
      const auto diff = ((ii + cycle_length) - jj) % cycle_length;
      if (diff == 1 || diff + 1 == cycle_length) {
        CHECK(edge_weight_opt);
        CHECK(edge_weight_opt.value() == ii + jj);
      } else {
        CHECK(!edge_weight_opt);
      }
    }
  }
  // Nonexistent vertices
  for (VertexWSM ii = 0; ii < cycle_length + 5; ++ii) {
    for (VertexWSM jj = cycle_length; jj < cycle_length + 10; ++jj) {
      CHECK(!ndata.get_edge_weight_opt(ii, jj));
    }
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
