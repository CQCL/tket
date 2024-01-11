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

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/Common/TemporaryRefactorCode.hpp>
#include <tkwsm/GraphTheoretic/NearNeighboursData.hpp>
#include <tkwsm/GraphTheoretic/NeighboursData.hpp>

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
  const auto number_of_vertices = ndata.get_number_of_nonisolated_vertices();
  CHECK(number_of_vertices == cycle_length);

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

  // Also, test near neighbours.
  NearNeighboursData near_neighbours_data(ndata);
  for (VertexWSM ii = 0; ii < cycle_length; ++ii) {
    for (unsigned distance = 0; distance <= 8; ++distance) {
      const auto count_within_d =
          near_neighbours_data.get_n_vertices_up_to_distance(ii, distance);
      switch (distance) {
        case 0:
          CHECK(count_within_d == 0);
          break;
        case 1:
          CHECK(count_within_d == 2);
          break;
        default:
          CHECK(count_within_d == 4);
          break;
      }
      if (distance < 2) {
        continue;
      }
      const auto& v_at_distance =
          near_neighbours_data.get_vertices_at_exact_distance(ii, distance);
      REQUIRE(v_at_distance.size() == number_of_vertices);

      if (distance > 2) {
        CHECK(v_at_distance.none());
        continue;
      }
      std::set<VertexWSM> v_set_at_distance;
      TemporaryRefactorCode::set_domain_from_bitset(
          v_set_at_distance, v_at_distance);
      CHECK(v_set_at_distance.size() == 2);
      for (auto v : v_set_at_distance) {
        const unsigned route1_dist = ((cycle_length + ii) - v) % cycle_length;
        const unsigned route2_dist =
            (cycle_length - route1_dist) % cycle_length;
        const unsigned shortest_distance = std::min(route1_dist, route2_dist);
        CHECK(shortest_distance == 2);
      }
    }
  }
}

static std::string to_string(const NeighboursData& ndata) {
  std::stringstream ss;
  const auto number_of_vertices = ndata.get_number_of_nonisolated_vertices();
  ss << number_of_vertices << " vertices. Neighbours and weights:";
  for (unsigned vv = 0; vv < number_of_vertices; ++vv) {
    ss << "\nv=" << vv << ": [ ";
    const std::vector<std::pair<VertexWSM, WeightWSM>>& data =
        ndata.get_neighbours_and_weights(vv);
    for (const std::pair<VertexWSM, WeightWSM>& entry : data) {
      ss << entry.first << ";" << entry.second << " ";
    }
    ss << "]";
  }
  return ss.str();
}

SCENARIO("neighbours_data with invalid and simple input data") {
  GraphEdgeWeights edge_weights;
  CHECK_THROWS_AS(NeighboursData(edge_weights), std::runtime_error);
  edge_weights[std::make_pair<VertexWSM, VertexWSM>(0, 0)] = 1;
  CHECK_THROWS_AS(NeighboursData(edge_weights), std::runtime_error);
  edge_weights.clear();
  edge_weights[get_edge(0, 1)] = 1;
  // v1>v2 is allowed...
  edge_weights[std::make_pair<VertexWSM, VertexWSM>(2, 0)] = 2;
  const NeighboursData ndata1(edge_weights);
  const auto ndata1_str = to_string(ndata1);
  CHECK(
      ndata1_str ==
      "3 vertices. Neighbours and weights:"
      "\nv=0: [ 1;1 2;2 ]"
      "\nv=1: [ 0;1 ]"
      "\nv=2: [ 0;2 ]");

  // Inconsistent edge weights are not allowed...
  edge_weights[std::make_pair<VertexWSM, VertexWSM>(0, 2)] = 3;
  REQUIRE(edge_weights.size() == 3);
  CHECK_THROWS_AS(NeighboursData(edge_weights), std::runtime_error);

  //...but duplicate data IS allowed, as long as it's not inconsistent
  edge_weights[std::make_pair<VertexWSM, VertexWSM>(0, 2)] = 2;
  const NeighboursData ndata2(edge_weights);
  CHECK(ndata1_str == to_string(ndata2));
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
