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
#include <numeric>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/GraphTheoretic/NeighboursData.hpp>
#include <tkwsm/InitPlacement/PrunedTargetEdges.hpp>

#include "TestWeightedGraphData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

static void test_validity_of_new_graph_data(
    const NeighboursData& pattern_ndata,
    const NeighboursData& explicit_target_ndata,
    WeightWSM implicit_target_weight,
    // Destructively modifies, i.e. sorts.
    std::vector<unsigned>& assigned_target_vertices,
    const GraphEdgeWeights& new_target_graph_data) {
  // Check that all weights and TV are valid.
  const unsigned number_of_tv =
      explicit_target_ndata.get_number_of_nonisolated_vertices();
  for (const auto& entry : new_target_graph_data) {
    const auto& edge = entry.first;
    REQUIRE(edge.first < edge.second);
    REQUIRE(edge.second < number_of_tv);
    const auto t_weight_opt =
        explicit_target_ndata.get_edge_weight_opt(edge.first, edge.second);
    if (t_weight_opt) {
      CHECK(entry.second == t_weight_opt.value());
    } else {
      CHECK(entry.second == implicit_target_weight);
    }
  }

  // Check that every original p-edge corresponds to a new t-edge.
  REQUIRE(
      assigned_target_vertices.size() ==
      pattern_ndata.get_number_of_nonisolated_vertices());
  for (unsigned pv = 0; pv < assigned_target_vertices.size(); ++pv) {
    const auto& p_neighbours_and_weights =
        pattern_ndata.get_neighbours_and_weights(pv);
    for (const auto& other_pv_weight_pair : p_neighbours_and_weights) {
      const auto& pv_other = other_pv_weight_pair.first;
      // INFO("pv=" << pv << ", joined to " << pv_other);
      REQUIRE(pv != pv_other);
      REQUIRE(pv_other < assigned_target_vertices.size());
      REQUIRE(
          new_target_graph_data.count(get_edge(
              assigned_target_vertices[pv],
              assigned_target_vertices[pv_other])) != 0);
    }
  }

  // Finally, check that every new t-edge has at least one originally used TV.
  // Here we cheat and sort the t-vertices list, for binary searches.
  std::sort(assigned_target_vertices.begin(), assigned_target_vertices.end());
  for (const auto& entry : new_target_graph_data) {
    const auto& edge = entry.first;
    if (!std::binary_search(
            assigned_target_vertices.cbegin(), assigned_target_vertices.cend(),
            edge.first)) {
      REQUIRE(std::binary_search(
          assigned_target_vertices.cbegin(), assigned_target_vertices.cend(),
          edge.second));
    }
  }
}

SCENARIO("Test adding unsued edges to simple random graphs with assignments") {
  TargetEdgePruningParameters parameters;
  parameters.max_additional_number_of_target_edges_factor_per_kilo = 800;
  parameters.min_implicit_unused_number_of_target_edges_factor_per_kilo = 400;

  // Just choose some "random" value; doesn't actually matter.
  const WeightWSM implicit_target_weight = 9999;

  std::vector<unsigned> assigned_target_vertices;
  RNG rng;

  {
    // Small test.
    const auto pattern_graph_data = get_graph_data(rng, 8, 20, 1000, 2000);
    const NeighboursData pattern_ndata(pattern_graph_data);
    {
      // Same number of target vertices. But be sneaky
      // and use fewer edges, thus FORCING the implicit completeness
      // to be used.
      const auto explicit_target_graph_data =
          get_graph_data(rng, 8, 10, 10, 100);

      const NeighboursData explicit_target_ndata(explicit_target_graph_data);
      assigned_target_vertices = {3, 2, 0, 4, 1, 6, 5, 7};
      const auto new_target_graph_data = get_new_target_graph_data(
          pattern_ndata, explicit_target_ndata, implicit_target_weight,
          assigned_target_vertices, parameters);

      CHECK(pattern_graph_data.size() == 15);
      CHECK(explicit_target_graph_data.size() == 10);
      CHECK(new_target_graph_data.size() == 18);

      test_validity_of_new_graph_data(
          pattern_ndata, explicit_target_ndata, implicit_target_weight,
          assigned_target_vertices, new_target_graph_data);
    }
    // Add some more target edges.
    {
      const auto explicit_target_graph_data =
          get_graph_data(rng, 8, 20, 10, 100);
      const NeighboursData explicit_target_ndata(explicit_target_graph_data);
      assigned_target_vertices = {2, 1, 5, 3, 7, 0, 4, 6};
      const auto new_target_graph_data = get_new_target_graph_data(
          pattern_ndata, explicit_target_ndata, implicit_target_weight,
          assigned_target_vertices, parameters);

      CHECK(explicit_target_graph_data.size() == 15);
      CHECK(new_target_graph_data.size() == 18);

      test_validity_of_new_graph_data(
          pattern_ndata, explicit_target_ndata, implicit_target_weight,
          assigned_target_vertices, new_target_graph_data);
    }
    // Add more target vertices.
    {
      const auto explicit_target_graph_data =
          get_graph_data(rng, 15, 50, 10, 100);
      const NeighboursData explicit_target_ndata(explicit_target_graph_data);
      assigned_target_vertices = {7, 3, 8, 0, 9, 2, 6, 13};
      const auto new_target_graph_data = get_new_target_graph_data(
          pattern_ndata, explicit_target_ndata, implicit_target_weight,
          assigned_target_vertices, parameters);

      CHECK(explicit_target_graph_data.size() == 37);
      CHECK(new_target_graph_data.size() == 26);

      test_validity_of_new_graph_data(
          pattern_ndata, explicit_target_ndata, implicit_target_weight,
          assigned_target_vertices, new_target_graph_data);
    }
  }

  // Finally, a bigger test.
  const auto pattern_graph_data = get_graph_data(rng, 20, 80, 1, 1000);

  const NeighboursData pattern_ndata(pattern_graph_data);

  // Notice that the target graph actually has quite low edge density,
  // lower than the pattern graph. This is probably quite common with real
  // applications:
  // connectivity in many real quantum computers (e.g. the IBM heavy hexagon
  // "brick wall pattern" machines) is quite low;
  // qubit intereactions in a real quantum circuit are presumably often
  // quite numerous, as you often have many gates and a lot of entanglement
  // in a useful circuit.
  const auto explicit_target_graph_data =
      get_graph_data(rng, 40, 150, 100, 500);

  const NeighboursData explicit_target_ndata(explicit_target_graph_data);
  assigned_target_vertices.resize(40);
  std::iota(
      assigned_target_vertices.begin(), assigned_target_vertices.end(), 0);
  rng.do_shuffle(assigned_target_vertices);

  // Cut down to size.
  assigned_target_vertices.resize(20);

  const auto new_target_graph_data = get_new_target_graph_data(
      pattern_ndata, explicit_target_ndata, implicit_target_weight,
      assigned_target_vertices, parameters);

  CHECK(pattern_graph_data.size() == 61);
  CHECK(explicit_target_graph_data.size() == 135);
  CHECK(new_target_graph_data.size() == 108);

  test_validity_of_new_graph_data(
      pattern_ndata, explicit_target_ndata, implicit_target_weight,
      assigned_target_vertices, new_target_graph_data);
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
