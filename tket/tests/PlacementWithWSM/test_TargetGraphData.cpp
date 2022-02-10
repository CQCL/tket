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

#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <sstream>

#include "../WeightSubgrMono/TestUtils/FixedArchitectures.hpp"
#include "../WeightSubgrMono/TestUtils/GraphGeneration.hpp"
#include "PlacementWithWSM/TargetGraphData.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

SCENARIO("basic cycle") {
  TargetGraphData::Parameters parameters;
  const auto graph = GraphGeneration::get_cycle(5);
  CHECK(
      str(graph) ==
      "5 edges with weights: [  (0,1: 1),  (0,4: 5),  (1,2: 2),  (2,3: 3),  "
      "(3,4: 4), ]\n"
      "5 vertices: {0 1 2 3 4 }\n");

  const TargetGraphData added_edges(graph, parameters);
  CHECK(
      str(added_edges.final_data) ==
      "10 edges with weights: [  (0,1: 1),  (0,2: 5),  (0,3: 17),  (0,4: 5),  "
      "(1,2: 2),"
      "  (1,3: 9),  (1,4: 8),  (2,3: 3),  (2,4: 13),  (3,4: 4), ]\n"
      "5 vertices: {0 1 2 3 4 }\n");

  parameters.max_edge_weight = 10;
  parameters.remove_high_edge_weights = false;
  {
    const TargetGraphData added_edges_capped_weights(graph, parameters);
    CHECK(
        added_edges_capped_weights.final_data.size() ==
        added_edges.final_data.size());
    for (const auto& entry : added_edges_capped_weights.final_data) {
      const auto uncapped_weight = added_edges.final_data.at(entry.first);
      const auto capped_weight =
          std::min(uncapped_weight, parameters.max_edge_weight.value());
      CHECK(entry.second == capped_weight);
    }
  }
  parameters.remove_high_edge_weights = true;
  {
    const TargetGraphData erased_high_weights(graph, parameters);
    for (const auto& entry : erased_high_weights.final_data) {
      CHECK(entry.second == added_edges.final_data.at(entry.first));
      CHECK(entry.second <= parameters.max_edge_weight.value());
    }
    // Everything <= the max should be present.
    for (const auto& entry : added_edges.final_data) {
      if (entry.second <= parameters.max_edge_weight.value()) {
        CHECK(erased_high_weights.final_data.count(entry.first) == 1);
      }
    }
  }
}

SCENARIO("larger fixed architectures") {
  const auto graph = FixedArchitectures::get_ibm_guadalupe_16_qubits();
  TargetGraphData::Parameters parameters;
  CHECK(
      str(graph) ==
      "16 edges with weights: [  (0,1: 1),  (1,2: 1),  "
      "(1,4: 1),  (2,3: 1),  (3,5: 1),  (4,7: 1),  (5,8: 1),  (6,7: 1),  "
      "(7,10: 1),  (8,9: 1),  (8,11: 1),  (10,12: 1),  (11,14: 1),  "
      "(12,13: 1),  (12,15: 1),  (13,14: 1), ]\n"
      "16 vertices: {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 }\n");

  const TargetGraphData added_edges(graph, parameters);
  CHECK(
      str(added_edges.final_data) ==
      "90 edges with weights: [  (0,1: 1),"
      "  (0,2: 4),  (0,3: 7),  (0,4: 4),  (0,5: 10),  (0,6: 10),  (0,7: 7),"
      "  (0,8: 13),  (0,10: 10),  (0,12: 13),  (1,2: 1),  (1,3: 4),  (1,4: 1),"
      "  (1,5: 7),  (1,6: 7),  (1,7: 4),  (1,8: 10),  (1,9: 13),  (1,10: 7),"
      "  (1,11: 13),  (1,12: 10),  (1,13: 13),  (1,15: 13),  (2,3: 1),"
      "  (2,4: 4),  (2,5: 4),  (2,6: 10),  (2,7: 7),  (2,8: 7),  (2,9: 10),"
      "  (2,10: 10),  (2,11: 10),  (2,12: 13),  (2,14: 13),  (3,5: 1),"
      "  (3,8: 4),  (3,9: 7),  (3,11: 7),  (3,13: 13),  (3,14: 10),"
      "  (4,6: 4),  (4,7: 1),  (4,10: 4),  (4,12: 7),  (4,13: 10),"
      "  (4,14: 13),  (4,15: 10),  (5,8: 1),  (5,9: 4),  (5,11: 4),"
      "  (5,12: 13),  (5,13: 10),  (5,14: 7),  (6,7: 1),  (6,10: 4),"
      "  (6,12: 7),  (6,13: 10),  (6,14: 13),  (6,15: 10),  (7,10: 1),"
      "  (7,11: 13),  (7,12: 4),  (7,13: 7),  (7,14: 10),  (7,15: 7),"
      "  (8,9: 1),  (8,10: 13),  (8,11: 1),  (8,12: 10),  (8,13: 7),"
      "  (8,14: 4),  (8,15: 13),  (9,11: 4),  (9,12: 13),  (9,13: 10),"
      "  (9,14: 7),  (10,11: 10),  (10,12: 1),  (10,13: 4),  (10,14: 7),"
      "  (10,15: 4),  (11,12: 7),  (11,13: 4),  (11,14: 1),  (11,15: 10),"
      "  (12,13: 1),  (12,14: 4),  (12,15: 1),  (13,14: 1),  (13,15: 4), ]\n"
      "16 vertices: {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 }\n");
}

}  // namespace tket
