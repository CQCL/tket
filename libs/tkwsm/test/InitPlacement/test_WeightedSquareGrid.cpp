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

#include <catch2/catch_test_macros.hpp>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

#include "TestUtilsIQP.hpp"
#include "WeightedSquareGrid.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

SCENARIO("small fixed square grid") {
  // For a square grid of 9 points, there are 12 edges.
  std::vector<WeightWSM> weights(12);
  for (unsigned ii = 0; ii < weights.size(); ++ii) {
    weights[ii] = (ii + 1) * 100;
  }
  // In a "random" order; but specifically include some horizontal,
  // vertical, and adjacent pairs.
  // The rows are [6 7 8], [3 4 5], [0 1 2].
  const std::vector<VertexWSM> vertices{5, 2, 8, 4, 0, 1, 6, 7, 3, 0, 8};

  // The vertices and weights layout should be
  //
  //    6 --500-- 7 --600-- 8
  //    |         |         |
  //   800      1000      1200
  //    |         |         |
  //    3 --300-- 4 --400-- 5
  //    |         |         |
  //   700       900      1100
  //    |         |         |
  //    0 --100-- 1 --200-- 2
  //
  WeightedSquareGrid grid(weights);
  CHECK(
      str(grid.get_graph_data()) ==
      "12 edges with weights: [  (0,1: 100),  (0,3: 700),  (1,2: 200),  "
      "(1,4: 900),  (2,5: 1100),  (3,4: 300),  (3,6: 800),  (4,5: 400),  "
      "(4,7: 1000),  (5,8: 1200),  (6,7: 500),  (7,8: 600), ]\n"
      "9 vertices: {0 1 2 3 4 5 6 7 8 }\n");

  // clang-format off
  // These have been manually checked! But paths are, of course,
  // not necessarily optimal.
  CHECK(get_many_paths_test_str(grid, vertices) ==
    "\nPath vertices: [ 5 2 ]"
    "\nEdge weights: [ 1100 ] (total weight 1100)"
    "\nPath vertices: [ 2 5 8 ]"
    "\nEdge weights: [ 1100 1200 ] (total weight 2300)"
    "\nPath vertices: [ 8 5 4 ]"
    "\nEdge weights: [ 1200 400 ] (total weight 1600)"
    "\nPath vertices: [ 4 1 0 ]"
    "\nEdge weights: [ 900 100 ] (total weight 1000)"
    "\nPath vertices: [ 0 1 ]"
    "\nEdge weights: [ 100 ] (total weight 100)"
    "\nPath vertices: [ 1 0 3 6 ]"
    "\nEdge weights: [ 100 700 800 ] (total weight 1600)"
    "\nPath vertices: [ 6 7 ]"
    "\nEdge weights: [ 500 ] (total weight 500)"
    "\nPath vertices: [ 7 4 3 ]"
    "\nEdge weights: [ 1000 300 ] (total weight 1300)"
    "\nPath vertices: [ 3 0 ]"
    "\nEdge weights: [ 700 ] (total weight 700)"
    "\nPath vertices: [ 0 1 2 5 8 ]"
    "\nEdge weights: [ 100 200 1100 1200 ] (total weight 2600)"
  );
  // clang-format on

  // Now try some token swaps. Let PV be multiples of 11.
  // Not a full placement!
  std::vector<std::pair<VertexWSM, VertexWSM>> initial_placement;
  for (unsigned ii = 2; ii <= 8; ++ii) {
    initial_placement.emplace_back(ii * 11, ii);
  }
  grid.initialise_with_qubit_placement(initial_placement);

  // Manually checked!
  CHECK(
      do_token_swaps_and_check_placements(
          std::vector<std::pair<VertexWSM, VertexWSM>>{
              {22, 66}, {88, 33}, {55, 77}},
          grid) ==
      // clang-format off
      "\n"
      "TOKEN Swap (22,66) between vertices 2 6; cost 3800\n"
      "Path vertices: [ 2 1 0 3 6 ]\n"
      "Edge weights: [ 200 100 700 800 ] (total weight 1800)\n"
      "NOW, placement: { 22->3 33->0 44->4 55->5 66->6 77->7 88->8 }\n"
      "\n"
      "TOKEN Swap (88,33) between vertices 8 0; cost 5400\n"
      "Path vertices: [ 8 5 2 1 0 ]\n"
      "Edge weights: [ 1200 1100 200 100 ] (total weight 2600)\n"
      "NOW, placement: { 22->3 33->5 44->4 55->2 66->6 77->7 88->8 }\n"
      "\n"
      "TOKEN Swap (55,77) between vertices 2 7; cost 4300\n"
      "Path vertices: [ 2 1 4 7 ]\n"
      "Edge weights: [ 200 900 1000 ] (total weight 2100)\n"
      "NOW, placement: { 22->3 33->5 44->1 55->4 66->6 77->7 88->8 }\n"
      // clang-format on
  );
}

SCENARIO("10x10 square grid") {
  // Do a TRICK with the weights: first half (horizontal weights)
  // are all 0 (mod 10); second half (vertical weights) 1 (mod 10).
  std::vector<WeightWSM> weights(180);
  RNG rng;
  for (unsigned ii = 0; ii < weights.size(); ++ii) {
    weights[ii] = 10 * (1 + rng.get_size_t(5));
    if (ii >= 90) {
      ++weights[ii];
    }
  }
  // In a "random" order.
  std::vector<VertexWSM> vertices{0};
  while (vertices.size() < 10) {
    for (;;) {
      const auto next_v = rng.get_size_t(99);
      if (next_v != vertices.back()) {
        vertices.push_back(next_v);
        break;
      }
    }
  }
  const WeightedSquareGrid grid(weights);

  // clang-format off
  // To go to the right, add 1; to go up, add 10.
  // Bottom row of grid is [0 1 2 ... 9].
  // Mostly manually checked...
  CHECK(get_many_paths_test_str(grid, vertices) ==
    "\nPath vertices: [ 0 10 11 12 13 14 15 25 35 45 55 65 75 85 ]"
    "\nEdge weights: [ 51 30 20 40 10 40 11 11 11 11 41 21 61 ] (total weight 358)"
    "\nPath vertices: [ 85 84 74 ]"
    "\nEdge weights: [ 60 11 ] (total weight 71)"
    "\nPath vertices: [ 74 64 65 55 45 35 36 ]"
    "\nEdge weights: [ 21 10 41 11 11 50 ] (total weight 144)"
    "\nPath vertices: [ 36 35 45 55 65 64 ]"
    "\nEdge weights: [ 50 11 11 41 10 ] (total weight 123)"
    "\nPath vertices: [ 64 65 66 67 57 58 ]"
    "\nEdge weights: [ 10 30 10 31 20 ] (total weight 101)"
    "\nPath vertices: [ 58 57 67 66 65 64 74 84 ]"
    "\nEdge weights: [ 20 31 10 30 10 21 11 ] (total weight 133)"
    "\nPath vertices: [ 84 74 64 54 44 43 42 41 31 21 11 10 0 ]"
    "\nEdge weights: [ 11 21 31 41 10 20 40 11 31 31 30 51 ] (total weight 328)"
    "\nPath vertices: [ 0 1 2 3 4 5 6 7 8 18 ]"
    "\nEdge weights: [ 50 20 50 60 10 30 20 10 11 ] (total weight 261)"
    "\nPath vertices: [ 18 19 29 ]"
    "\nEdge weights: [ 30 31 ] (total weight 61)"
  );
  // clang-format on
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
