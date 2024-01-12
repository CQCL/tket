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
#include <tkrng/RNG.hpp>
#include <tkwsm/InitPlacement/InputStructs.hpp>

#include "TestUtilsIQP.hpp"
#include "WeightedBinaryTree.hpp"
#include "WeightedSquareGrid.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

// These tests also record some results found with try_random_placements
// previously. (The code might have changed since these results,
// but I haven't bothered to update them.
// No point really; there's still plenty of experimentation to be done.
// This is more of a proof-of-concept that the WSM approach CAN give good
// solutions quickly).
//
// Results with binary trees are relatively faster, compared with IQP,
// than for square grids; but this is not surprising.
// There are several possibilities.
//
// (1) Path finding in binary trees is simpler than for square grids,
//    so maybe try_random_placements just churns through more choices.
//
// (2) Maybe the IQP from WSM approach is just inherently
//    not as good for trees as it is for square grids.
//
// (3) Maybe the MCCT parameters are just not tweaked very well.
//
// Results for square grids seem very good; for binary trees, not so good.
// But there's a lot of extra testing and experimentation to be done!
// Finding good default parameters, and maybe a better strategy for
// reintroducing deleted target edges, should give much better results.
//
// Not enough data for a clear conclusion; need more experiments!
//

SCENARIO("IQP for 10x10 square grid; 30 and 90 logical qubits") {
  std::vector<WeightWSM> weights;
  weights.reserve(180);
  RNG rng;
  while (weights.size() < 180) {
    weights.push_back(2 + rng.get_size_t(4));
  }
  WeightedSquareGrid square_grid(weights, 3);

  // Try two different placement problems on the same grid.
  // Of course, the actual number of gates will be slightly lower
  // (89 and 465, respecitvely)
  // due to random v1==v2 collisions which are simply skipped.
  std::vector<std::pair<VertexWSM, VertexWSM>> gates_30_pv_100_gates;
  std::vector<std::pair<VertexWSM, VertexWSM>> gates_90_pv_500_gates;
  {
    const unsigned num_pv = 30;
    const unsigned num_gates = 100;
    gates_30_pv_100_gates.reserve(num_gates);
    while (gates_30_pv_100_gates.size() < num_gates) {
      for (;;) {
        const VertexWSM first_v = rng.get_size_t(num_pv - 1);
        const VertexWSM second_v = rng.get_size_t(num_pv - 1);
        if (first_v != second_v) {
          // Make non contiguous PV.
          gates_30_pv_100_gates.emplace_back(
              10 + 2 * first_v, 10 + 2 * second_v);
          break;
        }
      }
    }
  }
  {
    const unsigned num_pv = 90;
    const unsigned num_gates = 500;
    gates_90_pv_500_gates.reserve(num_gates);
    while (gates_90_pv_500_gates.size() < num_gates) {
      for (;;) {
        const VertexWSM first_v = rng.get_size_t(num_pv - 1);
        const VertexWSM second_v = rng.get_size_t(num_pv - 1);
        if (first_v != second_v) {
          // Make non contiguous PV.
          gates_90_pv_500_gates.emplace_back(100 + first_v, 100 + second_v);
          break;
        }
      }
    }
  }
  // Uncomment to print out some results!
  // square_grid.try_random_placements(gates_30_pv_100_gates);
  // square_grid.try_random_placements(gates_90_pv_500_gates);

  // We must construct the pattern graph.
  // Obviously a lot more testing and experimentation should be done.
  PatternGraphDataInput pgd_input;
  pgd_input.initial_gate_weight = 1000;
  pgd_input.final_gate_weight = 10;
  pgd_input.method = PatternGraphDataInput::ReorderingMethod::ORIGINAL_ORDER;

  // Note that the MCCT solution (found in ~10ms) is optimal for the
  // internally generated WSM problem, but it takes ~13 seconds for
  // the main WSM routine to PROVE optimality.
  // Better than ANY solution found with try_random_placements!
  const PatternGraphData p_graph_data_30_pv(gates_30_pv_100_gates, pgd_input);
  run_end_to_end_IQP_and_check_solution(
      gates_30_pv_100_gates, p_graph_data_30_pv.pattern_graph_weights,
      square_grid, 1665);

  // Takes ~30 ms with MCCT. But try_random_placements needs >1 minute,
  // i.e. >1000x longer, to improve on it.
  const PatternGraphData p_graph_data_90_pv(gates_90_pv_500_gates, pgd_input);
  run_end_to_end_IQP_and_check_solution(
      gates_90_pv_500_gates, p_graph_data_90_pv.pattern_graph_weights,
      square_grid, 24526);

  // Now give some actual solutions found with try_random_placements.
  // Note that the above WSM-based solutions only take ~50 ms for the
  // MCCT part; in these cases the main WSM doesn't add anything
  // because the WSM problem set up by MCCT (by adding extra target edges)
  // does not seem to have a better solution than the MCCT one.
  // HOWEVER more experiments need to be done; maybe the strategy
  // for adding back target edges can be changed to give better WSM solutions.
  std::vector<CostedIQPSolution> costed_30_pv_solutions(2);

  // clang-format off
  // This took ~0.4 seconds to find with try_random_placements.
  costed_30_pv_solutions[0] = {3109, {{10,37}, {12,15}, {14,43}, {16,14},
    {18,44}, {20,5}, {22,21}, {24,29}, {26,55}, {28,30}, {30,28}, {32,27},
    {34,13}, {36,32}, {38,81}, {40,3}, {42,35}, {44,40}, {46,16}, {48,22},
    {50,12}, {52,33}, {54,45}, {56,6}, {58,58}, {60,38}, {62,17}, {64,11},
    {66,49}, {68,61}}};

  // Took ~23 seconds to find with try_random_placements.
  // This is still not as good as the IQP solution,
  // taking only ~50ms for the MCCT solution!
  costed_30_pv_solutions[1] = {1743, {{10,46}, {12,41}, {14,33}, {16,36},
    {18,42}, {20,56}, {22,11}, {24,16}, {26,5}, {28,4}, {30,6}, {32,25},
    {34,13}, {36,43}, {38,28}, {40,3}, {42,45}, {44,47}, {46,22}, {48,24},
    {50,34}, {52,32}, {54,37}, {56,14}, {58,23}, {60,38}, {62,27}, {64,35},
    {66,15}, {68,31}}};
  // clang-format on

  test_known_solutions(
      costed_30_pv_solutions, gates_30_pv_100_gates, square_grid);

  // Solutions for the 90 PV case. So, IQP currently can place
  // 90 logical qubits into 100 physical qubits, pretty well,
  // in only ~50 ms using MCCT.
  std::vector<CostedIQPSolution> costed_90_pv_solutions(5);

  // clang-format off
  // Took ~1 second to find with try_random_placements.
  costed_90_pv_solutions[0] = {26361, {{100,37}, {101,15}, {102,43}, {103,14},
    {104,44}, {105,5}, {106,1}, {107,29}, {108,55}, {109,35}, {110,80},
    {111,27}, {112,13}, {113,32}, {114,56}, {115,3}, {116,30}, {117,40},
    {118,16}, {119,22}, {120,20}, {121,33}, {122,2}, {123,79}, {124,58},
    {125,38}, {126,17}, {127,11}, {128,47}, {129,73}, {130,46}, {131,51},
    {132,21}, {133,61}, {134,18}, {135,0}, {136,49}, {137,6}, {138,75},
    {139,39}, {140,53}, {141,19}, {142,62}, {143,74}, {144,60}, {145,93},
    {146,87}, {147,8}, {148,59}, {149,89}, {150,78}, {151,96}, {152,25},
    {153,10}, {154,9}, {155,76}, {156,70}, {157,72}, {158,28}, {159,99},
    {160,65}, {161,77}, {162,91}, {163,63}, {164,50}, {165,66}, {166,82},
    {167,88}, {168,83}, {169,7}, {170,12}, {171,31}, {172,41}, {173,90},
    {174,54}, {175,67}, {176,57}, {177,85}, {178,42}, {179,48}, {180,4},
    {181,94}, {182,81}, {183,92}, {184,23}, {185,34}, {186,95}, {187,26},
    {188,69}, {189,24}}};

  // Took ~2 minutes.
  costed_90_pv_solutions[1] = {24091, {{100,28}, {101,25}, {102,74}, {103,41},
    {104,62}, {105,26}, {106,0}, {107,18}, {108,55}, {109,29}, {110,80},
    {111,16}, {112,13}, {113,34}, {114,64}, {115,3}, {116,71}, {117,32},
    {118,27}, {119,14}, {120,76}, {121,33}, {122,2}, {123,59}, {124,58},
    {125,38}, {126,17}, {127,11}, {128,36}, {129,72}, {130,1}, {131,61},
    {132,39}, {133,5}, {134,8}, {135,81}, {136,48}, {137,6}, {138,75},
    {139,21}, {140,53}, {141,19}, {142,99}, {143,40}, {144,60}, {145,91},
    {146,63}, {147,35}, {148,78}, {149,89}, {150,79}, {151,95}, {152,22},
    {153,10}, {154,9}, {155,94}, {156,70}, {157,73}, {158,37}, {159,86},
    {160,67}, {161,77}, {162,93}, {163,84}, {164,23}, {165,87}, {166,66},
    {167,98}, {168,83}, {169,7}, {170,12}, {171,51}, {172,20}, {173,90},
    {174,54}, {175,85}, {176,57}, {177,56}, {178,52}, {179,49}, {180,4},
    {181,44}, {182,46}, {183,92}, {184,31}, {185,50}, {186,96}, {187,45},
    {188,68}, {189,24}}};

  // Took ~3 minutes.
  costed_90_pv_solutions[2] = {23053, {{100,57}, {101,25}, {102,74}, {103,41},
    {104,62}, {105,35}, {106,0}, {107,18}, {108,26}, {109,29}, {110,80},
    {111,16}, {112,13}, {113,34}, {114,77}, {115,3}, {116,42}, {117,32},
    {118,27}, {119,14}, {120,87}, {121,33}, {122,2}, {123,59}, {124,58},
    {125,38}, {126,17}, {127,11}, {128,36}, {129,61}, {130,1}, {131,73},
    {132,39}, {133,5}, {134,8}, {135,93}, {136,48}, {137,6}, {138,75},
    {139,30}, {140,53}, {141,19}, {142,88}, {143,51}, {144,60}, {145,91},
    {146,63}, {147,46}, {148,78}, {149,89}, {150,69}, {151,95}, {152,22},
    {153,10}, {154,9}, {155,82}, {156,70}, {157,71}, {158,37}, {159,86},
    {160,67}, {161,76}, {162,81}, {163,96}, {164,23}, {165,66}, {166,64},
    {167,98}, {168,83}, {169,7}, {170,12}, {171,40}, {172,20}, {173,90},
    {174,65}, {175,85}, {176,56}, {177,28}, {178,52}, {179,49}, {180,4},
    {181,44}, {182,45}, {183,92}, {184,31}, {185,50}, {186,84}, {187,55},
    {188,79}, {189,24}}};

  // Took ~10 minutes.
  costed_90_pv_solutions[3] = {22376, {{100,46}, {101,75}, {102,73}, {103,42},
    {104,63}, {105,28}, {106,0}, {107,18}, {108,25}, {109,8}, {110,80},
    {111,67}, {112,13}, {113,33}, {114,76}, {115,3}, {116,41}, {117,32},
    {118,5}, {119,15}, {120,87}, {121,31}, {122,2}, {123,59}, {124,27},
    {125,37}, {126,17}, {127,11}, {128,55}, {129,72}, {130,1}, {131,77},
    {132,29}, {133,4}, {134,19}, {135,84}, {136,48}, {137,6}, {138,74},
    {139,30}, {140,53}, {141,9}, {142,88}, {143,51}, {144,60}, {145,61},
    {146,62}, {147,45}, {148,78}, {149,89}, {150,69}, {151,95}, {152,22},
    {153,10}, {154,39}, {155,82}, {156,70}, {157,65}, {158,26}, {159,86},
    {160,56}, {161,54}, {162,81}, {163,96}, {164,23}, {165,66}, {166,64},
    {167,98}, {168,94}, {169,7}, {170,12}, {171,40}, {172,20}, {173,90},
    {174,58}, {175,85}, {176,16}, {177,38}, {178,52}, {179,49}, {180,36},
    {181,44}, {182,35}, {183,92}, {184,43}, {185,50}, {186,83}, {187,57},
    {188,79}, {189,24}}};

  // A previous IQP solution found in ~50 ms with MCCT
  // by tweaking input parameters.
  // (Better than the current solution - but it needs extensive tweaking
  // and experimentation to find the best default parameters - or maybe
  // calculate good parameters from input data).
  costed_90_pv_solutions[4] = {23578, {{100,60}, {101,35}, {102,55}, {103,9},
    {104,78}, {105,40}, {106,33}, {107,83}, {108,73}, {109,84}, {110,49},
    {111,43}, {112,31}, {113,25}, {114,13}, {115,58}, {116,7}, {117,41},
    {118,15}, {119,52}, {120,67}, {121,79}, {122,62}, {123,50}, {124,91},
    {125,21}, {126,70}, {127,97}, {128,56}, {129,5}, {130,82}, {131,94},
    {132,27}, {133,23}, {134,4}, {135,34}, {136,17}, {137,32}, {138,37},
    {139,39}, {140,14}, {141,53}, {142,92}, {143,98}, {144,29}, {145,61},
    {146,28}, {147,63}, {148,93}, {149,45}, {150,46}, {151,36}, {152,64},
    {153,85}, {154,30}, {155,51}, {156,66}, {157,89}, {158,20}, {159,44},
    {160,96}, {161,95}, {162,8}, {163,6}, {164,22}, {165,59}, {166,88},
    {167,57}, {168,80}, {169,86}, {170,72}, {171,42}, {172,68}, {173,77},
    {174,54}, {175,16}, {176,69}, {177,18}, {178,75}, {179,74}, {180,24},
    {181,65}, {182,38}, {183,81}, {184,76}, {185,48}, {186,71}, {187,26},
    {188,47}, {189,87}}};
  // clang-format on

  test_known_solutions(
      costed_90_pv_solutions, gates_90_pv_500_gates, square_grid);
}

SCENARIO("Binary tree with ~30 vertices; almost full embedding") {
  RNG rng;
  std::vector<WeightWSM> weights{0, 0};
  weights.resize(30);
  for (unsigned nn = 2; nn < weights.size(); ++nn) {
    weights[nn] = 1 + rng.get_size_t(5);
  }
  WeightedBinaryTree tree(weights, 3);

  std::vector<std::pair<VertexWSM, VertexWSM>> gates(50);

  // Vertex 0 doesn't exist in the binary tree; only {1,2,3,...}
  const unsigned max_token = tree.get_max_vertex_number() - 1;

  for (unsigned ii = 0; ii < gates.size();) {
    gates[ii].first = rng.get_size_t(max_token);
    gates[ii].second = rng.get_size_t(max_token);
    if (gates[ii].first == gates[ii].second) {
      continue;
    }
    ++ii;
  }
  CHECK(gates.size() == 50);

  // Uncomment to print out some results!
  // tree.try_random_placements(gates);

  // Solve with IQP.
  PatternGraphDataInput pgd_input;
  pgd_input.initial_gate_weight = 100;
  pgd_input.final_gate_weight = 20;
  pgd_input.method = PatternGraphDataInput::ReorderingMethod::ORIGINAL_ORDER;
  const PatternGraphData p_graph_data(gates, pgd_input);
  CHECK(get_number_of_vertices(p_graph_data.pattern_graph_weights) == 27);

  run_end_to_end_IQP_and_check_solution(
      gates, p_graph_data.pattern_graph_weights, tree, 1157);

  std::vector<CostedIQPSolution> costed_solutions(4);

  // clang-format off
  // This took ~0.4 seconds to find with try_random_placements;
  costed_solutions[0] = {1131, {{0,7}, {1,11}, {2,6}, {3,17}, {4,26}, {6,8},
    {7,5}, {8,20}, {9,9}, {10,21}, {12,4}, {13,2}, {14,24}, {15,19}, {16,15},
    {17,1}, {18,23}, {19,27}, {20,18}, {21,10}, {22,3}, {23,13}, {24,25},
    {25,12}, {26,22}, {27,16}, {28,14}}};

  // ~0.5 seconds.
  costed_solutions[1] = {1068, {{0,12}, {1,11}, {2,7}, {3,17}, {4,6}, {6,4},
    {7,5}, {8,20}, {9,18}, {10,23}, {12,8}, {13,2}, {14,24}, {15,9}, {16,15},
    {17,1}, {18,10}, {19,19}, {20,13}, {21,21}, {22,3}, {23,27}, {24,25},
    {25,26}, {26,22}, {27,16}, {28,14}}};

  // Took ~1 second to find.
  costed_solutions[2] = {874, {{0,6}, {1,20}, {2,3}, {3,18}, {4,13}, {6,4},
    {7,5}, {8,22}, {9,17}, {10,23}, {12,8}, {13,11}, {14,12}, {15,9}, {16,15},
    {17,1}, {18,10}, {19,25}, {20,24}, {21,21}, {22,7}, {23,27}, {24,19},
    {25,26}, {26,2}, {27,16}, {28,14}}};

  // This took ~2.5 seconds  and
  // it didn't improve after >1 minute.
  costed_solutions[3] = {722, {{0,6}, {1,10}, {2,3}, {3,24}, {4,13}, {6,4},
    {7,5}, {8,23}, {9,16}, {10,22}, {12,9}, {13,11}, {14,14}, {15,8}, {16,15},
    {17,1}, {18,21}, {19,26}, {20,28}, {21,20}, {22,7}, {23,27}, {24,19},
    {25,25}, {26,2}, {27,18}, {28,12}}};
  // clang-format on
  test_known_solutions(costed_solutions, gates, tree);
}

SCENARIO("Binary tree with ~100 vertices, ~30 logical qubits") {
  RNG rng;
  std::vector<WeightWSM> weights{0, 0};
  weights.resize(100);
  for (unsigned nn = 2; nn < weights.size(); ++nn) {
    std::size_t a10 = rng.get_size_t(10);
    std::size_t a5 = rng.get_size_t(5);
    weights[nn] = 1 + a5 + 10 * (a10 / 8);
  }
  WeightedBinaryTree tree(weights, 4);

  std::vector<std::pair<VertexWSM, VertexWSM>> gates(200);
  const unsigned max_token = 30;
  for (unsigned ii = 0; ii < gates.size();) {
    gates[ii].first = rng.get_size_t(max_token);
    gates[ii].second = rng.get_size_t(max_token);
    if (gates[ii].first == gates[ii].second) {
      continue;
    }
    ++ii;
  }
  PatternGraphDataInput pgd_input;
  pgd_input.initial_gate_weight = 100;
  pgd_input.final_gate_weight = 20;

  // We'll try time-slicing and see what difference it makes.
  pgd_input.method =
      PatternGraphDataInput::ReorderingMethod::TIME_SLICES_OF_PARALLEL_GATES;
  const PatternGraphData p_graph_data(gates, pgd_input);
  CHECK(gates.size() == 200);
  CHECK(get_number_of_vertices(p_graph_data.pattern_graph_weights) == 31);

  // First, test the gates in the original order, i.e. the "wrong" order!
  // I.e., we've deliberately constructed a WSM problem for gates
  // occurring in time-sliced order, but then run through with the original
  // gates - an easy mistake to make!

  // To repeat, "gates" is DELIBERATELY the wrong order, since
  // p_graph_data was constructed by REORDERING the input gates
  // into time slices of parallel gates.
  run_end_to_end_IQP_and_check_solution(
      gates, p_graph_data.pattern_graph_weights, tree, 13193);

  // Now, some results found with try_random_placements.
  std::vector<CostedIQPSolution> costed_solutions_original_gates(3);

  // clang-format off
  // Quite a bit worse than WSM, took ~0.4 seconds.
  costed_solutions_original_gates[0] = {17119, {{0,2}, {1,23}, {2,19}, {3,86},
    {4,84}, {5,20}, {6,34}, {7,35}, {8,55}, {9,25}, {10,28}, {11,93}, {12,8},
    {13,88}, {14,5}, {15,41}, {16,31}, {17,6}, {18,75}, {19,22}, {20,82},
    {21,33}, {22,16}, {23,4}, {24,39}, {25,32}, {26,71}, {27,12}, {28,78},
    {29,53}, {30,14}}};

  // Took ~1.6 seconds.
  costed_solutions_original_gates[1] = {13244, {{0,2}, {1,42}, {2,93}, {3,86},
    {4,84}, {5,20}, {6,71}, {7,35}, {8,11}, {9,4}, {10,81}, {11,19}, {12,8},
    {13,52}, {14,5}, {15,45}, {16,87}, {17,6}, {18,77}, {19,32}, {20,82},
    {21,33}, {22,16}, {23,85}, {24,43}, {25,22}, {26,34}, {27,12}, {28,37},
    {29,53}, {30,14}}};

  // Good, but took >40 seconds.
  costed_solutions_original_gates[2] = {8202, {{0,9}, {1,21}, {2,17}, {3,5},
    {4,33}, {5,10}, {6,68}, {7,69}, {8,40}, {9,43}, {10,1}, {11,8}, {12,71},
    {13,3}, {14,22}, {15,7}, {16,86}, {17,6}, {18,4}, {19,16}, {20,41},
    {21,20}, {22,42}, {23,46}, {24,35}, {25,23}, {26,83}, {27,11}, {28,2},
    {29,70}, {30,34}}};
  // clang-format on
  test_known_solutions(costed_solutions_original_gates, gates, tree);

  // Now, use the REORDERED gates, i.e. the "correct" order.
  REQUIRE(p_graph_data.reordered_gates.size() == gates.size());
  std::vector<std::pair<VertexWSM, VertexWSM>> time_sliced_reordered_gates;
  time_sliced_reordered_gates.reserve(gates.size());
  for (const auto& entry : p_graph_data.reordered_gates) {
    time_sliced_reordered_gates.emplace_back(entry.gate);
  }

  run_end_to_end_IQP_and_check_solution(
      time_sliced_reordered_gates, p_graph_data.pattern_graph_weights, tree,
      13098);

  std::vector<CostedIQPSolution> costed_solutions_reordered_gates(3);

  // clang-format off
  // Took ~300 ms.
  costed_solutions_reordered_gates[0] = {17578, {{0,2}, {1,23}, {2,78}, {3,57},
    {4,45}, {5,20}, {6,34}, {7,30}, {8,55}, {9,35}, {10,28}, {11,48}, {12,7},
    {13,88}, {14,5}, {15,10}, {16,13}, {17,6}, {18,75}, {19,24}, {20,25},
    {21,33}, {22,16}, {23,59}, {24,39}, {25,32}, {26,71}, {27,12}, {28,19},
    {29,53}, {30,14}}};

  // Took ~4 seconds.
  costed_solutions_reordered_gates[1] = {13020, {{0,2}, {1,23}, {2,9}, {3,1},
    {4,83}, {5,21}, {6,17}, {7,15}, {8,96}, {9,48}, {10,7}, {11,3}, {12,86},
    {13,52}, {14,8}, {15,84}, {16,43}, {17,6}, {18,18}, {19,34}, {20,25},
    {21,5}, {22,70}, {23,59}, {24,29}, {25,4}, {26,82}, {27,12}, {28,19},
    {29,10}, {30,16}}};

  // After ~25 seconds.
  costed_solutions_reordered_gates[2] = {10868, {{0,97}, {1,11}, {2,14}, {3,3},
    {4,20}, {5,21}, {6,35}, {7,4}, {8,7}, {9,13}, {10,9}, {11,48}, {12,86},
    {13,24}, {14,8}, {15,42}, {16,2}, {17,41}, {18,18}, {19,17}, {20,25},
    {21,31}, {22,34}, {23,59}, {24,10}, {25,5}, {26,22}, {27,1}, {28,12},
    {29,43}, {30,16}}};
  // clang-format on
  test_known_solutions(
      costed_solutions_reordered_gates, time_sliced_reordered_gates, tree);
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
