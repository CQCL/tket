// Copyright Quantinuum
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
#include <chrono>
#include <iostream>
#include <numeric>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/GraphTheoretic/NeighboursData.hpp>
#include <tkwsm/GraphTheoretic/VertexRelabelling.hpp>
#include <tkwsm/InitPlacement/MonteCarloCompleteTargetSolution.hpp>
#include <tkwsm/InitPlacement/UtilsIQP.hpp>

#include "TestWeightedGraphData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

typedef std::chrono::steady_clock Clock;

namespace {
struct MCCTSolution {
  std::vector<unsigned> tv_assignments;
  WeightWSM scalar_product;
};
}  // namespace

static void test_random_graph_data(
    const GraphEdgeWeights& pattern_graph_data,
    const GraphEdgeWeights& explicit_target_graph_data,
    // The last element should be the one actually produced by the run.
    const std::vector<MCCTSolution>& solutions,
    unsigned final_number_of_iterations, bool verbose = false) {
  REQUIRE(!solutions.empty());
  WeightWSM implicit_target_weight =
      2 * get_max_weight(explicit_target_graph_data);

  const NeighboursData pattern_ndata(pattern_graph_data);
  const NeighboursData target_ndata(explicit_target_graph_data);

  // Check that the given solutions are actually all valid.
  for (const MCCTSolution& solution : solutions) {
    CHECK(
        get_scalar_product_with_complete_target(
            pattern_ndata, target_ndata, implicit_target_weight,
            solution.tv_assignments) == solution.scalar_product);
  }

  const auto start_time = Clock::now();
  const MonteCarloCompleteTargetSolution calc_solution(
      pattern_ndata, target_ndata, implicit_target_weight, 1000000);
  const auto end_time = Clock::now();
  if (verbose) {
    const auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_time - start_time)
                             .count();

    std::cerr << "\nProblem: "
              << pattern_ndata.get_number_of_nonisolated_vertices() << " pv, "
              << target_ndata.get_number_of_nonisolated_vertices() << " tv, "
              << calc_solution.iterations() << " iterations; took " << time_ms
              << " ms.\n";
  }
  CHECK(
      calc_solution.get_best_scalar_product() ==
      get_scalar_product_with_complete_target(
          pattern_ndata, target_ndata, implicit_target_weight,
          calc_solution.get_best_assignments()));
  CHECK(calc_solution.iterations() == final_number_of_iterations);
  CHECK(
      calc_solution.get_best_assignments() == solutions.back().tv_assignments);
  CHECK(
      calc_solution.get_best_scalar_product() ==
      solutions.back().scalar_product);
}

SCENARIO(
    "Monte Carlo solutions for random complete target graph: small and "
    "medium") {
  std::vector<MCCTSolution> solutions(4);

  // The last element in "solutions" is always the actual calculated solution.
  // Obviously these will change as the MCCT algorithm and default parameters
  // are changed.
  {
    // Small test.
    RNG rng;
    const auto pattern_graph_data = get_graph_data(rng, 6, 8, 1000, 2000);
    const auto explicit_target_graph_data =
        get_graph_data(rng, 10, 15, 10, 100);
    // clang-format off
    solutions[0].tv_assignments = { 1, 5, 7, 0, 4, 2};
    solutions[0].scalar_product = 1090000;
    solutions[1].tv_assignments = {5, 2, 4, 0, 7, 1};
    solutions[1].scalar_product = 675000;
    solutions[2].tv_assignments = { 7, 1, 5, 4, 9, 2 };
    solutions[2].scalar_product = 320000;
    //solutions[3].tv_assignments = { 9, 2, 5, 8, 1, 0 };
    solutions[3].tv_assignments = { 8, 5, 1, 3, 9, 2 };
    solutions[3].scalar_product = 290000;
    // clang-format on
    test_random_graph_data(
        pattern_graph_data, explicit_target_graph_data, solutions, 1000000);
  }
  {
    // Medium test.
    RNG rng;
    const auto pattern_graph_data = get_graph_data(rng, 10, 20, 1000, 2000);
    const auto explicit_target_graph_data =
        get_graph_data(rng, 20, 30, 10, 100);
    // clang-format off
    solutions[0].tv_assignments = { 19, 1, 5, 15, 13, 11, 7, 0, 4, 2};
    solutions[0].scalar_product = 3380384;
    solutions[1].tv_assignments = {19, 1, 15, 7, 12, 5, 11, 0, 4, 10};
    solutions[1].scalar_product = 1682984;
    solutions[2].tv_assignments = { 4, 19, 12, 10, 13, 18, 15, 5, 9, 0 };
    solutions[2].scalar_product = 1359440;
    solutions[2].tv_assignments = { 4, 19, 12, 10, 13, 18, 15, 5, 9, 0 };
    solutions[2].scalar_product = 1359440;
    //solutions[3].tv_assignments = { 19, 1, 5, 4, 0, 9, 14, 18, 16, 8 };
    //solutions[3].scalar_product = 1114518;
    solutions[3].tv_assignments = { 18, 15, 12, 13, 19, 0, 4, 5, 10, 1 };
    solutions[3].scalar_product = 1304714;
    //test_random_graph_data(pattern_graph_data, explicit_target_graph_data,
    //          solutions, 124594);
    //  clang-format on
    test_random_graph_data(
        pattern_graph_data, explicit_target_graph_data, solutions, 125673);
  }
}

SCENARIO("Monte Carlo solutions for random complete target graph: large") {
  RNG rng;
  std::vector<MCCTSolution> solutions(3);
  const auto pattern_graph_data = get_graph_data(rng, 50, 300, 1000, 2000);
  const auto explicit_target_graph_data = get_graph_data(rng, 60, 500, 10, 100);
  // clang-format off
  solutions[0].tv_assignments = {19, 25, 1, 5, 32, 49, 41, 29, 15, 26, 13, 24,
    20, 11, 53, 37, 47, 7, 0, 46, 4, 40, 58, 56, 2, 44, 35, 43, 52, 34, 17, 3,
    55, 28, 59, 23, 22, 10, 38, 16, 31, 8, 14, 50, 39, 6, 21, 33, 30, 18};
  solutions[0].scalar_product = 59155064;
  solutions[1].tv_assignments = {14, 26, 1, 36, 29, 18, 46, 57, 15, 35, 13,
    40, 22, 48, 49, 27, 19, 34, 17, 20, 53, 31, 56, 39, 50, 4, 7, 58, 32, 52,
    43, 6, 55, 30, 9, 51, 23, 5, 59, 25, 45, 44, 24, 0, 33, 3, 10, 16, 47, 21};
  solutions[1].scalar_product = 42161896;
  solutions[2].tv_assignments = { 18, 8, 39, 15, 3, 45, 16, 43, 14, 20, 4, 51,
    36, 44, 6, 0, 29, 25, 53, 42, 54, 22, 30, 58, 47, 1, 13, 57, 55, 52, 31,
    56, 21, 41, 5, 9, 33, 46, 17, 19, 7, 48, 50, 40, 23, 27, 49, 11, 26, 12};
  solutions[2].scalar_product = 36316240;
  // clang-format on
  test_random_graph_data(
      pattern_graph_data, explicit_target_graph_data, solutions, 1000000);
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
