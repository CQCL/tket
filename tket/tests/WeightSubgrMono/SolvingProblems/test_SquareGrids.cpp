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

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/SquareGridGeneration.hpp"
#include "../TestUtils/TestSettings.hpp"

// It's easy to prove that embedding a grid (a rectangle in the 2D integer
// lattice Z^2, sides parallel to the axes) into another grid
// can only be done in the obvious ways (reflections; rotations by
// 90, 180, 270 degrees; translations),
// [provided that neither grid degenerates into a line],
// so we can easily calculate optimal solutions by brute force.

namespace tket {
namespace WeightedSubgraphMonomorphism {

static std::vector<SquareGrid> get_test_grids() {
  std::vector<SquareGrid> grids;
  size_t total_edges = 0;
  auto r_engine = GraphGeneration::get_r_engine();

  const std::vector<std::pair<unsigned, unsigned>> grid_sizes{
      {1, 1}, {1, 1}, {2, 1}, {1, 2},   {1, 3},   {3, 1},   {4, 1},
      {4, 2}, {5, 3}, {5, 5}, {10, 10}, {15, 15}, {20, 20},
  };

  for (int test_counter = 0; test_counter < 1; ++test_counter) {
    for (const auto& grid_size : grid_sizes) {
      grids.emplace_back();
      grids.back().width = grid_size.first;
      grids.back().height = grid_size.second;
      grids.back().fill_weights(r_engine);
      total_edges += grids.back().horiz_weights.size();
      total_edges += grids.back().vert_weights.size();
    }
  }
  const auto& os = TestSettings::get().os;
  os << "\n\n###########\n\n######### generated " << grids.size()
     << " square grids, total " << total_edges << " edges: [";
  for (const auto& size : grid_sizes) {
    os << size.first << "x" << size.second << " ";
  }
  os << "]";
  return grids;
}

SCENARIO("test searching with square grids") {
  const unsigned timeout = 20000;
  const auto grids = get_test_grids();
  const auto& os = TestSettings::get().os;
  os << "\n### Timeout " << timeout << "\n";

  // Let's try aggressively pruning the total weight, see if it speeds things
  // up. i.e., we restrict the max weight to be the given percentage of the
  // actual known optimal solution. To prevent overflow, don't actually
  // calculate with percentages this big.
  const unsigned max_actual_percentage = 10000;

  // All the (i,j) pairs for which G(i) --> G(j) routinely takes >0.5 seconds.
  // (But still < 8 seconds).
  const std::set<std::pair<unsigned, unsigned>> harder_problems{
      {7, 12}, {8, 11},  {8, 12},  {9, 11},
      {9, 12}, {10, 11}, {10, 12}, {11, 12}};
  // const bool skip_harder_problems = true;
  const bool skip_harder_problems = false;
  unsigned skipped_problems_count = 0;
  /*
  // Try different weight pruning factors.
  const std::vector<unsigned> scalar_product_percentages {
      max_actual_percentage, 200, 150, 100, 50 };
  //*/
  const std::vector<unsigned> scalar_product_percentages{max_actual_percentage};

  for (unsigned scalar_product_percentage : scalar_product_percentages) {
    size_t success_count = 0;
    size_t timeout_count = 0;
    size_t failure_count = 0;
    long long total_time_init = 0;
    long long total_time_search = 0;
    std::vector<GraphEdgeWeights> gdata;
    gdata.reserve(grids.size());
    for (const auto& grid : grids) {
      gdata.emplace_back(grid.get_graph_edge_weights());
    }
    os << "\n#### now testing all against all:";
    if (scalar_product_percentage < max_actual_percentage) {
      os << " AGGRESSIVE weight squeezing: " << scalar_product_percentage
         << "%";
    }

    CheckedSolution::Statistics stats;

    for (unsigned ii = 0; ii < grids.size(); ++ii) {
      for (unsigned jj = 0; jj < grids.size(); ++jj) {
        if (skip_harder_problems && harder_problems.count({ii, jj}) != 0) {
          ++skipped_problems_count;
          continue;
        }

        const WeightWSM optimal_solution =
            grids[ii].get_subgraph_isomorphism_min_scalar_product(grids[jj]);
        CheckedSolution::ProblemInformation info;
        MainSolver::Parameters solver_params;

        if (optimal_solution == 0) {
          info.existence = CheckedSolution::ProblemInformation::
              SolutionsExistence::KNOWN_TO_BE_INSOLUBLE;
          // Where no square grid embedding exists, it's trivial to prove;
          // just counting vertices is enough.
          // (All large grids are square; rectangles would not be so easy!
          // E.g. embedding 5x1 into 4x4 is impossible, but you need widths
          // and heights to see that easily; counting vertices and edges is
          // insufficient. But, only small grids are non-square here).
          solver_params.timeout_ms = 100;
          CheckedSolution(gdata[ii], gdata[jj], info, solver_params, stats);
          continue;
        }

        // There is a known optimal solution.
        solver_params.timeout_ms = timeout;
        info.known_optimal_solution = optimal_solution;

        os << "\n#### g" << ii << " (" << grids[ii].width << "x"
           << grids[ii].height << ") -> g" << jj << " (" << grids[jj].width
           << "x" << grids[jj].height << ")";

        if (scalar_product_percentage < max_actual_percentage) {
          const auto weight_constraint =
              (scalar_product_percentage * optimal_solution) / 100;
          solver_params.weight_upper_bound_constraint = weight_constraint;
          os << " : SQUEEZE " << weight_constraint;
        }
        CheckedSolution(gdata[ii], gdata[jj], info, solver_params, stats);
      }
    }
    os << "\n\n### FINAL time (ms): " << stats.total_init_time_ms << "+"
       << stats.total_search_time_ms << "; " << stats.success_count
       << " success; " << stats.failure_count << " failures; "
       << stats.timeout_count << " timeouts.";

    if (skipped_problems_count > 0) {
      os << " Skipped " << skipped_problems_count << " problems.";
    }
    if (skip_harder_problems) {
      CHECK(harder_problems.size() == skipped_problems_count);
    } else {
      CHECK(skipped_problems_count == 0);
    }
    CHECK(stats.success_count == 169 - skipped_problems_count);
    CHECK(stats.timeout_count == 0);
    CHECK(stats.failure_count == 0);
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
