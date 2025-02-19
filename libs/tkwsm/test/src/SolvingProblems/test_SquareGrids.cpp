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
  unsigned total_edges = 0;
  RNG rng;

  const std::vector<std::pair<unsigned, unsigned>> grid_sizes{
      {1, 1}, {1, 1}, {2, 1}, {1, 2},   {1, 3},   {3, 1},   {4, 1},
      {4, 2}, {5, 3}, {5, 5}, {10, 10}, {15, 15}, {20, 20},
  };

  for (int test_counter = 0; test_counter < 1; ++test_counter) {
    for (const auto& grid_size : grid_sizes) {
      grids.emplace_back();
      grids.back().width = grid_size.first;
      grids.back().height = grid_size.second;
      grids.back().fill_weights(rng);
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

static const std::vector<SquareGrid>& get_test_grids_ref() {
  static const auto square_grids = get_test_grids();
  return square_grids;
}

static std::vector<GraphEdgeWeights> get_graph_data() {
  const auto& grids = get_test_grids_ref();
  std::vector<GraphEdgeWeights> graph_data;
  graph_data.reserve(grids.size());
  for (const auto& grid : grids) {
    graph_data.emplace_back(grid.get_graph_edge_weights());
  }
  return graph_data;
}

static const std::vector<GraphEdgeWeights>& get_graph_data_ref() {
  static const auto graph_data = get_graph_data();
  return graph_data;
}

static std::map<std::pair<unsigned, unsigned>, WeightWSM>
get_known_solutions_map() {
  std::map<std::pair<unsigned, unsigned>, WeightWSM> map;
  const auto& grids = get_test_grids_ref();
  for (unsigned ii = 0; ii < grids.size(); ++ii) {
    for (unsigned jj = 0; jj < grids.size(); ++jj) {
      map[std::make_pair(ii, jj)] =
          grids[ii].get_subgraph_isomorphism_min_scalar_product(grids[jj]);
    }
  }
  return map;
}

static const std::map<std::pair<unsigned, unsigned>, WeightWSM>&
get_known_solutions_map_ref() {
  static const auto map = get_known_solutions_map();
  return map;
}

// g[i] -> g[j] problems which routinely take >200 ms,
// but still <1 second.
static std::set<std::pair<unsigned, unsigned>> get_slower_problems() {
  return {{6, 12}, {7, 12}, {8, 11}, {8, 12}, {9, 11}, {9, 12}, {12, 12}};
}

// Each g[i] -> g[j] problem taking >1 second.
static std::set<std::pair<unsigned, unsigned>> get_monster_problems() {
  return {{10, 11}, {10, 12}, {11, 12}};
}

// We're trying to split the tests up into shorter and longer tests.
// However, just testing shorter problems is not enough to make a short
// test, because there may be more shorter problems.
static void test(
    const std::set<std::pair<unsigned, unsigned>>& problems_to_skip,
    const std::set<std::pair<unsigned, unsigned>>& problems_to_run_if_nonempty,
    unsigned initial_problems_to_skip, unsigned max_problems_count,
    unsigned expected_problems_count, unsigned expected_total_scalar_product,
    unsigned timeout) {
  const auto& os = TestSettings::get().os;
  os << "\n\n### Square grids: expecting to test " << expected_problems_count
     << " problems;";
  if (!problems_to_skip.empty()) {
    os << " skipping " << problems_to_skip.size() << " problems;";
    REQUIRE(problems_to_run_if_nonempty.empty());
  }
  if (!problems_to_run_if_nonempty.empty()) {
    os << " running " << problems_to_run_if_nonempty.size() << " problems;";
    REQUIRE(problems_to_skip.empty());
    REQUIRE(expected_problems_count == problems_to_run_if_nonempty.size());
    REQUIRE(max_problems_count >= expected_problems_count);
  }

  const auto& grids = get_test_grids_ref();
  const auto& graph_data = get_graph_data_ref();

  unsigned skipped_problems_count = 0;
  unsigned tested_problems_count = 0;
  unsigned problems_count = 0;
  unsigned success_count = 0;
  unsigned timeout_count = 0;
  unsigned failure_count = 0;
  long long total_time_init = 0;
  long long total_time_search = 0;

  CheckedSolution::Statistics stats("square grids");
  bool break_out = false;

  for (unsigned ii = 0; ii < grids.size(); ++ii) {
    if (break_out) {
      break;
    }
    for (unsigned jj = 0; jj < grids.size(); ++jj) {
      ++problems_count;
      if (problems_count <= initial_problems_to_skip) {
        continue;
      }
      if (tested_problems_count >= max_problems_count) {
        break_out = true;
        break;
      }
      const auto pair = std::make_pair(ii, jj);
      if (problems_to_skip.count(pair) != 0) {
        ++skipped_problems_count;
        continue;
      }
      if (!problems_to_run_if_nonempty.empty() &&
          problems_to_run_if_nonempty.count(pair) == 0) {
        ++skipped_problems_count;
        continue;
      }
      ++tested_problems_count;
      const auto optimal_solution = get_known_solutions_map_ref().at(pair);
      CheckedSolution::ProblemInformation info;
      MainSolverParameters solver_params;

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
        CheckedSolution(
            graph_data[ii], graph_data[jj], info, solver_params, stats);
        continue;
      }

      // There is a known optimal solution.
      solver_params.timeout_ms = timeout;
      info.known_optimal_solution = optimal_solution;

      os << "\n#### g" << ii << " (" << grids[ii].width << "x"
         << grids[ii].height << ") -> g" << jj << " (" << grids[jj].width << "x"
         << grids[jj].height << ")";

      CheckedSolution(
          graph_data[ii], graph_data[jj], info, solver_params, stats);
    }
  }
  stats.finish();

  if (skipped_problems_count > 0) {
    os << "Skipped " << skipped_problems_count << " problems.\n";
  }
  // CHECK(skipped_problems_count + tested_problems_count == problems_count);
  CHECK(expected_problems_count == tested_problems_count);
  CHECK(max_problems_count >= tested_problems_count);
  CHECK(stats.timeout_count == 0);
  CHECK(stats.failure_count == 0);
}

// These are still individually very short: ~100 ms
std::set<std::pair<unsigned, unsigned>> get_tiny_pattern_longer_problems() {
  return {{0, 10}, {0, 11}, {0, 12}, {1, 10}, {1, 11}, {1, 12}};
}

// Embedding 1x1, 2x1, etc. is fast, even into large 20x20 targets.
SCENARIO("Trivial problems: tiny pattern square grids") {
  const std::set<std::pair<unsigned, unsigned>> empty_problems;
  test(
      get_tiny_pattern_longer_problems(), empty_problems, 0, 20, 20, 1000,
      // test coverage takes longer than normal running,
      // as does Valgrind;
      // 10ms timeout would be sufficient in normal runs.
      100 * 10);
}

SCENARIO("Easy problems: tiny pattern square grids") {
  const std::set<std::pair<unsigned, unsigned>> empty_problems;
  test(
      empty_problems, get_tiny_pattern_longer_problems(), 0, 20, 6, 1000,
      // Longer timeout for Valgrind.
      10 * 500);
}

SCENARIO(
    "short-to-medium problems: reasonable size pattern and target square "
    "grids") {
  const std::set<std::pair<unsigned, unsigned>> empty_problems;
  auto problems_to_skip = get_tiny_pattern_longer_problems();
  for (auto problem : get_slower_problems()) {
    CHECK(problems_to_skip.insert(problem).second);
  }
  for (auto problem : get_monster_problems()) {
    CHECK(problems_to_skip.insert(problem).second);
  }
  // Skip the first few problems, which are too easy
  test(problems_to_skip, empty_problems, 20, 1000, 136, 5000, 5000);
}

SCENARIO("Medium square grid problems only") {
  const std::set<std::pair<unsigned, unsigned>> empty_problems;
  test(empty_problems, get_slower_problems(), 0, 7, 7, 9999, 20000);
}

SCENARIO("Monster square grid problems only") {
  const std::set<std::pair<unsigned, unsigned>> empty_problems;
  test(empty_problems, get_monster_problems(), 0, 3, 3, 9999, 100000);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
