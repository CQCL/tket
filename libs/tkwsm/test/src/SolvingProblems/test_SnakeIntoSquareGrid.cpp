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
#include <map>
#include <sstream>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/GraphGeneration.hpp"
#include "../TestUtils/ProblemGeneration.hpp"
#include "../TestUtils/ResumedSolutionChecker.hpp"
#include "../TestUtils/TestSettings.hpp"

/*
Let's try embedding paths (lines) of length 2,3,4,5,... into
5x5 square grids, with each line edge having weight 1,
to give some fixed problems for testing/benchmarking.

Note that vertex local pruning/filtering is (almost)
completely useless for this,
because (almost) every p-vertex (in the line graph) can be mapped to
every t-vertex (in the square grid).

[Fun exercise: for the 9 points (x,y) with x,y in {0,1,2},
joined with horiz/vert grid edges,
no snake starting at (1,0) can cover every point.
What happens for general WxH grids?!]

Thus, the times depend heavily on WEIGHT-based pruning.
*/

namespace tket {
namespace WeightedSubgraphMonomorphism {

// KEY: is the problem name
// VALUE: the collection of solved problems.
//    In each problem, the very last entry of the vector gives
//    the solution value cutoff point for deciding between short/long
//    tests (roughly, as soon as a single problem takes ~50ms).
//    (Although, this means that a "long" test could be shorter than
//    a "short" test, because it has fewer problems).
static std::map<std::string, ProblemGeneration::EncodedSquareGrid> get_data() {
  return {
      // clang-format off
      {"Uniform1, small weights",
       {0x1093fb7292ecde4, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14,
        16, 18, 16}},
      {"Uniform2, small weights",
       {0x9372a0ee562901cc, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15,
        17, 18, 15}},
      {"Uniform3, small weights",
       {0x196df104e143cde2, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15,
        17, 18, 14}},
      {"Uniform4, small weights",
       {0x4e1bc8532fd80f73, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15,
        17, 18, 17}},
      {"Uniform5, small weights",
       {0xadf9bf4ee6c8c7a0, 2, 3, 4, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 17,
        18, 21, 17}},
      {"Uniform6, small weights",
       {0x9372a0ee562901cc, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15,
        17, 18, 15}},
      {"Uniform1, large weights",
       {0x1093fb7292ecde4, 10, 100, 1000, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 21,
        22, 32, 123, 133, 233, 234, 343, 344, 123}},
      {"Uniform2, large weights",
       {0x9372a0ee562901cc, 4, 9, 30, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15,
        19, 23, 24, 38, 23}},
      {"Uniform3, large weights",
       {0x196df104e143cde2, 4, 9, 30, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 18, 19,
        23, 26, 28, 34, 35, 50, 51, 19}},
      {"Uniform4, large weights",
       {0x4e1bc8532fd80f73, 4, 9, 30, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 17, 18,
        19, 23, 24, 28, 32, 33, 42, 43, 23}},
      {"Uniform5, large weights",
       {0xadf9bf4ee6c8c7a0, 4, 9, 30, 1, 2, 3, 4, 5, 9, 10, 11, 12, 16, 17, 18,
        27, 28, 37, 38, 50, 54, 63, 72, 27}},
      {"Uniform6, large weights",
       {0x9372a0ee562901cc, 4, 9, 30, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15,
        19, 23, 24, 38, 23}},
      // clang-format on
  };
}

static void test(
    bool short_test, unsigned expected_total_solution_value,
    unsigned number_of_targets_to_test, unsigned number_of_targets_to_skip) {
  const auto solved_problems_map = get_data();
  ProblemGeneration::EncodedSquareGrid calc_problems;
  unsigned timeout = short_test ? 1000 : 10000;
  CheckedSolution::ProblemInformation info;
  info.existence = CheckedSolution::ProblemInformation::SolutionsExistence::
      KNOWN_TO_BE_SOLUBLE;

  const MainSolverParameters solver_params(timeout);
  const auto& os = TestSettings::get().os;
  unsigned total_solution_value = 0;
  unsigned target_graph_count = 0;
  unsigned target_graph_test_count = 0;
  unsigned problem_count = 0;

  const std::string short_test_str = short_test ? "SHORT" : "LONG";
  std::stringstream ss;
  ss << "Embedding snakes into square grids; skipping "
     << number_of_targets_to_skip << " initial targets; testing max "
     << number_of_targets_to_test << " targets; " << short_test_str
     << " test problems";

  CheckedSolution::Statistics statistics(ss.str());
  ResumedSolutionChecker resumption_checker;

  for (const auto& entry : solved_problems_map) {
    ++target_graph_count;
    if (target_graph_count <= number_of_targets_to_skip) {
      continue;
    }
    if (target_graph_test_count >= number_of_targets_to_test) {
      break;
    }
    ++target_graph_test_count;
    // What are the acceptable solution values? We only consider problems
    // where the solution value lies within this interval.
    const unsigned min_solution_value = short_test ? 0 : entry.second.back();
    const unsigned max_solution_value =
        short_test ? entry.second.back() - 1
                   : entry.second.at(entry.second.size() - 2);

    os << "\nEmbedding snakes: '" << entry.first << "', square grid target "
       << std::hex << entry.second[0] << std::dec << "; " << short_test_str
       << " test; timeout=" << timeout << "; only values in ["
       << min_solution_value << "," << max_solution_value << "]";

    // The solution values have a monotonic property: obvious,
    // since the line pattern graphs all have edge weight 1.
    // Also, the cutoff point should be an actual value,
    // and not the last.
    {
      unsigned cutoff_point_count = 0;
      for (unsigned jj = 5; jj + 1 < entry.second.size(); ++jj) {
        REQUIRE(entry.second[jj - 1] < entry.second[jj]);
        if (entry.second[jj] == entry.second.back()) {
          ++cutoff_point_count;
        }
      }
      REQUIRE(cutoff_point_count == 1);
      REQUIRE(entry.second.back() < entry.second.at(entry.second.size() - 2));
    }

    const auto t_graph_data =
        ProblemGeneration::get_target_graph_for_encoded_square_grid(
            entry.second);
    calc_problems.clear();

    // Get the line graphs.
    for (unsigned ii = 0; ii + 1 < entry.second.size(); ++ii) {
      if (ii < 4 || (!(entry.second[ii] >= min_solution_value &&
                       entry.second[ii] <= max_solution_value))) {
        // The first few entries encode the grid weights,
        // NOT the expected final scalar product.
        // Otherwise, if it's outside the range, don't test it.
        calc_problems.push_back(entry.second[ii]);
        continue;
      }
      ++problem_count;
      total_solution_value += entry.second[ii];
      const auto line_graph = GraphGeneration::get_line(ii - 2, false);
      for (const auto& line_entry : line_graph) {
        REQUIRE(line_entry.second == 1);
      }
      const CheckedSolution checked_solution(
          line_graph, t_graph_data, info, solver_params, statistics);

      resumption_checker.check(
          checked_solution, line_graph, t_graph_data, solver_params);

      // Should be no timeouts, and a complete solution.
      CHECK(checked_solution.finished);
      calc_problems.push_back(checked_solution.scalar_product);
    }
    // Remember the final entry, used to split up the tests by time
    calc_problems.push_back(entry.second.back());
    CHECK(entry.second == calc_problems);
  }
  statistics.finish();
  CHECK(expected_total_solution_value == total_solution_value);
}

SCENARIO("embedding paths into square grids - quicker problems, fewer tests") {
  test(true, 150, 1, 0);
}

SCENARIO("embedding paths into square grids - quicker problems, more tests") {
  // We test only quicker problems, but MORE of them;
  // so this is a "long" test.
  test(true, 1048, 100, 1);
}

SCENARIO("embedding paths into square grids - all slower problems") {
  // Test all slower problems, but MORE of them;
  // so this is a "long" test.
  test(false, 2729, 100, 0);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
