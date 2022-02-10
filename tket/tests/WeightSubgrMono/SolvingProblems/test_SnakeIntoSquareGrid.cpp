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
#include <map>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/GraphGeneration.hpp"
#include "../TestUtils/ProblemGeneration.hpp"
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

// KEY: is the problem name  VALUE: the collection of solved problems.
static std::map<std::string, ProblemGeneration::EncodedSquareGrid> get_data() {
  return {
      {"Uniform1, small weights",
       {0x1093fb7292ecde4, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14,
        16, 18}},

      {"Uniform2, small weights",
       {0x9372a0ee562901cc, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15,
        17, 18}},
      {"Uniform3, small weights",
       {0x196df104e143cde2, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15,
        17, 18}},
      {"Uniform4, small weights",
       {0x4e1bc8532fd80f73, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15,
        17, 18}},
      {"Uniform5, small weights",
       {0xadf9bf4ee6c8c7a0, 2, 3, 4, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 17,
        18, 21}},
      {"Uniform6, small weights",
       {0x9372a0ee562901cc, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15,
        17, 18}},

      {"Uniform1, large weights",
       {0x1093fb7292ecde4,
        10,
        100,
        1000,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        20,
        21,
        22,
        32,
        123,
        133,
        233,
        234,
        343,
        344}},
      {"Uniform2, large weights",
       {0x9372a0ee562901cc,
        4,
        9,
        30,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        12,
        13,
        14,
        15,
        19,
        23,
        24,
        38}},
      {"Uniform3, large weights",
       {0x196df104e143cde2,
        4,
        9,
        30,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        13,
        14,
        18,
        19,
        23,
        26,
        28,
        34,
        35,
        50,
        51}},
      {"Uniform4, large weights",
       {0x4e1bc8532fd80f73,
        4,
        9,
        30,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        11,
        12,
        13,
        17,
        18,
        19,
        23,
        24,
        28,
        32,
        33,
        42,
        43}},
      {"Uniform5, large weights",
       {0xadf9bf4ee6c8c7a0,
        4,
        9,
        30,
        1,
        2,
        3,
        4,
        5,
        9,
        10,
        11,
        12,
        16,
        17,
        18,
        27,
        28,
        37,
        38,
        50,
        54,
        63,
        72}},
      {"Uniform6, large weights",
       {
           0x9372a0ee562901cc,
           4,
           9,
           30,
           1,
           2,
           3,
           4,
           5,
           6,
           7,
           8,
           12,
           13,
           14,
           15,
           19,
           23,
           24,
           38,
       }},

  };
}

SCENARIO("embedding paths into square grids") {
  const auto solved_problems_map = get_data();
  ProblemGeneration::EncodedSquareGrid calc_problems;
  unsigned timeout = 3000;
  CheckedSolution::ProblemInformation info;
  info.existence = CheckedSolution::ProblemInformation::SolutionsExistence::
      KNOWN_TO_BE_SOLUBLE;

  CheckedSolution::Statistics statistics;
  const MainSolver::Parameters solver_params(timeout);
  const auto& os = TestSettings::get().os;

  for (const auto& entry : solved_problems_map) {
    os << "\nTesting '" << entry.first << "', square grid target " << std::hex
       << entry.second[0] << std::dec
       << ", embedding shakes of length <= " << entry.second.size() - 3
       << ", timeout=" << timeout << ":";
    const auto t_graph_data =
        ProblemGeneration::get_target_graph_for_encoded_square_grid(
            entry.second);
    calc_problems.clear();

    // Get the line graphs.
    for (unsigned ii = 0; ii < entry.second.size(); ++ii) {
      if (ii < 4) {
        // The first few entries encode the grid weights,
        // NOT the expected final scalar product.
        calc_problems.push_back(entry.second[ii]);
        continue;
      }
      const auto line_graph = GraphGeneration::get_line(ii - 2, false);
      for (const auto& line_entry : line_graph) {
        REQUIRE(line_entry.second == 1);
      }

      const CheckedSolution checked_solution(
          line_graph, t_graph_data, info, solver_params, statistics);
      // Should be no timeouts, and a complete solution.
      CHECK(checked_solution.finished);
      CHECK(checked_solution.complete_solution_weight);
      if (checked_solution.complete_solution_weight) {
        calc_problems.push_back(
            checked_solution.complete_solution_weight.value());
      }
    }
    CHECK(entry.second == calc_problems);
  }
  os << "\nFIN snakes into grids: total time " << statistics.total_init_time_ms
     << "+" << statistics.total_search_time_ms << " ms.";
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
