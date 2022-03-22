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
#include <catch2/catch.hpp>
#include <set>
#include <sstream>
#include <string>

#include "WeightSubgrMono/EndToEndWrappers/MainSolver.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("Single fixed graph: multiple self embeddings") {
  // Draw it! (0,3) and (4,5) can be swapped, but everything else is fixed.
  const std::vector<std::pair<unsigned, unsigned>> edges{
      {0, 1}, {1, 3}, {0, 3}, {1, 2}, {2, 4}, {2, 5}};

  GraphEdgeWeights pattern_graph;
  for (const auto& edge : edges) {
    pattern_graph[get_edge(edge.first, edge.second)] = 1;
  }
  MainSolverParameters parameters;

  // Actually, <<1ms should be enough
  parameters.timeout_ms = 10;

  // Element[i] is the number of iterations for max_number = i
  // (platform independent, but should always be increasing).
  const unsigned number_of_tests = 10;
  std::stringstream iterations_ss;
  iterations_ss << "[ ";
  const std::set<std::string> all_embeddings{
      "012345", "012354", "312045", "312054"};
  std::set<std::string> calc_solution_strings;

  for (parameters.for_multiple_full_solutions_the_max_number_to_obtain = 0;
       parameters.for_multiple_full_solutions_the_max_number_to_obtain < 10;
       ++parameters.for_multiple_full_solutions_the_max_number_to_obtain) {
    MainSolver solver;
    solver.solve(pattern_graph, pattern_graph, parameters);
    const auto& stored_solutions = solver.get_some_full_solutions();
    calc_solution_strings.clear();
    for (const auto& solution : stored_solutions) {
      std::stringstream ss;
      for (unsigned ii = 0; ii < solution.size(); ++ii) {
        CHECK(solution[ii].first == ii);
        ss << solution[ii].second;
      }
      calc_solution_strings.insert(ss.str());
      CHECK(all_embeddings.count(ss.str()) != 0);
    }
    CHECK(calc_solution_strings.size() == stored_solutions.size());
    if (parameters.for_multiple_full_solutions_the_max_number_to_obtain <= 4) {
      CHECK(
          calc_solution_strings.size() ==
          parameters.for_multiple_full_solutions_the_max_number_to_obtain);
    } else {
      CHECK(calc_solution_strings.size() == 4);
    }
    iterations_ss << solver.get_solution_statistics().iterations << " ";
  }
  iterations_ss << "]";
  CHECK(iterations_ss.str() == "[ 4 0 1 2 3 4 4 4 4 4 ]");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
