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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <set>
#include <sstream>
#include <string>
#include <tkwsm/EndToEndWrappers/MainSolver.hpp>

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

  const unsigned number_of_tests = 10;

  const std::set<std::string> all_embeddings{
      "012345", "012354", "312045", "312054"};
  std::set<std::string> calc_solution_strings;

  // Record the number of iterations also.
  std::stringstream iterations_ss;
  iterations_ss << "[ ";

  for (parameters.for_multiple_full_solutions_the_max_number_to_obtain = 0;
       parameters.for_multiple_full_solutions_the_max_number_to_obtain < 10;
       ++parameters.for_multiple_full_solutions_the_max_number_to_obtain) {
    const MainSolver solver(pattern_graph, pattern_graph, parameters);
    const auto& solution_data = solver.get_solution_data();
    iterations_ss << solution_data.iterations << " ";

    CHECK(solution_data.trivial_weight_lower_bound == edges.size());
    CHECK(solution_data.trivial_weight_initial_upper_bound == edges.size());

    const auto& solutions = solution_data.solutions;
    if (parameters.for_multiple_full_solutions_the_max_number_to_obtain == 0) {
      CHECK(solutions.size() == 1);
      CHECK(solution_data.finished);
    } else {
      if (parameters.for_multiple_full_solutions_the_max_number_to_obtain <=
          4) {
        CHECK(
            solutions.size() ==
            parameters.for_multiple_full_solutions_the_max_number_to_obtain);
        CHECK(!solution_data.finished);
      } else {
        CHECK(solutions.size() == 4);
        CHECK(solution_data.finished);
      }
    }

    calc_solution_strings.clear();

    for (const SolutionWSM& solution : solutions) {
      CHECK(solution.scalar_product == edges.size());
      CHECK(solution.total_p_edges_weight == edges.size());
      std::stringstream ss;
      for (unsigned ii = 0; ii < solution.assignments.size(); ++ii) {
        CHECK(solution.assignments[ii].first == ii);
        ss << solution.assignments[ii].second;
      }
      const auto solution_str = ss.str();
      calc_solution_strings.insert(solution_str);
      CHECK(all_embeddings.count(solution_str) != 0);
    }
    // Every solution should be different.
    CHECK(calc_solution_strings.size() == solutions.size());
  }
  iterations_ss << "]";
  CHECK(iterations_ss.str() == "[ 1 1 2 3 4 5 5 5 5 5 ]");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
