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

#include "ResumedSolutionChecker.hpp"

#include <catch2/catch_test_macros.hpp>

#include "CheckedSolution.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

void ResumedSolutionChecker::check(
    const CheckedSolution& solution, const GraphEdgeWeights& pdata,
    const GraphEdgeWeights& tdata, MainSolverParameters solver_params) {
  if (m_remaining_problems == 0 ||
      !(solution.iterations >= m_min_number_of_iterations &&
        solution.iterations <= m_max_number_of_iterations)) {
    return;
  }
  --m_remaining_problems;
  REQUIRE(solution.iterations <= solver_params.iterations_timeout);
  const std::size_t iterations_addition =
      solution.iterations / m_number_of_chunks;
  solver_params.iterations_timeout = iterations_addition;
  MainSolver solver(pdata, tdata, solver_params);
  const auto& paused_data = solver.get_solution_data();
  CHECK(paused_data.iterations == solver_params.iterations_timeout);
  CHECK(!paused_data.finished);

  for (;;) {
    solver_params.iterations_timeout += iterations_addition;
    bool should_break = false;
    // Once we get close to the number of iterations, fix on it exactly.
    if (solver_params.iterations_timeout + 2 > solution.iterations) {
      solver_params.iterations_timeout = solution.iterations;
      should_break = true;
    }
    solver.solve(solver_params);
    if (should_break) {
      break;
    }
  }
  const auto& final_data = solver.get_solution_data();
  CHECK(final_data.finished == solution.finished);
  CHECK(final_data.iterations == solution.iterations);
  CHECK(solution.assignments.empty() == final_data.solutions.empty());
  if (!solution.assignments.empty() && !final_data.solutions.empty()) {
    // This should be the best.
    const auto& new_best_solution = final_data.solutions[0];
    CHECK(new_best_solution.assignments == solution.assignments);
    CHECK(new_best_solution.get_errors(pdata, tdata) == "");
    CHECK(new_best_solution.scalar_product == solution.scalar_product);
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
