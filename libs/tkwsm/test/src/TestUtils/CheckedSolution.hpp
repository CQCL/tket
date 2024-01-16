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

#pragma once
#include <optional>
#include <tkwsm/EndToEndWrappers/MainSolver.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

// Full end-to-end solve and automatic checking of the solution.
// Only for use with single solutions, i.e.
// for_multiple_full_solutions_the_max_number_to_obtain should be set to 0.
struct CheckedSolution {
  // Extra information about the problem to be solved.
  struct ProblemInformation {
    std::optional<WeightWSM> known_optimal_solution;
    std::optional<WeightWSM> known_lower_bound;
    std::optional<WeightWSM> known_upper_bound;

    enum class SolutionsExistence {
      KNOWN_TO_BE_SOLUBLE,
      KNOWN_TO_BE_INSOLUBLE,
      UNKNOWN
    };
    SolutionsExistence existence = SolutionsExistence::UNKNOWN;
  };

  struct Statistics {
    explicit Statistics(const std::string& test_name);
    Statistics(const std::string& test_name, std::size_t number_of_graphs);

    enum class Expectation {
      ALL_SUCCESS,
      ALL_SUCCESS_OR_TIMEOUT,
      SUCCESS_FAILURE_TIMEOUTS_ALL_ALLOWED
    };
    void finish(Expectation expectation = Expectation::ALL_SUCCESS) const;

    unsigned success_count;
    unsigned failure_count;
    unsigned timeout_count;
    long long total_init_time_ms;
    long long total_search_time_ms;
    unsigned long long total_iterations;
  };

  // Did the solver time out?
  bool finished = false;

  std::size_t iterations = 0;

  // The best solution found. If nonempty, it is complete,
  // and has been checked to be valid.
  std::vector<std::pair<VertexWSM, VertexWSM>> assignments;

  // If assignments is nonempty, the reported scalar product of the solution
  // (automatically checked to be correct).
  WeightWSM scalar_product = 0;

  // This is rare, but occasionally an impossible TV is detected.
  std::vector<VertexWSM> impossible_target_vertices;

  // Solve the given problem, updating the statistics.
  CheckedSolution(
      const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
      ProblemInformation info, const MainSolverParameters& solver_params,
      Statistics& stats);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
