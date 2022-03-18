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

#pragma once
#include <optional>

#include "WeightSubgrMono/EndToEndWrappers/MainSolver.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

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
    unsigned success_count = 0;
    unsigned failure_count = 0;
    unsigned timeout_count = 0;
    long long total_init_time_ms = 0;
    long long total_search_time_ms = 0;
  };

  // We don't really care about incomplete solutions.
  // If non-null, it means that the solver found a complete solution
  // (however, it may have timed out).
  std::optional<WeightWSM> complete_solution_weight;
  bool finished = false;

  std::size_t iterations = 0;

  // The best solution found. This may or may not be complete.
  // If "complete_solution_weight" is not null, then it IS complete,
  // and has been checked to be valid.
  std::vector<std::pair<VertexWSM, VertexWSM>> assignments;

  // Solve the given problem, updating the statistics.
  CheckedSolution(
      const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
      ProblemInformation info, const MainSolverParameters& solver_params,
      Statistics& stats,
      const std::vector<std::pair<VertexWSM, VertexWSM>>& suggested_assignments = {});
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
