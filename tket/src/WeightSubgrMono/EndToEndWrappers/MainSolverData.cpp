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

#include "WeightSubgrMono/EndToEndWrappers/MainSolverData.hpp"

#include <algorithm>
#include <chrono>

#include "Utils/Assert.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

typedef std::chrono::steady_clock Clock;

ReductionResult MainSolverData::initialise(
    const GraphEdgeWeights& pattern_edges,
    const GraphEdgeWeights& target_edges) {
  statistics = SolutionStatistics();
  initialised = true;
  previous_upper_bound_constraint = {};
  shared_data_ptr.reset();

  if (!fixed_data.initialise(pattern_edges, target_edges)) {
    return ReductionResult::FAILURE;
  }
  shared_data_ptr = std::make_unique<SharedData>(fixed_data);
  if (!shared_data_ptr) {
    return ReductionResult::FAILURE;
  }
  return shared_data_ptr->initialise(branch);
}


void MainSolverData::do_one_solve_iteration_with_suggestion(
      const std::vector<std::pair<VertexWSM, VertexWSM>>& suggested_assignments) {
  TKET_ASSERT(initialised);
  TKET_ASSERT(shared_data_ptr);
  const auto reduction_result = shared_data_ptr->search_with_suggestion(
        branch, var_ordering, val_ordering, suggested_assignments);
  if (reduction_result == ReductionResult::FINISHED) {
    statistics.finished = true;
  }
}


void MainSolverData::solve_loop_after_initialisation(const MainSolverParameters& parameters) {
  TKET_ASSERT(initialised);
  TKET_ASSERT(shared_data_ptr);
  const auto solve_start_time = Clock::now();
  const auto solve_timeout_time =
      solve_start_time + std::chrono::milliseconds(parameters.timeout_ms);

  const auto& solution =
      shared_data_ptr->solution_storage.best_solution();

  // We don't want to compute time twice
  // (even though it's hopefully cheap; but it's not free!)
  bool time_added = false;

  for (; statistics.iterations < parameters.max_iterations;
       statistics.iterations += 1) {

    const auto reduction_result = shared_data_ptr->search(
        branch, var_ordering, val_ordering);

    if (reduction_result == ReductionResult::FINISHED) {
      statistics.finished = true;
      break;
    }
    if (parameters.terminate_with_first_full_solution && solution.complete) {
      break;
    }
    const auto next_time = Clock::now();
    if (next_time >= solve_timeout_time) {
      statistics.finished = false;
      statistics.search_time_ms +=
          std::chrono::duration_cast<std::chrono::milliseconds>(
              next_time - solve_start_time)
              .count();
      time_added = true;
      break;
    }
  }

  if (!time_added) {
    statistics.search_time_ms +=
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - solve_start_time)
            .count();
  }
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
