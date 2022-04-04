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

#include "WeightSubgrMono/EndToEndWrappers/MainSolver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <sstream>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/EndToEndWrappers/CheckedWeightBounds.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchBranch.hpp"
#include "WeightSubgrMono/Searching/SharedData.hpp"
#include "WeightSubgrMono/Searching/ValueOrdering.hpp"
#include "WeightSubgrMono/Searching/VariableOrdering.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

MainSolver::MainSolver(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges,
      const MainSolverParameters& parameters) 
    : MainSolver(pattern_edges, target_edges) {
  solve(parameters);
}


MainSolver::MainSolver(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges,
      const std::vector<std::pair<VertexWSM, VertexWSM>>&
          suggested_assignments) 
    : MainSolver(pattern_edges, target_edges) {
  do_one_solve_iteration_with_suggestion(suggested_assignments);
}

typedef std::chrono::steady_clock Clock;

const SolutionStatistics& MainSolver::get_solution_statistics() const {
  return m_data.statistics;
}


MainSolver::MainSolver(
    const GraphEdgeWeights& pattern_edges,
    const GraphEdgeWeights& target_edges) {
  const auto init_start = Clock::now();
  const auto init_result = m_data.initialise(pattern_edges, target_edges);

  if (init_result == ReductionResult::FAILURE) {
    // No solution is possible.
    m_data.statistics.finished = true;
  }

  if (init_result == ReductionResult::FINISHED) {
    // The initial domains and vertex filtering have led to
    // a unique solution, which is valid; and stored in the shared data.
    m_data.statistics.finished = true;
    TKET_ASSERT(m_data.shared_data_ptr);

    // TODO: Check for int overflow and warn...although doesn't actually matter,
    // since the solution is unique anyway...!

    const auto& unique_solution =
        m_data.shared_data_ptr->solution_storage.best_solution();

    m_data.statistics.trivial_weight_lower_bound =
        unique_solution.total_scalar_product_weight;
    m_data.statistics.trivial_weight_initial_upper_bound =
        unique_solution.total_scalar_product_weight;
  }

  if (init_result == ReductionResult::SUCCESS) {
    // So far as we know, a solution may be possible; we need to search.
    // But we haven't found one yet.
    m_data.statistics.finished = false;
    TKET_ASSERT(m_data.shared_data_ptr);

    // The cheapest check.
    const CheckedWeightBounds checked_bounds(
        m_data.fixed_data, pattern_edges, target_edges, 10);

    if (!checked_bounds.other_inconsistency_occurred &&
        checked_bounds.lower_bound && checked_bounds.upper_bound &&
        checked_bounds.lower_bound.value() <=
            checked_bounds.upper_bound.value()) {
      m_data.statistics.trivial_weight_lower_bound =
          checked_bounds.lower_bound.value();
      m_data.statistics.trivial_weight_initial_upper_bound =
          checked_bounds.upper_bound.value();
    } else {
      // We're finished, but with no solutions.
      m_data.statistics.finished = true;
      m_data.shared_data_ptr.reset();
    }
  }

  m_data.statistics.initialisation_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          Clock::now() - init_start)
          .count();
}

void MainSolver::do_one_solve_iteration_with_suggestion(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& suggested_assignments) {
  if (m_data.statistics.finished) {
    return;
  }
  m_data.do_one_solve_iteration_with_suggestion(suggested_assignments);
}

void MainSolver::solve(const MainSolverParameters& parameters) {
  if (m_data.statistics.finished) {
    return;
  }

  TKET_ASSERT(m_data.initialised);
  TKET_ASSERT(m_data.shared_data_ptr);
  if (parameters.weight_upper_bound_constraint) {
    const WeightWSM weight_constraint =
        parameters.weight_upper_bound_constraint.value();
    if (m_data.previous_upper_bound_constraint) {
      TKET_ASSERT(
          m_data.previous_upper_bound_constraint.value() >= weight_constraint);
      m_data.previous_upper_bound_constraint = weight_constraint;
    }
    m_data.shared_data_ptr->solution_storage.set_pruning_weight(
        weight_constraint);
  }
  m_data.solve_loop_after_initialisation(parameters);
}

const SolutionWSM& MainSolver::get_best_solution() const {
  if (m_data.shared_data_ptr) {
    return m_data.shared_data_ptr->solution_storage.best_solution();
  }
  return m_data.empty_solution;
}

const MainSolver::FullSolutionsList& MainSolver::get_some_full_solutions()
    const {
  TKET_ASSERT(m_data.shared_data_ptr);
  return m_data.shared_data_ptr->solution_storage.get_some_full_solutions();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
