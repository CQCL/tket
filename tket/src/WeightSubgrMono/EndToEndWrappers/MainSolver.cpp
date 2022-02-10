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

typedef std::chrono::steady_clock Clock;

MainSolver::Parameters::Parameters()
    : timeout_ms(5000), terminate_with_first_full_solution(false) {
  set_maximum(max_iterations);
}

MainSolver::Parameters::Parameters(long long t_out_ms) : Parameters() {
  timeout_ms = t_out_ms;
}

struct MainSolver::Impl {
  bool initialised;
  FixedData fixed_data;
  std::unique_ptr<SharedData> shared_data_ptr;
  SearchBranch branch;
  VariableOrdering var_ordering;
  ValueOrdering val_ordering;
  SolutionStatistics statistics;
  std::optional<WeightWSM> previous_upper_bound_constraint;

  // Hack: in case we finish upon initialisation; an empty solution.
  // Otherwise, use the one stored in shared_data_ptr.
  const SolutionWSM empty_solution;

  Impl() : initialised(false) {}

  ~Impl() {}

  // Sets "fixed_data" and, if successful, shared_data_ptr.
  ReductionResult initialise(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges);

  void solve(const MainSolver::Parameters& parameters);
};

ReductionResult MainSolver::Impl::initialise(
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

MainSolver::MainSolver() : m_pimpl(std::make_unique<Impl>()) {}

MainSolver::~MainSolver() {}

const SolutionStatistics& MainSolver::solve(long long timeout_ms) {
  Parameters parameters;
  parameters.timeout_ms = timeout_ms;
  return solve(parameters);
}

const SolutionStatistics& MainSolver::solve(
    const GraphEdgeWeights& pattern_edges, const GraphEdgeWeights& target_edges,
    long long timeout_ms) {
  Parameters parameters;
  parameters.timeout_ms = timeout_ms;
  return solve(pattern_edges, target_edges, parameters);
}

const SolutionStatistics& MainSolver::solve(
    const GraphEdgeWeights& pattern_edges, const GraphEdgeWeights& target_edges,
    const Parameters& parameters) {
  initialise(pattern_edges, target_edges);
  return solve(parameters);
}

const SolutionStatistics& MainSolver::initialise(
    const GraphEdgeWeights& pattern_edges,
    const GraphEdgeWeights& target_edges) {
  TKET_ASSERT(m_pimpl);

  const auto init_start = Clock::now();
  const auto init_result = m_pimpl->initialise(pattern_edges, target_edges);

  if (init_result == ReductionResult::FAILURE) {
    // No solution is possible.
    m_pimpl->shared_data_ptr.reset();
    m_pimpl->statistics.finished = true;
  }

  if (init_result == ReductionResult::FINISHED) {
    // The initial domains and vertex filtering have led to
    // a unique solution, which is valid; and stored in the shared data.
    m_pimpl->statistics.finished = true;
    TKET_ASSERT(m_pimpl->shared_data_ptr);

    // TODO: Check for int overflow and warn...although doesn't actually matter,
    // since the solution is unique anyway...!

    const auto& unique_solution =
        m_pimpl->shared_data_ptr->solution_storage.best_solution();

    m_pimpl->statistics.trivial_weight_lower_bound =
        unique_solution.total_scalar_product_weight;
    m_pimpl->statistics.trivial_weight_initial_upper_bound =
        unique_solution.total_scalar_product_weight;
  }

  if (init_result == ReductionResult::SUCCESS) {
    // So far as we know, a solution may be possible; we need to search.
    // But we haven't found one yet.
    m_pimpl->statistics.finished = false;
    TKET_ASSERT(m_pimpl->shared_data_ptr);

    // The cheapest check.
    const CheckedWeightBounds checked_bounds(
        m_pimpl->fixed_data, pattern_edges, target_edges, 10);

    if (!checked_bounds.other_inconsistency_occurred &&
        checked_bounds.lower_bound && checked_bounds.upper_bound &&
        checked_bounds.lower_bound.value() <=
            checked_bounds.upper_bound.value()) {
      m_pimpl->statistics.trivial_weight_lower_bound =
          checked_bounds.lower_bound.value();
      m_pimpl->statistics.trivial_weight_initial_upper_bound =
          checked_bounds.upper_bound.value();
    } else {
      // We're finished, but with no solutions.
      m_pimpl->statistics.finished = true;
      m_pimpl->shared_data_ptr.reset();
    }
  }

  m_pimpl->statistics.initialisation_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          Clock::now() - init_start)
          .count();

  return m_pimpl->statistics;
}

const SolutionStatistics& MainSolver::solve(
    const MainSolver::Parameters& parameters) {
  if (m_pimpl->statistics.finished) {
    return m_pimpl->statistics;
  }

  const auto solve_start_time = Clock::now();
  const auto solve_timeout_time =
      solve_start_time + std::chrono::milliseconds(parameters.timeout_ms);
  TKET_ASSERT(m_pimpl);
  TKET_ASSERT(m_pimpl->initialised);
  TKET_ASSERT(m_pimpl->shared_data_ptr);
  if (parameters.weight_upper_bound_constraint) {
    const WeightWSM weight_constraint =
        parameters.weight_upper_bound_constraint.value();
    if (m_pimpl->previous_upper_bound_constraint) {
      TKET_ASSERT(
          m_pimpl->previous_upper_bound_constraint.value() >=
          weight_constraint);
      m_pimpl->previous_upper_bound_constraint = weight_constraint;
    }
    m_pimpl->shared_data_ptr->solution_storage.set_pruning_weight(
        weight_constraint);
  }

  const auto& solution =
      m_pimpl->shared_data_ptr->solution_storage.best_solution();

  bool time_added = false;
  for (; m_pimpl->statistics.iterations < parameters.max_iterations;
       m_pimpl->statistics.iterations += 1) {
    const bool search_finished = !m_pimpl->shared_data_ptr->search(
        m_pimpl->branch, m_pimpl->var_ordering, m_pimpl->val_ordering);

    if (search_finished) {
      m_pimpl->statistics.finished = true;
      break;
    }
    if (parameters.terminate_with_first_full_solution && solution.complete) {
      break;
    }
    const auto next_time = Clock::now();
    if (next_time >= solve_timeout_time) {
      m_pimpl->statistics.finished = false;
      m_pimpl->statistics.search_time_ms +=
          std::chrono::duration_cast<std::chrono::milliseconds>(
              next_time - solve_start_time)
              .count();
      time_added = true;
      break;
    }
  }
  if (!time_added) {
    m_pimpl->statistics.search_time_ms +=
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - solve_start_time)
            .count();
  }
  return m_pimpl->statistics;
}

const SolutionWSM& MainSolver::get_best_solution() const {
  TKET_ASSERT(m_pimpl);
  if (m_pimpl->shared_data_ptr) {
    return m_pimpl->shared_data_ptr->solution_storage.best_solution();
  }
  return m_pimpl->empty_solution;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
