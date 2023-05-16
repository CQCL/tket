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

#include "tkwsm/InitPlacement/MonteCarloManager.hpp"

#include "tkwsm/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

MonteCarloManager::MonteCarloManager(
    const MonteCarloManagerParameters& parameters)
    : m_parameters(parameters),
      m_runs_without_record_breaking(0),
      m_runs_without_progress(0),
      m_next_iteration_to_reset_if_no_progress(0),
      m_next_iteration_to_reset_if_no_record_breaker(0) {
  set_maximum(m_best_cost);
}

void MonteCarloManager::update_after_reset(unsigned iteration) {
  unsigned max_extra_iters_without_record_breaker =
      (iteration *
       m_parameters
           .per_kilo_fraction_of_allowed_extra_iterations_without_record_breakers) /
      1024;

  m_next_iteration_to_reset_if_no_record_breaker =
      iteration + std::max(
                      max_extra_iters_without_record_breaker,
                      m_parameters.min_iterations_for_change);
  update_after_weak_progress(iteration);
}

void MonteCarloManager::update_after_weak_progress(unsigned iteration) {
  unsigned max_extra_iters_for_no_progress =
      (iteration *
       m_parameters
           .per_kilo_fraction_of_allowed_extra_iterations_without_weak_progress) /
      1024;

  m_next_iteration_to_reset_if_no_progress =
      iteration + std::max(
                      max_extra_iters_for_no_progress,
                      m_parameters.min_iterations_for_change);
}

MonteCarloManager::Action MonteCarloManager::register_progress(
    WeightWSM new_cost, unsigned iteration) {
  m_runs_without_progress = 0;
  if (new_cost < m_best_cost) {
    // It's a record breaker! Treat the same as a reset.
    m_runs_without_record_breaking = 0;
    update_after_reset(iteration);
    m_best_cost = new_cost;
    return Action::CONTINUE_WITH_CURRENT_SOLUTION;
  }
  // It's not a record breaker, but it is still progress.
  if (iteration < m_next_iteration_to_reset_if_no_record_breaker) {
    update_after_weak_progress(iteration);
    return Action::CONTINUE_WITH_CURRENT_SOLUTION;
  }
  // Too long without a new record. Reset!
  ++m_runs_without_record_breaking;
  if (m_runs_without_record_breaking >
      m_parameters.max_runs_without_record_breaking) {
    return Action::TERMINATE;
  }
  update_after_reset(iteration);
  return Action::RESET_TO_NEW_SOLUTION;
}

MonteCarloManager::Action MonteCarloManager::register_failure(
    unsigned iteration) {
  if (iteration >= m_next_iteration_to_reset_if_no_record_breaker ||
      iteration >= m_next_iteration_to_reset_if_no_progress) {
    // We're going to reset! (Or terminate).
    ++m_runs_without_progress;
    ++m_runs_without_record_breaking;
    if (m_runs_without_progress > m_parameters.max_runs_without_progress ||
        m_runs_without_record_breaking >
            m_parameters.max_runs_without_record_breaking) {
      return Action::TERMINATE;
    }
    update_after_reset(iteration);
    return Action::RESET_TO_NEW_SOLUTION;
  }
  return Action::CONTINUE_WITH_CURRENT_SOLUTION;
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
