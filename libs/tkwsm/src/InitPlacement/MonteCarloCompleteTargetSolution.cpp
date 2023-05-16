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

#include "tkwsm/InitPlacement/MonteCarloCompleteTargetSolution.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/InitPlacement/UtilsIQP.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

MonteCarloCompleteTargetSolution::MonteCarloCompleteTargetSolution(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    // Any target edge not mentioned in target_ndata is given this weight.
    WeightWSM implicit_target_weight,
    // If left unspecified, i.e. set to zero,
    // will choose a reasonable default.
    unsigned max_iterations)
    : m_implicit_target_weight(implicit_target_weight),

      m_iterations(0),
      m_max_iterations(max_iterations),
      m_solution_jumper(pattern_ndata, target_ndata, implicit_target_weight) {
  const unsigned number_of_pv =
      pattern_ndata.get_number_of_nonisolated_vertices();
  if (m_max_iterations == 0) {
    // Nothing special about these magic numbers; need to experiment!
    m_max_iterations =
        1000 + 100 * (number_of_pv +
                      target_ndata.get_number_of_nonisolated_vertices());
  }

  m_random_bits_and_tv.resize(
      target_ndata.get_number_of_nonisolated_vertices());
  for (unsigned tv = 0; tv < m_random_bits_and_tv.size(); ++tv) {
    m_random_bits_and_tv[tv].second = tv;
  }
  set_maximum(m_best_scalar_product);
  m_number_of_random_bits = 20;
  if (target_ndata.get_number_of_nonisolated_vertices() > 1000) {
    m_number_of_random_bits = 30;
  }

  reset_target_vertices();
  const unsigned max_pv = number_of_pv - 1;

  for (; m_iterations < m_max_iterations; ++m_iterations) {
    const unsigned pv = m_rng.get_size_t(max_pv);
    const unsigned tv = m_rng.get_element(m_random_bits_and_tv).second;

    const auto scalar_product_decrease_opt =
        m_solution_jumper.perform_move_and_get_scalar_product_decrease(
            pv, tv, 1);

    MonteCarloManager::Action action;

    if (scalar_product_decrease_opt) {
      TKET_ASSERT(
          scalar_product_decrease_opt.value() <= m_current_scalar_product);
      m_current_scalar_product -= scalar_product_decrease_opt.value();

      TKET_ASSERT(
          m_current_scalar_product == get_scalar_product_with_complete_target(
                                          pattern_ndata, target_ndata,
                                          implicit_target_weight,
                                          m_solution_jumper.get_assignments()));

      new_solution_is_record_breaker();
      action =
          m_manager.register_progress(m_current_scalar_product, m_iterations);
    } else {
      // No decrease yet.
      action = m_manager.register_failure(m_iterations);
    }
    if (action == MonteCarloManager::Action::TERMINATE) {
      return;
    }
    if (action == MonteCarloManager::Action::CONTINUE_WITH_CURRENT_SOLUTION) {
      continue;
    }
    TKET_ASSERT(action == MonteCarloManager::Action::RESET_TO_NEW_SOLUTION);
    reset_target_vertices();
  }
}

void MonteCarloCompleteTargetSolution::reset_target_vertices() {
  for (auto& entry : m_random_bits_and_tv) {
    entry.first =
        m_fast_random_bits.get_random_bits(m_rng, m_number_of_random_bits);
  }

  auto& assignments = m_solution_jumper.get_assignments_to_overwrite();
  std::sort(m_random_bits_and_tv.begin(), m_random_bits_and_tv.end());

  for (unsigned pv = 0; pv < assignments.size(); ++pv) {
    assignments[pv] = m_random_bits_and_tv[pv].second;
  }
  m_current_scalar_product =
      m_solution_jumper.reset_and_get_new_scalar_product();

  TKET_ASSERT(
      m_current_scalar_product == get_scalar_product_with_complete_target(
                                      m_solution_jumper.get_pattern_ndata(),
                                      m_solution_jumper.get_target_ndata(),
                                      m_implicit_target_weight, assignments));

  // Check if we've got a new best solution and update.
  new_solution_is_record_breaker();
}

bool MonteCarloCompleteTargetSolution::new_solution_is_record_breaker() {
  if (m_current_scalar_product >= m_best_scalar_product) {
    return false;
  }
  m_best_scalar_product = m_current_scalar_product;
  m_best_assignments = m_solution_jumper.get_assignments();
  return true;
}

const std::vector<unsigned>&
MonteCarloCompleteTargetSolution::get_best_assignments() const {
  return m_best_assignments;
}

WeightWSM MonteCarloCompleteTargetSolution::get_best_scalar_product() const {
  return m_best_scalar_product;
}

unsigned MonteCarloCompleteTargetSolution::iterations() const {
  return m_iterations;
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
