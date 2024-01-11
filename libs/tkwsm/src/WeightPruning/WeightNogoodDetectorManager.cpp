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

#include "tkwsm/WeightPruning/WeightNogoodDetectorManager.hpp"

#include <algorithm>

#include "tkwsm/Common/DyadicFraction.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

WeightNogoodDetectorManager::WeightNogoodDetectorManager(
    WeightWSM total_p_edge_weights)
    : m_total_p_edge_weights(total_p_edge_weights) {}

template <typename T>
void clamp(T& val, T low, T high) {
  val = std::clamp(val, low, high);
}

namespace {
struct TempCounter {
  unsigned good_count = 0;
  unsigned ok_count = 0;
  unsigned bad_count = 0;
  unsigned old_total = 0;

  unsigned get_total() const { return good_count + ok_count + bad_count; }
};
}  // namespace

static TempCounter& get_counter() {
  static TempCounter counter;
  return counter;
}

bool WeightNogoodDetectorManager::should_activate_detector(
    WeightWSM current_weight, WeightWSM max_weight,
    WeightWSM current_sum_of_p_edge_weights, std::size_t n_assigned_vertices,
    std::size_t n_unassigned_vertices) {
  const auto total_n_vertices = n_assigned_vertices + n_unassigned_vertices;

  // Note that vertex numbers are sensibly small, so won't overflow.
  // Check state resetting.
  // (More detail: as we move up and down the search tree, nearby search nodes
  // are quite similar to each other, so our dynamic feedback strategy
  // of, roughly, "do more checks if it's successful a lot; reduce
  // the amount of checking if it's failing a lot" seems sensible.
  // But, if conditions change a lot, treat it as a fresh search).
  if (m_state.can_reset &&
      (n_assigned_vertices <= 2 ||
       1024 * n_assigned_vertices <=
           m_parameters.drop_below_n_assigned_vertices_reset_pk *
               total_n_vertices)) {
    // The number of assigned vertices is very low,
    // so we've backtracked far enough that we'll reset
    // (treat it as a whole new search).
    m_state = State();
  } else {
    if (!m_state.can_reset &&
        (n_unassigned_vertices <= 2 ||
         1024 * n_assigned_vertices >=
             m_parameters.rise_above_n_assigned_vertices_allow_reset_pk *
                 total_n_vertices)) {
      // We haven't reset for a while; we've now got enough vertices
      // (i.e., the search has progressed far enough that we'll
      // now allow reset once we backtrack far enough).
      m_state.can_reset = true;
    }
  }
  if (m_state.remaining_skips > 0) {
    --m_state.remaining_skips;
    return false;
  }

  // Are there enough assignments (and remaining unassigned vertices -
  // if you're very close to a solution, just do it,
  // rather than estimating)?
  if (n_assigned_vertices <= 2 || n_unassigned_vertices <= 2 ||
      current_weight == 0 || m_total_p_edge_weights == 0) {
    // Not enough assignments (or remaining unassigned vertices)
    // for a useful estimate.
    return false;
  }

  // Note that we are careful to avoid int oveflow.

  if (current_weight < max_weight / 1024 ||
      DyadicFraction(current_weight) <
          DyadicFraction(max_weight)
              .mult_n_over_k(m_state.min_weight_pk_to_activate)) {
    // We haven't got enough weight to start estimating.
    return false;
  }

  // If we very crudely estimate the total weight after all
  // assignments are made (simply based on the fraction of
  // pattern graph assignments so far:
  // either from number of vertices, or total p-weight),
  // is it large enough for a detection to be worthwhile
  // (because it has more chance of success, i.e. getting
  // a large estimate)?
  if (
      // The fraction of pattern graph assignments
      // based upon total p-weight.
      DyadicFraction(current_weight).mult(m_total_p_edge_weights) <
          DyadicFraction(max_weight)
              .mult(current_sum_of_p_edge_weights)
              .mult_n_over_k(m_state.final_weight_estimate_pk_to_activate) &&

      DyadicFraction(current_weight).mult(total_n_vertices) <
          DyadicFraction(max_weight)
              .mult(n_assigned_vertices)
              .mult_n_over_k(m_state.final_weight_estimate_pk_to_activate)) {
    return false;
  }
  return true;
}

void WeightNogoodDetectorManager::register_success() {
  get_counter().good_count += 1;
  m_state.final_weight_estimate_pk_to_activate *=
      m_parameters.final_weight_estimate_pk_success_growth_pk;
  m_state.final_weight_estimate_pk_to_activate /= 1024;

  m_state.min_weight_pk_to_activate *=
      m_parameters.success_pk_to_activate_growth_factor_pk;
  m_state.min_weight_pk_to_activate /= 1024;
  clamp_state_values();
  m_state.remaining_skips = 0;
}

void WeightNogoodDetectorManager::register_lower_bound_failure(
    WeightWSM current_weight, WeightWSM max_weight,
    WeightWSM extra_weight_lower_bound) {
  if (current_weight + 2 * extra_weight_lower_bound < max_weight) {
    // Bad failure.
    get_counter().bad_count += 1;
    m_state.final_weight_estimate_pk_to_activate *=
        m_parameters.final_weight_estimate_pk_bad_failure_growth_pk;
    m_state.min_weight_pk_to_activate *=
        m_parameters.bad_failure_pk_to_activate_growth_factor_pk;
    m_state.remaining_skips = m_parameters.bad_failure_skip;
  } else {
    get_counter().ok_count += 1;
    // "OK" failure; the lower bound wasn't that far off the amount needed.
    m_state.final_weight_estimate_pk_to_activate *=
        m_parameters.final_weight_estimate_pk_ok_failure_growth_pk;
    m_state.min_weight_pk_to_activate *=
        m_parameters.ok_failure_pk_to_activate_growth_factor_pk;
    m_state.remaining_skips = m_parameters.ok_failure_skip;
  }
  m_state.final_weight_estimate_pk_to_activate /= 1024;
  m_state.min_weight_pk_to_activate /= 1024;
  clamp_state_values();
}

void WeightNogoodDetectorManager::clamp_state_values() {
  clamp(
      m_state.final_weight_estimate_pk_to_activate,
      m_parameters.min_final_weight_estimate_pk,
      m_parameters.max_final_weight_estimate_pk);

  clamp(
      m_state.min_weight_pk_to_activate, m_parameters.min_weight_pk_to_activate,
      m_parameters.max_weight_pk_to_activate);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
