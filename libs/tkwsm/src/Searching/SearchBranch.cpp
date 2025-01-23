// Copyright Quantinuum
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

#include "tkwsm/Searching/SearchBranch.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/EndToEndWrappers/SolutionData.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Reducing/DistancesReducer.hpp"
#include "tkwsm/WeightPruning/WeightChecker.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SearchBranch::SearchBranch(
    const DomainInitialiser::InitialDomains& initial_domains,
    const NeighboursData& pattern_ndata, NearNeighboursData& pattern_near_ndata,
    const NeighboursData& target_ndata, NearNeighboursData& target_near_ndata,
    unsigned max_distance_reduction_value, ExtraStatistics& extra_statistics)
    : m_pattern_ndata(pattern_ndata),
      m_target_ndata(target_ndata),
      m_extra_statistics(extra_statistics),
      m_derived_graphs_reducer(m_pattern_ndata, m_target_ndata),
      m_nodes_raw_data_wrapper(
          initial_domains, m_target_ndata.get_number_of_nonisolated_vertices()),
      m_domains_accessor(m_nodes_raw_data_wrapper),
      m_node_list_traversal(m_nodes_raw_data_wrapper) {
  m_extra_statistics.number_of_pattern_vertices =
      m_pattern_ndata.get_number_of_nonisolated_vertices();
  m_extra_statistics.number_of_target_vertices =
      m_target_ndata.get_number_of_nonisolated_vertices();

  m_extra_statistics.initial_number_of_possible_assignments = 0;
  for (const boost::dynamic_bitset<>& domain : initial_domains) {
    m_extra_statistics.initial_number_of_possible_assignments += domain.count();
  }

  // In what order should we do reduction/checks?
  // The simplest/cheapest first? Most powerful?
  // Seems a difficult question...
  m_reducer_wrappers.reserve(max_distance_reduction_value + 1);
  /// m_reducer_wrappers.emplace_back(m_neighbours_reducer);
  m_reducer_wrappers.emplace_back(m_derived_graphs_reducer);
  m_distance_reducers.reserve(max_distance_reduction_value);

  for (unsigned distance = 1; distance <= max_distance_reduction_value;
       ++distance) {
    m_distance_reducers.emplace_back(
        pattern_near_ndata, target_near_ndata, distance);
  }
  // Now that all the reducer objects are stored
  // (the vector will not be resized),
  // we can safely take references to them.
  for (DistancesReducer& reducer : m_distance_reducers) {
    m_reducer_wrappers.emplace_back(reducer);
  }
}

const DomainsAccessor& SearchBranch::get_domains_accessor() const {
  return m_domains_accessor;
}
DomainsAccessor& SearchBranch::get_domains_accessor_nonconst() {
  return m_domains_accessor;
}

boost::dynamic_bitset<> SearchBranch::get_used_target_vertices() const {
  return m_node_list_traversal.get_used_target_vertices();
}

const ExtraStatistics& SearchBranch::get_updated_extra_statistics() {
  if (m_weight_checker_ptr) {
    const auto weight_checker_data_opt =
        m_weight_checker_ptr->get_tv_data_opt();
    if (weight_checker_data_opt) {
      m_extra_statistics.n_tv_initially_passed_to_weight_nogood_detector =
          weight_checker_data_opt->initial_number_of_tv;
      m_extra_statistics.n_tv_still_valid_in_weight_nogood_detector =
          weight_checker_data_opt->final_number_of_tv;
    }
  }
  m_extra_statistics.total_number_of_assignments_tried = 0;
  for (const auto& entry : m_checked_assignments) {
    m_extra_statistics.total_number_of_assignments_tried += entry.second.size();
  }
  return m_extra_statistics;
}

void SearchBranch::activate_weight_checker(WeightWSM total_p_edge_weights) {
  m_weight_checker_ptr = std::make_unique<WeightChecker>(
      m_pattern_ndata, m_target_ndata, *this, total_p_edge_weights,
      m_impossible_target_vertices);
  TKET_ASSERT(m_weight_checker_ptr);
}

bool SearchBranch::perform_single_assignment_checks_in_reduce_loop(
    std::size_t num_assignments_alldiff_processed) {
  const std::vector<std::pair<VertexWSM, VertexWSM>>& new_assignments =
      m_domains_accessor.get_new_assignments();
  for (auto ii = num_assignments_alldiff_processed; ii < new_assignments.size();
       ++ii) {
    if (m_checked_assignments[new_assignments[ii].first].count(
            new_assignments[ii].second) != 0) {
      continue;
    }
    for (auto& reducer_wrapper : m_reducer_wrappers) {
      if (!reducer_wrapper.check(new_assignments[ii])) {
        // The whole point is, the check for PV->TV does NOT depend on
        // the other domains; since it's failed the check, it was ALWAYS
        // invalid, so we remove it completely from ALL data.
        m_node_list_traversal.erase_impossible_assignment(new_assignments[ii]);
        // The current node is a nogood, so we'll move up from it shortly;
        // no point in further checks or reductions, this node is doomed!
        ++m_extra_statistics.total_number_of_impossible_assignments;
        return false;
      }
    }
  }
  return true;
}

bool SearchBranch::check_and_update_scalar_product_in_reduce_loop(
    const ReductionParameters& parameters,
    std::size_t num_assignments_alldiff_processed) {
  const auto weight_result_opt = m_weight_calculator(
      m_pattern_ndata, m_target_ndata, m_domains_accessor,
      num_assignments_alldiff_processed, parameters.max_weight);

  if (!weight_result_opt) {
    return false;
  }
  m_domains_accessor.set_scalar_product(weight_result_opt->scalar_product)
      .set_total_p_edge_weights(
          weight_result_opt->total_extra_p_edge_weights +
          m_domains_accessor.get_total_p_edge_weights());
  return true;
}

bool SearchBranch::perform_weight_nogood_check_in_reduce_loop(
    const ReductionParameters& parameters) {
  if (is_maximum(parameters.max_weight)) {
    // No scalar product constraint, so nothing to check.
    return true;
  }
  const auto scalar_product = m_domains_accessor.get_scalar_product();
  if (scalar_product > parameters.max_weight) {
    return false;
  }

  if (m_weight_checker_ptr) {
    m_impossible_target_vertices.clear();

    const WeightWSM max_extra_scalar_product =
        parameters.max_weight - scalar_product;

    bool node_is_valid = m_weight_checker_ptr->check(
        m_domains_accessor, max_extra_scalar_product);

    const auto number_of_pv =
        m_domains_accessor.get_number_of_pattern_vertices();

    for (VertexWSM tv : m_impossible_target_vertices) {
      // This is rare. Crudely treat as a list of assignments.
      // We want to pass them all in now, though, even if we're at a nogood;
      // they might not be detected again.
      m_extra_statistics.impossible_target_vertices.emplace_back(tv);
      for (unsigned pv = 0; pv < number_of_pv; ++pv) {
        m_node_list_traversal.erase_impossible_assignment(
            std::make_pair(pv, tv));
      }
      // Detecting an invalid TV is separate from whether
      // the current node is a nogood.
      node_is_valid &= m_domains_accessor.current_node_is_valid();
    }
    return node_is_valid;
  }
  return true;
}

ReductionResult SearchBranch::perform_reducers_in_reduce_loop() {
  for (auto& reducer_wrapper : m_reducer_wrappers) {
    const auto reduction_result =
        reducer_wrapper.reduce(m_domains_accessor, m_work_set_for_reducers);
    if (reduction_result != ReductionResult::SUCCESS) {
      return reduction_result;
    }
  }
  return ReductionResult::SUCCESS;
}

bool SearchBranch::perform_main_reduce_loop(
    const ReductionParameters& parameters) {
  // Within the loop, we reduce; we break out when it stops changing,
  // so is valid and fully reduced.

  // At the start of each new loop, this tells us how many initial elements
  // of new_assignments have been processed by the alldiff reducer.
  std::size_t num_assignments_alldiff_processed = 0;

  // We will not change this directly. The reference will remain valid,
  // but it will change indirectly due to other actions (reducing, etc.).
  const std::vector<std::pair<VertexWSM, VertexWSM>>& new_assignments =
      m_domains_accessor.get_new_assignments();

  for (;;) {
    if (!(m_domains_accessor.alldiff_reduce_current_node(
              num_assignments_alldiff_processed) &&
          // Checks/updates which don't alter domains follow.
          perform_single_assignment_checks_in_reduce_loop(
              num_assignments_alldiff_processed) &&
          check_and_update_scalar_product_in_reduce_loop(
              parameters, num_assignments_alldiff_processed))) {
      return false;
    }

    num_assignments_alldiff_processed = new_assignments.size();

    if (!perform_weight_nogood_check_in_reduce_loop(parameters)) {
      return false;
    }
    // It's rare, but possible, that the weight nogood check found an impossible
    // TV, in which case it was erased, and the domains DID change...
    if (num_assignments_alldiff_processed != new_assignments.size()) {
      // Jump back to start!
      continue;
    }

    // Now do REDUCTIONS, which do alter domains.
    // The scalar product updates etc. have all been completed.
    // If we detect that a new assignment has occurred,
    // we jump back to the start of the loop. (Since the node reductions
    // and assignment checking should be faster
    // than these more complicated reductions).
    {
      const auto reducers_result = perform_reducers_in_reduce_loop();
      if (reducers_result == ReductionResult::NOGOOD) {
        return false;
      }
      if (reducers_result == ReductionResult::NEW_ASSIGNMENTS) {
        // Jump back to start of loop.
        TKET_ASSERT(num_assignments_alldiff_processed < new_assignments.size());
        continue;
      }
      TKET_ASSERT(num_assignments_alldiff_processed == new_assignments.size());
    }

    // It is reasonable that the final reduction is the Hall set reducer,
    // because it alone benefits from smaller domains
    // (i.e., is more likely to detect a possible reduction).
    // The standard reducers ONLY work with actual ASSIGNMENTS;
    // it makes no difference if a domain has size 2 or 100.
    //*
    const auto hall_set_result =
        m_hall_set_reduction.reduce(m_domains_accessor);
    if (hall_set_result == ReductionResult::NOGOOD) {
      return false;
    }

    TKET_ASSERT(num_assignments_alldiff_processed <= new_assignments.size());
    if (num_assignments_alldiff_processed == new_assignments.size()) {
      TKET_ASSERT(hall_set_result == ReductionResult::SUCCESS);
      // No further assignments from reductions.
      break;
    }
    // Will jump back to start of loop.
    TKET_ASSERT(hall_set_result == ReductionResult::NEW_ASSIGNMENTS);
    TKET_ASSERT(num_assignments_alldiff_processed < new_assignments.size());
  }
  return true;
}

bool SearchBranch::reduce_current_node(const ReductionParameters& parameters) {
  for (auto& reducer_wrapper : m_reducer_wrappers) {
    reducer_wrapper.clear();
  }
  m_hall_set_reduction.clear();

  if (!perform_main_reduce_loop(parameters)) {
    return false;
  }
  // Now that we've processed all the assignments,
  // we don't need them.
  m_domains_accessor.clear_new_assignments();
  return true;
}

bool SearchBranch::backtrack(const ReductionParameters& parameters) {
  do {
    if (!m_node_list_traversal.move_up()) {
      return false;
    }
  } while (!reduce_current_node(parameters));
  return true;
}

void SearchBranch::move_down(VertexWSM p_vertex, VertexWSM t_vertex) {
  m_node_list_traversal.move_down(p_vertex, t_vertex);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
