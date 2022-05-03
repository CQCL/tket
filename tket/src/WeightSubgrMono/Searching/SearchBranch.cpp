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

#include "WeightSubgrMono/Searching/SearchBranch.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Reducing/DistancesReducer.hpp"
#include "WeightSubgrMono/WeightPruning/WeightChecker.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SearchBranch::SearchBranch(
    const PossibleAssignments& initial_pattern_v_to_possible_target_v,
    const NeighboursData& pattern_ndata, NearNeighboursData& pattern_near_ndata,
    const NeighboursData& target_ndata, NearNeighboursData& target_near_ndata,
    unsigned max_distance_reduction_value)
    : m_pattern_ndata(pattern_ndata),
      m_target_ndata(target_ndata),
      m_derived_graphs_reducer(m_pattern_ndata, m_target_ndata),
      m_neighbours_reducer(m_pattern_ndata, m_target_ndata),
      m_nodes_raw_data_wrapper(initial_pattern_v_to_possible_target_v),
      m_domains_accessor(m_nodes_raw_data_wrapper),
      m_node_list_traversal(m_nodes_raw_data_wrapper) {
  // In what order should we do reduction/checks?
  // The simplest/cheapest first? Most powerful?
  // Seems a difficult question...
  m_reducer_wrappers.emplace_back(m_neighbours_reducer);
  m_reducer_wrappers.emplace_back(m_derived_graphs_reducer);

  if (max_distance_reduction_value >= 2) {
    m_distance_reducers.reserve(max_distance_reduction_value - 1);
    m_reducer_wrappers.reserve(max_distance_reduction_value - 1);
    for (unsigned distance = 2; distance <= max_distance_reduction_value;
         ++distance) {
      m_distance_reducers.emplace_back(
          m_pattern_ndata, pattern_near_ndata, m_target_ndata,
          target_near_ndata, distance);
    }
    // Now that all the reducer objects are stored
    // (the vector will not be resized),
    // we can safely take references to them.
    for (auto& reducer : m_distance_reducers) {
      m_reducer_wrappers.emplace_back(reducer);
    }
  }
}

const DomainsAccessor& SearchBranch::get_domains_accessor() const {
  return m_domains_accessor;
}
DomainsAccessor& SearchBranch::get_domains_accessor_nonconst() {
  return m_domains_accessor;
}

std::set<VertexWSM> SearchBranch::get_used_target_vertices() const {
  return m_node_list_traversal.get_used_target_vertices();
}

void SearchBranch::activate_weight_checker(WeightWSM total_p_edge_weights) {
  m_weight_checker_ptr = std::make_unique<WeightChecker>(
      m_pattern_ndata, m_target_ndata, *this, total_p_edge_weights);
  TKET_ASSERT(m_weight_checker_ptr);
}

bool SearchBranch::perform_single_assignment_checks_in_reduce_loop(
    std::size_t num_assignments_alldiff_processed) {
  const auto& new_assignments = m_domains_accessor.get_new_assignments();
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
        // invalid, so we remove it completely from all data.
        m_node_list_traversal.erase_impossible_assignment(
            new_assignments[ii], NodeListTraversal::ImpossibleAssignmentAction::
                                     PROCESS_CURRENT_NODE);
        // The current node is a nogood, so we'll move up from it shortly;
        // no point in further checks or reductions, this node is doomed!
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
      num_assignments_alldiff_processed, parameters.max_weight,
      m_domains_accessor.get_candidate_vertices_for_assignment_nonconst());

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
    const WeightWSM max_extra_scalar_product =
        parameters.max_weight - scalar_product;

    auto weight_check_result = m_weight_checker_ptr->operator()(
        m_domains_accessor, max_extra_scalar_product);

    if (weight_check_result.invalid_t_vertex) {
      // This is rare. Crudely treat as a list of assignments.
      // We want to pass them all in, though, even if we're at a nogood;
      // they might not be detected again.
      const VertexWSM invalid_tv = weight_check_result.invalid_t_vertex.value();
      auto action = weight_check_result.nogood
                        ? NodeListTraversal::ImpossibleAssignmentAction::
                              PROCESS_CURRENT_NODE
                        : NodeListTraversal::ImpossibleAssignmentAction::
                              IGNORE_CURRENT_NODE;

      for (VertexWSM pv : m_domains_accessor.get_pattern_vertices()) {
        if (!m_node_list_traversal.erase_impossible_assignment(
                std::make_pair(pv, invalid_tv), action)) {
          // We've found the current node to be invalid.
          // No point in checking the node further. HOWEVER we don't return yet,
          // so that TV is completely removed from all nodes.
          weight_check_result.nogood = true;
          action = NodeListTraversal::ImpossibleAssignmentAction::
              PROCESS_CURRENT_NODE;
        }
      }
    }
    if (weight_check_result.nogood) {
      return false;
    }
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
  const auto& new_assignments = m_domains_accessor.get_new_assignments();

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
    const auto hall_set_result = m_hall_set_reduction.reduce(
        m_domains_accessor, m_work_set_for_reducers);
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
  // Avoid the candidate vertices building up too much.
  // Will be useful later (when we need to choose a new assignment to make).
  auto& candidate_vertices =
      m_domains_accessor.get_candidate_vertices_for_assignment_nonconst();
  for (const auto& entry : m_domains_accessor.get_new_assignments()) {
    candidate_vertices.erase(entry.first);
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
