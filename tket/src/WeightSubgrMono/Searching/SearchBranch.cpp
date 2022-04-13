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
    PossibleAssignments initial_pattern_v_to_possible_target_v,
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    DistancesReducer& distances_reducer)
    : m_pattern_ndata(pattern_ndata),
      m_target_ndata(target_ndata),
      m_distances_reducer(distances_reducer),
      m_derived_graphs_reducer(m_pattern_ndata, m_target_ndata),
      m_assignment_checker(m_derived_graphs_reducer, m_distances_reducer),
      m_enriched_nodes_index(0) {
  m_enriched_nodes.resize(1);
  m_enriched_nodes[0].node.set_possible_assignments(
      initial_pattern_v_to_possible_target_v);
  m_enriched_nodes[0].superficially_valid = true;
}

std::set<VertexWSM> SearchBranch::get_used_target_vertices() const {
  std::set<VertexWSM> vertices;
  for (std::size_t index = 0; index <= m_enriched_nodes_index; ++index) {
    const auto& domains =
        m_enriched_nodes[index].node.get_possible_assignments();
    for (const auto& entry : domains) {
      for (VertexWSM tv : entry.second) {
        vertices.insert(tv);
      }
    }
  }
  return vertices;
}

void SearchBranch::activate_weight_checker(WeightWSM total_p_edge_weights) {
  m_weight_checker_ptr = std::make_unique<WeightChecker>(
      m_pattern_ndata, m_target_ndata, *this, total_p_edge_weights);
  TKET_ASSERT(m_weight_checker_ptr);
}

bool SearchBranch::attempt_to_clear_outstanding_impossible_data() {
  // Invalid single TVs are very rare;
  // do a simple hack of adding them to assignments.
  while (!m_outstanding_impossible_target_vertices.empty()) {
    const VertexWSM tv = m_outstanding_impossible_target_vertices.back();
    m_outstanding_impossible_target_vertices.pop_back();
    for (std::size_t index = 0; index <= m_enriched_nodes_index; ++index) {
      for (const auto& entry :
           m_enriched_nodes[index].node.get_possible_assignments()) {
        if (entry.second.count(tv) != 0) {
          m_outstanding_impossible_assignments.emplace_back(entry.first, tv);
        }
      }
    }
  }

  auto& node = get_current_node_nonconst();

  // DON'T increment ii, as we erase by swapping with the back.
  for (unsigned ii = 0; ii < m_outstanding_impossible_assignments.size();) {
    auto& impossible_assignment = m_outstanding_impossible_assignments[ii];
    const auto erasure_result = node.erase_assignment(impossible_assignment);
    if (!erasure_result.valid) {
      return false;
    }
    // Erase from previous nodes also, IF safe to do so.
    bool was_erased = true;
    for (unsigned index = 0; index < m_enriched_nodes_index; ++index) {
      const auto erase_attempt_result =
          m_enriched_nodes[index].node.attempt_to_erase_assignment(
              impossible_assignment);

      if (erase_attempt_result == NodeWSM::SoftErasureResult::TV_REMAINS) {
        was_erased = false;
        break;
      }
    }
    if (was_erased) {
      // We can discard it; swap with the back!
      if (ii + 1 < m_outstanding_impossible_assignments.size()) {
        impossible_assignment = m_outstanding_impossible_assignments.back();
      }
      m_outstanding_impossible_assignments.pop_back();
    } else {
      // We cannot erase yet.
      ++ii;
    }
  }
  return true;
}


bool SearchBranch::reduce_current_node(const ReductionParameters& parameters) {
  if (!attempt_to_clear_outstanding_impossible_data()) {
    return false;
  }
  m_distances_reducer.reset(parameters.max_distance_reduction_value);
  m_derived_graphs_reducer.clear();

  auto& node = get_current_node_nonconst();
  const auto& new_assignments = node.get_new_assignments();

  // Within the loop, we reduce; we break out when it stops changing,
  // so is valid and fully reduced.

  // At the start of each new loop, this tells us how many initial elements
  // of new_assignments have been processed by the alldiff reducer,
  // and the standard non-altering filter/checking routines.
  std::size_t n_assignments_processed = 0;

  m_hall_set_reducer.clear();

  for (;;) {
    if (!node.alldiff_reduce(n_assignments_processed)) {
      return false;
    }
    // First, do all checks/updates which don't alter domains.

    for (auto ii = n_assignments_processed; ii < new_assignments.size(); ++ii) {
      if (!m_assignment_checker(
              new_assignments[ii], parameters.max_distance_reduction_value)) {
        register_impossible_assignment(new_assignments[ii]);
        return false;
      }
    }
    {
      const auto weight_result_opt = m_weight_updater(
          m_pattern_ndata, m_target_ndata, node.get_possible_assignments(),
          new_assignments, n_assignments_processed, node.get_scalar_product(),
          parameters.max_weight,
          m_enriched_nodes[m_enriched_nodes_index]
              .pvs_adjacent_to_newly_assigned_vertices);

      if (!weight_result_opt) {
        return false;
      }
      node.set_scalar_product(weight_result_opt->scalar_product);
      node.set_total_pattern_edge_weights(
          node.get_total_pattern_edge_weights() +
          weight_result_opt->total_extra_p_edge_weights);
    }
    {
      if (!is_maximum(parameters.max_weight)) {
        const WeightWSM current_scalar_product = node.get_scalar_product();
        if (current_scalar_product > parameters.max_weight) {
          return false;
        }
        const WeightWSM max_extra_scalar_product =
            parameters.max_weight - current_scalar_product;
        if (m_weight_checker_ptr) {
          const auto weight_check_result =
              m_weight_checker_ptr->operator()(node, max_extra_scalar_product);
          if (weight_check_result.invalid_t_vertex) {
            m_outstanding_impossible_target_vertices.emplace_back(
                weight_check_result.invalid_t_vertex.value());
          }
          if (weight_check_result.nogood) {
            return false;
          }
        }
      }
    }

    // Now begin REDUCTIONS, which do alter domains.
    // If we detect that a new assignment has occurred,
    // we jump back to the start of the loop. (Since the node reductions
    // and assignment checking should be faster
    // than these more complicated reductions).

    n_assignments_processed = new_assignments.size();
    {
      const auto distance_reduction_result = m_distances_reducer(node);

      if (distance_reduction_result.impossible_assignment) {
        register_impossible_assignment(
            distance_reduction_result.impossible_assignment.value());
        return false;
      }

      if (distance_reduction_result.nogood_found) {
        return false;
      }
      if (distance_reduction_result.new_assignments_created) {
        // Jump back to start of loop.
        continue;
      }
    }
    {
      const auto hall_set_result = m_hall_set_reducer(node);
      if (hall_set_result == HallSetReducer::Result::NEW_ASSIGNMENTS) {
        // Jump back to start of loop.
        continue;
      }
      if (hall_set_result == HallSetReducer::Result::NOGOOD) {
        return false;
      }
    }
    {
      const auto derived_graphs_result = m_derived_graphs_reducer.reduce(node);
      if (derived_graphs_result ==
          DerivedGraphsReducer::ReductionResult::NEW_ASSIGNMENT) {
        continue;
      }
      if (derived_graphs_result ==
          DerivedGraphsReducer::ReductionResult::NOGOOD) {
        return false;
      }
      TKET_ASSERT(
          derived_graphs_result ==
          DerivedGraphsReducer::ReductionResult::SUCCESS);
    }

    if (n_assignments_processed == new_assignments.size()) {
      // No further assignments from reductions.
      break;
    }
    TKET_ASSERT(false);
  }

  // Avoid the candidate vertices building up too much.
  // Will be useful later (when we need to choose a new assignment to make).
  for (const auto& entry : new_assignments) {
    m_enriched_nodes[m_enriched_nodes_index]
        .pvs_adjacent_to_newly_assigned_vertices.erase(entry.first);
  }

  // Now that we've processed all the assignments,
  // we don't need them.
  node.clear_new_assignments();
  return true;
}

bool SearchBranch::backtrack(const ReductionParameters& parameters) {
  do {
    if (m_enriched_nodes_index == 0) {
      return false;
    }
    --m_enriched_nodes_index;
  } while (!m_enriched_nodes[m_enriched_nodes_index].superficially_valid ||
           !reduce_current_node(parameters));
  return true;
}

void SearchBranch::move_down(VertexWSM p_vertex, VertexWSM t_vertex) {
  // We shouldn't be moving down unless we're FULLY reduced.
  TKET_ASSERT(get_current_node().get_new_assignments().empty());

  ++m_enriched_nodes_index;
  if (m_enriched_nodes_index >= m_enriched_nodes.size()) {
    m_enriched_nodes.resize(m_enriched_nodes_index + 1);
  }
  auto& prev_enriched_node = m_enriched_nodes.at(m_enriched_nodes_index - 1);
  auto& current_enriched_node = m_enriched_nodes.at(m_enriched_nodes_index);

  // Remove the assignment from the previous node.
  const auto assignment = std::make_pair(p_vertex, t_vertex);
  const auto erasure_result =
      prev_enriched_node.node.erase_assignment(assignment);

  prev_enriched_node.superficially_valid = erasure_result.valid;
  TKET_ASSERT(erasure_result.assignment_was_possible);

  // Make the new node, with the new assignment.
  current_enriched_node.node = prev_enriched_node.node;
  current_enriched_node.node.clear_new_assignments();
  current_enriched_node.node.force_assignment(assignment);

  // Deal with the enriched data not in "node".
  current_enriched_node.superficially_valid = true;
  current_enriched_node.pvs_adjacent_to_newly_assigned_vertices =
      prev_enriched_node.pvs_adjacent_to_newly_assigned_vertices;
}

void SearchBranch::register_impossible_assignment(
    const std::pair<VertexWSM, VertexWSM>& assignment) {
  m_outstanding_impossible_assignments.emplace_back(assignment);
}

const NodeWSM& SearchBranch::get_current_node() const {
  return m_enriched_nodes.at(m_enriched_nodes_index).node;
}
NodeWSM& SearchBranch::get_current_node_nonconst() {
  return m_enriched_nodes.at(m_enriched_nodes_index).node;
}
std::set<VertexWSM>& SearchBranch::get_current_node_candidate_variables() {
  return m_enriched_nodes.at(m_enriched_nodes_index)
      .pvs_adjacent_to_newly_assigned_vertices;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
