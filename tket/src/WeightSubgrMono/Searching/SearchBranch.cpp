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

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SharedData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SearchBranch::SearchBranch() : m_level(0) {}

ReductionResult SearchBranch::initialise(SharedData& shared_data) {
  m_level = 0;
  m_move_down_has_been_called = false;
  if (m_enriched_nodes.empty()) {
    m_enriched_nodes.resize(1);
  }
  auto& enriched_node = m_enriched_nodes[0];
  enriched_node.node_wrapper =
      SearchNodeWrapper(shared_data.fixed_data.initial_node);
  enriched_node.clear_enriched_data();
  m_assignments.clear();
  WeightWSM max_weight;
  set_maximum(max_weight);
  return reduce_current_node(shared_data, max_weight);
}

const Assignments& SearchBranch::get_assignments() const {
  return m_assignments;
}

Assignments& SearchBranch::get_assignments_mutable() { return m_assignments; }

bool SearchBranch::erase_assignment(VertexWSM pv, VertexWSM tv) {
  for (std::size_t level = 0; level <= m_level; ++level) {
    auto& domains = m_enriched_nodes[level]
                        .node_wrapper.get_mutable()
                        .pattern_v_to_possible_target_v;
    auto iter = domains.find(pv);
    if (iter == domains.end()) {
      return false;
    }
    auto& domain = iter->second;
    if (domain.size() <= 2) {
      return false;
    }
    auto tv_iter = domain.find(tv);
    if (tv_iter == domain.end()) {
      return false;
    }
    domain.erase(tv_iter);
  }
  return true;
}

ReductionResult SearchBranch::reduce_current_node(
    SharedData& shared_data, WeightWSM max_weight) {
  // E.g., we might have variable domains
  //  Dom(u) = Dom(v) = {a,b},
  // and later filters/reductions might reduce them to Dom(u) = Dom(v) = {a}.
  // We use this set to check for this; the assignments are not checked
  // when they occur.
  m_values_assigned_in_this_node.clear();

  // We will break out of this loop if and only if
  // we have finished (assigned all vertices);
  // otherwise, we shall return.
  for (;;) {
    auto& enriched_node = m_enriched_nodes[m_level];
    const auto assignments_processed_in_this_node =
        enriched_node.n_assignments_processed_by_all_diff_propagator;

    {
      const auto& chosen_assignments =
          enriched_node.node_wrapper.get().chosen_assignments;
      for (auto assignment_index = m_values_assigned_in_this_node.size();
           assignment_index < chosen_assignments.size(); ++assignment_index) {
        const auto& new_tv = chosen_assignments[assignment_index].second;
        if (!m_values_assigned_in_this_node.insert(new_tv).second) {
          // A duplicate TV.
          return ReductionResult::FAILURE;
        }

        const auto& pv = chosen_assignments[assignment_index].first;
        if (!shared_data.fixed_data.target_is_complete) {
          if (!shared_data.derived_graphs_filter.is_compatible(
                  pv, new_tv, shared_data.fixed_data)) {
            erase_assignment(pv, new_tv);
            return ReductionResult::FAILURE;
          }
        }
      }
    }
    if (!shared_data.fixed_data.alldiff_propagator.reduce(
            m_assignments, enriched_node.node_wrapper,
            enriched_node.n_assignments_processed_by_all_diff_propagator)) {
      return ReductionResult::FAILURE;
    }

    // We tried putting an extra derived_graphs_filter check here,
    // but it only slowed things down slightly;
    // this is not completely stupid, as m_assignments was smaller
    // when we checked above.

    if (!shared_data.fixed_data.weight_updater(
            shared_data.fixed_data, m_assignments, enriched_node.node_wrapper,
            assignments_processed_in_this_node, max_weight)) {
      return ReductionResult::FAILURE;
    }

    if (enriched_node.node_wrapper.get()
            .pattern_v_to_possible_target_v.empty()) {
      // If we are here, everything is at least CORRECT,
      // regardless of the other filters: all vertices have been assigned,
      // and all new edges checked.
      break;
    }

    // Now, more filtering.
    // QUESTION: in which ORDER should we apply different filters?
    // It should not affect correctness, but it very definitely
    // could affect speed.
    //
    // Very unclear...each filter potentially could reduce domain sizes,
    // and/or make new assignments, maybe enabling other filters
    // to reduce further...
    //
    // Estimates of both filter calculation time
    // and power would be helpful,
    // but even then, a bit unclear.

    // NOTE: the reducers may intersect current domains,
    // and thus reduce them. If reduced to empty, it's a nogood;
    // but if reduced to size 1, it's a new assignment.
    // HOWEVER, the assignments are not checked for all diff propagation;
    // that must be done above, at the start of this containing loop.
    const auto current_n_chosen_assignments =
        enriched_node.node_wrapper.get().chosen_assignments.size();

    {
      const auto current_number_of_assignments = m_assignments.size();

      if (!shared_data.close_vertices_filter.reduce(
              shared_data.fixed_data, m_assignments,
              assignments_processed_in_this_node, enriched_node.node_wrapper)) {
        return ReductionResult::FAILURE;
      }
      if (current_number_of_assignments != m_assignments.size()) {
        // Propagate the new assignments and edge weights immediately.
        TKET_ASSERT(current_number_of_assignments < m_assignments.size());
        continue;
      }
      TKET_ASSERT(
          current_n_chosen_assignments ==
          enriched_node.node_wrapper.get().chosen_assignments.size());
    }
    {
      const auto current_number_of_assignments = m_assignments.size();
      if (!shared_data.fixed_data.hall_set_reducer.reduce(
              enriched_node.node_wrapper, *this)) {
        return ReductionResult::FAILURE;
      }
      if (current_number_of_assignments != m_assignments.size()) {
        TKET_ASSERT(current_number_of_assignments < m_assignments.size());
        continue;
      }
      TKET_ASSERT(
          current_n_chosen_assignments ==
          enriched_node.node_wrapper.get().chosen_assignments.size());
    }
    /*
    // TODO: a very annoying intermittent Heisenbug appears to be
    // in here somewhere...but no time to fix it now, will come back to it...
    {
      const auto current_number_of_assignments = m_assignments.size();
      if (!shared_data.derived_graphs_reducer.reduce_domains(
              shared_data.fixed_data, m_assignments,
              assignments_processed_in_this_node, enriched_node.node_wrapper,
              shared_data.derived_graphs_filter.get_derived_pattern_graphs(),
              shared_data.derived_graphs_filter.get_derived_target_graphs())) {
        return ReductionResult::FAILURE;
      }
      if (current_number_of_assignments != m_assignments.size()) {
        TKET_ASSERT(current_number_of_assignments < m_assignments.size());
        continue;
      }
      TKET_ASSERT(
          current_n_chosen_assignments ==
          enriched_node.node_wrapper.get().chosen_assignments.size());
    }
    //*/

    // Nogood detectors; these don't REDUCE anything, they just try to detect
    // if we're already at a dead end (we just don't know it yet...)

    const auto& node = enriched_node.node_wrapper.get();
    if (!shared_data.fixed_data.problem_is_unweighted &&
        m_weight_nogood_detector_manager.should_activate_detector(
            node.current_scalar_product, max_weight, node.total_p_edge_weights,
            shared_data.fixed_data.total_p_edge_weights, m_assignments.size(),
            node.pattern_v_to_possible_target_v.size())) {
      const WeightWSM max_extra_weight =
          max_weight - node.current_scalar_product;
      const auto weight_nogood_result =
          shared_data.fixed_data.weight_nogood_detector(
              shared_data.fixed_data, node.pattern_v_to_possible_target_v,
              m_assignments, max_extra_weight);

      if (weight_nogood_result.assignment_with_invalid_t_vertex) {
        const auto& p_vertices_map =
            shared_data.fixed_data.pattern_neighbours_data.get_map();
        const auto& t_vertex =
            weight_nogood_result.assignment_with_invalid_t_vertex->second;
        for (const auto& entry : p_vertices_map) {
          const auto& p_vertex = entry.first;
          erase_assignment(p_vertex, t_vertex);
        }
        return ReductionResult::FAILURE;
      }

      if (!weight_nogood_result.extra_weight_lower_bound) {
        m_weight_nogood_detector_manager.register_success();
        return ReductionResult::FAILURE;
      }
      const auto extra_weight_lower_bound =
          weight_nogood_result.extra_weight_lower_bound.value();
      TKET_ASSERT(extra_weight_lower_bound <= max_extra_weight);

      m_weight_nogood_detector_manager.register_lower_bound_failure(
          node.current_scalar_product, max_weight, extra_weight_lower_bound);
    }

    // If we've reached here, then all the reducers/filters reduced
    // the domain sizes and searched for inconsistencies, but didn't find any;
    // thus this node is now fully reduced, ready for more searching.
    return ReductionResult::SUCCESS;
  }

  // If we reach here, we've broken out of the loop;
  // this can only happen if all vertices have been assigned,
  // AND they are all valid (all edges are assigned also);
  // this has been checked.
  TKET_ASSERT(
      m_assignments.size() == shared_data.fixed_data.pattern_neighbours_data
                                  .get_number_of_nonisolated_vertices());

  return ReductionResult::FINISHED;
}

bool SearchBranch::backtrack() {
  if (m_level == 0) {
    return false;
  }
  for (const auto& entry :
       m_enriched_nodes[m_level].node_wrapper.get().chosen_assignments) {
    if (m_assignments.count(entry.first) != 0) {
      TKET_ASSERT(m_assignments.at(entry.first) == entry.second);
    }
    m_assignments.erase(entry.first);
  }
  --m_level;
  return true;
}

const SearchBranch::EnrichedNodes& SearchBranch::get_data(
    std::size_t& level) const {
  level = m_level;
  return m_enriched_nodes;
}

const SearchNodeWrapper& SearchBranch::get_current_node_wrapper() const {
  return m_enriched_nodes.at(m_level).node_wrapper;
}

bool SearchBranch::move_down_has_been_called() const {
  return m_move_down_has_been_called;
}

void SearchBranch::move_down(VertexWSM p_vertex, VertexWSM t_vertex) {
  m_move_down_has_been_called = true;
  {
    auto& node = m_enriched_nodes.at(m_level).node_wrapper.get_mutable();
    auto& domain_data = node.pattern_v_to_possible_target_v.at(p_vertex);
    TKET_ASSERT(domain_data.erase(t_vertex) == 1);
    // Now we can move down; our current choice has been erased,
    // so will not be repeated in future.
  }

  ++m_level;
  if (m_level >= m_enriched_nodes.size()) {
    m_enriched_nodes.resize(m_level + 1);
  }
  m_enriched_nodes[m_level]
      .node_wrapper.get_mutable()
      .initialise_from_assignment(
          p_vertex, t_vertex, m_enriched_nodes[m_level - 1].node_wrapper.get());

  // We're maybe NOT fully reduced, but that's OK.
  m_enriched_nodes[m_level].clear_enriched_data();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
