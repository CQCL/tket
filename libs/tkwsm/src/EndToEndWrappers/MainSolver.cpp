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

#include "tkwsm/EndToEndWrappers/MainSolver.hpp"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/EndToEndWrappers/PreSearchComponents.hpp"
#include "tkwsm/EndToEndWrappers/SearchComponents.hpp"
#include "tkwsm/GraphTheoretic/DomainInitialiser.hpp"
#include "tkwsm/WeightPruning/WeightChecker.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

typedef std::chrono::steady_clock Clock;

MainSolver::MainSolver(
    const GraphEdgeWeights& pattern_edges, const GraphEdgeWeights& target_edges,
    const MainSolverParameters& parameters)
    : m_pattern_vertex_relabelling(pattern_edges),
      m_target_vertex_relabelling(target_edges),
      m_pattern_neighbours_data(
          m_pattern_vertex_relabelling.new_edges_and_weights),
      m_target_neighbours_data(
          m_target_vertex_relabelling.new_edges_and_weights) {
  const auto num_p_vertices =
      m_pattern_neighbours_data.get_number_of_nonisolated_vertices();
  if (num_p_vertices == 0) {
    m_solution_data.trivial_weight_lower_bound = 0;
    m_solution_data.trivial_weight_initial_upper_bound = 0;
    m_solution_data.finished = true;
    return;
  }
  const auto num_t_vertices =
      m_target_neighbours_data.get_number_of_nonisolated_vertices();
  {
    const auto number_of_possible_t_edges =
        (num_t_vertices * (num_t_vertices - 1)) / 2;
    m_solution_data.target_is_complete =
        number_of_possible_t_edges ==
        m_target_neighbours_data.get_number_of_edges();
  }

  // Start off assuming that it's impossible. So L = +inf, U = 0 make
  // sense mathematically (infimum over empty set is +infinity, etc. etc.)
  m_solution_data.trivial_weight_initial_upper_bound = 0;
  set_maximum(m_solution_data.trivial_weight_lower_bound);

  if (m_pattern_neighbours_data.get_number_of_edges() >
          m_target_neighbours_data.get_number_of_edges() ||
      num_p_vertices > num_t_vertices) {
    m_solution_data.finished = true;
    return;
  }

  const auto init_start = Clock::now();

  m_pre_search_components_ptr = std::make_unique<PreSearchComponents>(
      m_pattern_neighbours_data, m_target_neighbours_data);
  TKET_ASSERT(m_pre_search_components_ptr);

  {
    DomainInitialiser::InitialDomains initial_domains;
    const bool initialisation_succeeded =
        DomainInitialiser::full_initialisation(
            initial_domains, m_pattern_neighbours_data,
            m_pre_search_components_ptr->pattern_near_ndata,
            m_target_neighbours_data,
            m_pre_search_components_ptr->target_near_ndata,
            parameters.max_distance_for_domain_initialisation_distance_filter);

    if (initialisation_succeeded) {
      m_search_components_ptr = std::make_unique<SearchComponents>();
      TKET_ASSERT(m_search_components_ptr);

      m_search_branch_ptr = std::make_unique<SearchBranch>(
          initial_domains, m_pattern_neighbours_data,
          m_pre_search_components_ptr->pattern_near_ndata,
          m_target_neighbours_data,
          m_pre_search_components_ptr->target_near_ndata,
          parameters.max_distance_for_distance_reduction_during_search,
          m_solution_data.extra_statistics);
    }
    m_solution_data.initialisation_time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - init_start)
            .count();

    if (!initialisation_succeeded) {
      m_solution_data.finished = true;
      return;
    }
  }

  // The problem is not obviously impossible, so fill in more information.
  {
    std::vector<WeightWSM> p_weights =
        m_pattern_neighbours_data.get_weights_expensive();
    std::sort(p_weights.begin(), p_weights.end());

    std::vector<WeightWSM> t_weights =
        m_target_neighbours_data.get_weights_expensive();
    std::sort(t_weights.begin(), t_weights.end());

    TKET_ASSERT(
        p_weights.size() == m_pattern_neighbours_data.get_number_of_edges());
    TKET_ASSERT(
        t_weights.size() == m_target_neighbours_data.get_number_of_edges());
    TKET_ASSERT(p_weights.size() <= t_weights.size());
    m_solution_data.total_p_edge_weights =
        std::accumulate(p_weights.cbegin(), p_weights.cend(), WeightWSM(0));

    // Now get trivial lower, upper bounds.
    // When considering   sum a[i].b[p(i)],  where p can be any permutation
    // and a,b are sequences of real numbers,
    // the min value arises when p is such that a,b are in opposite order.
    // The max value arises when a,b are in the same order.
    m_solution_data.trivial_weight_lower_bound = 0;

    // Lower bound: p-weights increasing, t-weights decreasing;
    // BUT we only need use the smallest t-weights, for the minimum.
    for (unsigned ii = 0; ii < p_weights.size(); ++ii) {
      const auto product = get_product_or_throw(
          p_weights[ii], t_weights[(p_weights.size() - 1) - ii]);
      m_solution_data.trivial_weight_lower_bound =
          get_sum_or_throw(m_solution_data.trivial_weight_lower_bound, product);
    }

    // Now, for the MAXIMUM, we use the LARGEST t-weights,
    // as well as summing both vectors in increasing order.
    const unsigned offset = t_weights.size() - p_weights.size();
    m_solution_data.trivial_weight_initial_upper_bound = 0;

    for (unsigned ii = 0; ii < p_weights.size(); ++ii) {
      const auto product =
          get_product_or_throw(p_weights[ii], t_weights[ii + offset]);
      m_solution_data.trivial_weight_initial_upper_bound = get_sum_or_throw(
          m_solution_data.trivial_weight_initial_upper_bound, product);
    }
  }

  if (m_solution_data.trivial_weight_lower_bound !=
      m_solution_data.trivial_weight_initial_upper_bound) {
    // It's not an unweighted problem, so it's worth checking for weights.
    m_search_branch_ptr->activate_weight_checker(
        m_solution_data.total_p_edge_weights);
  }

  if (m_solution_data.initialisation_time_ms >= parameters.timeout_ms) {
    return;
  }

  const auto search_start_time = Clock::now();
  const auto desired_search_end_time =
      search_start_time + std::chrono::milliseconds(parameters.timeout_ms);

  if (parameters.iterations_timeout != 0) {
    internal_solve(
        parameters, parameters.iterations_timeout, desired_search_end_time);
  }

  m_solution_data.search_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          Clock::now() - search_start_time)
          .count();
}

MainSolver::~MainSolver() {}

const SolutionData& MainSolver::get_solution_data() const {
  if (m_search_branch_ptr) {
    m_search_branch_ptr->get_updated_extra_statistics();
  }
  if (m_pattern_vertex_relabelling.new_to_old_vertex_labels.empty() &&
      m_target_vertex_relabelling.new_to_old_vertex_labels.empty()) {
    return m_solution_data;
  }
  // We must relabel some vertices.
  m_solution_data_original_vertices = m_solution_data;
  for (SolutionWSM& solution : m_solution_data_original_vertices.solutions) {
    if (!m_pattern_vertex_relabelling.new_to_old_vertex_labels.empty()) {
      for (std::pair<VertexWSM, VertexWSM>& assignment : solution.assignments) {
        assignment.first =
            m_pattern_vertex_relabelling.new_to_old_vertex_labels.at(
                assignment.first);
      }
    }
    if (!m_target_vertex_relabelling.new_to_old_vertex_labels.empty()) {
      for (std::pair<VertexWSM, VertexWSM>& assignment : solution.assignments) {
        assignment.second =
            m_target_vertex_relabelling.new_to_old_vertex_labels.at(
                assignment.second);
      }
    }
  }
  return m_solution_data_original_vertices;
}

void MainSolver::solve(const MainSolverParameters& parameters) {
  if (m_solution_data.finished) {
    return;
  }
  const auto search_start_time = Clock::now();
  const auto desired_search_end_time =
      search_start_time + std::chrono::milliseconds(parameters.timeout_ms);
  const auto max_iterations_opt = get_checked_sum(
      m_solution_data.iterations, parameters.iterations_timeout);
  decltype(m_solution_data.iterations) max_iterations;
  if (max_iterations_opt) {
    max_iterations = max_iterations_opt.value();
  } else {
    set_maximum(max_iterations);
  }
  internal_solve(parameters, max_iterations, desired_search_end_time);
  m_solution_data.search_time_ms +=
      std::chrono::duration_cast<std::chrono::milliseconds>(
          Clock::now() - search_start_time)
          .count();
}

static bool terminate_with_enough_full_solutions(
    const MainSolverParameters& parameters, const SolutionData& solution_data) {
  if (parameters.for_multiple_full_solutions_the_max_number_to_obtain == 0) {
    return parameters.terminate_with_first_full_solution &&
           !solution_data.solutions.empty();
  }
  return solution_data.solutions.size() >=
         parameters.for_multiple_full_solutions_the_max_number_to_obtain;
}

void MainSolver::internal_solve(
    const MainSolverParameters& parameters, std::size_t max_iterations,
    const std::chrono::steady_clock::time_point& desired_end_time) {
  if (m_solution_data.finished ||
      terminate_with_enough_full_solutions(parameters, m_solution_data)) {
    return;
  }
  TKET_ASSERT(m_pre_search_components_ptr);
  TKET_ASSERT(m_search_branch_ptr);

  SearchBranch::ReductionParameters reduction_parameters;
  decltype(reduction_parameters.max_weight) initial_weight_upper_bound;

  if (parameters.weight_upper_bound_constraint) {
    initial_weight_upper_bound =
        parameters.weight_upper_bound_constraint.value();
  } else {
    set_maximum(initial_weight_upper_bound);
  }

  while (m_solution_data.iterations < max_iterations) {
    // Set the maximum weight.
    if (m_solution_data.solutions.empty()) {
      // We have no solution yet.
      reduction_parameters.max_weight = initial_weight_upper_bound;
    } else {
      if (parameters.for_multiple_full_solutions_the_max_number_to_obtain > 0) {
        // Because we want to store multiple solutions,
        // we allow equally good, or worse, solutions.
        reduction_parameters.max_weight = initial_weight_upper_bound;
      } else {
        // We only want a single solution, so only allow a strictly better
        // solution. Force exactly one solution to be stored - the best. We
        // should only have one, but previously we've stored more. So just
        // choose the best.
        {
          auto best_scalar_product =
              m_solution_data.solutions[0].scalar_product;
          unsigned best_index = 0;
          for (unsigned index = 1; index < m_solution_data.solutions.size();
               ++index) {
            if (m_solution_data.solutions[index].scalar_product <
                best_scalar_product) {
              best_scalar_product =
                  m_solution_data.solutions[index].scalar_product;
              best_index = index;
            }
          }
          if (best_index > 0) {
            m_solution_data.solutions[0] =
                m_solution_data.solutions[best_index];
          }
        }
        m_solution_data.solutions.resize(1);
        reduction_parameters.max_weight =
            m_solution_data.solutions[0].scalar_product;
        if (reduction_parameters.max_weight == 0) {
          // We can't do better than zero!
          m_solution_data.finished = true;
          return;
        }
        // Make it strictly better.
        --reduction_parameters.max_weight;
        reduction_parameters.max_weight = std::min(
            reduction_parameters.max_weight, initial_weight_upper_bound);
      }
    }

    if (reduction_parameters.max_weight <
        m_solution_data.trivial_weight_lower_bound) {
      m_solution_data.finished = true;
      return;
    }

    // Now we can search!
    ++m_solution_data.iterations;

    // On the first move ONLY, we don't backtrack; but we also haven't reduced.
    if (m_solution_data.iterations == 1) {
      if (!m_search_branch_ptr->reduce_current_node(reduction_parameters)) {
        m_solution_data.finished = true;
        return;
      }
    } else {
      if (!m_search_branch_ptr->backtrack(reduction_parameters)) {
        m_solution_data.finished = true;
        return;
      }
    }
    if (move_down_from_reduced_node(reduction_parameters)) {
      // We've GOT a complete solution! Note that it MUST be good enough
      // to add, since we've set the max weight already.
      // We also already checked that we haven't yet got too many,
      // if we're storing more than one.
      add_solution_from_final_node(parameters, reduction_parameters);
      if (terminate_with_enough_full_solutions(parameters, m_solution_data)) {
        return;
      }
    }
    if (Clock::now() >= desired_end_time) {
      return;
    }
  }
}

void MainSolver::add_solution_from_final_node(
    const MainSolverParameters& parameters,
    const SearchBranch::ReductionParameters& reduction_parameters) {
  TKET_ASSERT(m_pre_search_components_ptr);
  TKET_ASSERT(m_search_branch_ptr);

  const DomainsAccessor& accessor = m_search_branch_ptr->get_domains_accessor();
  const WeightWSM scalar_product = accessor.get_scalar_product();

  TKET_ASSERT(
      accessor.get_total_p_edge_weights() ==
      m_solution_data.total_p_edge_weights);

  TKET_ASSERT(scalar_product <= reduction_parameters.max_weight);

  // We'll overwrite the solution into back().
  if (parameters.for_multiple_full_solutions_the_max_number_to_obtain > 0 ||
      m_solution_data.solutions.empty()) {
    m_solution_data.solutions.emplace_back();
  }
  std::vector<std::pair<VertexWSM, VertexWSM>>& assignments =
      m_solution_data.solutions.back().assignments;
  const auto number_of_pv = accessor.get_number_of_pattern_vertices();
  assignments.clear();
  assignments.reserve(number_of_pv);

  for (unsigned pv = 0; pv < number_of_pv; ++pv) {
    const BitsetInformation bitset_information(accessor.get_domain(pv));
    TKET_ASSERT(bitset_information.single_element);
    assignments.emplace_back(pv, bitset_information.single_element.value());
  }
  m_solution_data.solutions.back().scalar_product = scalar_product;

  m_solution_data.solutions.back().total_p_edges_weight =
      m_solution_data.total_p_edge_weights;
}

bool MainSolver::move_down_from_reduced_node(
    const SearchBranch::ReductionParameters& reduction_parameters) {
  TKET_ASSERT(m_search_components_ptr);
  TKET_ASSERT(m_search_branch_ptr);

  for (;;) {
    const VariableOrdering::Result next_var_result =
        m_search_components_ptr->variable_ordering.get_variable(
            m_search_branch_ptr->get_domains_accessor_nonconst(),
            m_search_components_ptr->rng);

    if (next_var_result.empty_domain) {
      return false;
    }
    if (!next_var_result.variable_opt) {
      // If no variable to choose, it means we've got a full solution
      // (every PV is assigned).
      // Before we reached here, we already checked that any new pattern edges
      // joining any newly assigned PV to an existing assigned PV
      // ARE indeed mapped to valid target edges,
      // so we don't need any further validity check.
      break;
    }
    // We've chosen a variable (i.e., PV) to assign:
    const VertexWSM& next_pv = next_var_result.variable_opt.value();

    // Now choose a value (i.e., some TV in Domain(PV)).
    // Thus the new assignment will be next_pv -> next_tv.
    const VertexWSM next_tv =
        m_search_components_ptr->value_ordering.get_target_value(
            m_search_branch_ptr->get_domains_accessor().get_domain(next_pv),
            m_target_neighbours_data, m_search_components_ptr->rng);

    m_search_branch_ptr->move_down(next_pv, next_tv);
    if (!m_search_branch_ptr->reduce_current_node(reduction_parameters)) {
      return false;
    }
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
