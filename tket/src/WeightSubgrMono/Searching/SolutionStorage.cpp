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

#include "WeightSubgrMono/Searching/SolutionStorage.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

#include "Utils/Assert.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static void print_assignments(const SolutionWSM& solution) {
  std::cerr << "\nAssigned " << solution.assignments.size() << " vars: ";
  for (const auto& entry : solution.assignments) {
    std::cerr << " " << entry.first << ":" << entry.second << " ";
  }
  std::cerr << "\n";
}

static void copy_solution_into_storage(
    const SearchBranch& branch, WeightWSM new_p_edges_weight,
    WeightWSM new_scalar_prod, SolutionWSM& solution) {
  solution.assignments.clear();

  for (const auto& entry : branch.get_assignments()) {
    solution.assignments.emplace_back(entry);
  }
  solution.total_p_edges_weight = new_p_edges_weight;
  solution.total_scalar_product_weight = new_scalar_prod;
}

SolutionStorage::SolutionStorage() : m_log_level(0), m_pruning_weight(0) {}

const SolutionWSM& SolutionStorage::best_solution() const { return m_solution; }

SolutionStorage& SolutionStorage::set_pruning_weight(WeightWSM weight) {
  TKET_ASSERT(weight > 0);
  m_pruning_weight = weight;
  return *this;
}

SolutionStorage& SolutionStorage::set_log_level(unsigned log_level) {
  m_log_level = log_level;
  return *this;
}

std::optional<WeightWSM> SolutionStorage::get_acceptable_scalar_product()
    const {
  if (m_pruning_weight == 0) {
    if (m_solution.complete) {
      return m_solution.total_scalar_product_weight;
    }
    return {};
  }
  auto weight = m_pruning_weight;
  if (m_solution.complete) {
    if (m_solution.total_scalar_product_weight <= 1) {
      return 0;
    }
    WeightWSM smaller_weight = m_solution.total_scalar_product_weight;
    --smaller_weight;
    weight = std::min(weight, smaller_weight);
  }
  return weight;
}

bool SolutionStorage::add_full_solution(const SearchBranch& branch) {
  const auto& node = branch.get_current_node_wrapper().get();
  // Everything should be assigned, if it really is a full solution.
  TKET_ASSERT(node.pattern_v_to_possible_target_v.empty());

  const auto new_p_edges_weight = node.total_p_edge_weights;
  const auto new_scalar_prod = node.current_scalar_product;

  TKET_ASSERT(new_p_edges_weight >= m_solution.total_p_edges_weight);

  if (m_solution.complete) {
    TKET_ASSERT(new_p_edges_weight == m_solution.total_p_edges_weight);
  }
  const auto weight_optional = get_acceptable_scalar_product();
  if (weight_optional && new_scalar_prod > weight_optional.value()) {
    return false;
  }

  const auto current_number_of_assignments = m_solution.assignments.size();
  copy_solution_into_storage(
      branch, new_p_edges_weight, new_scalar_prod, m_solution);
  TKET_ASSERT(current_number_of_assignments <= m_solution.assignments.size());
  if (m_solution.complete) {
    TKET_ASSERT(current_number_of_assignments == m_solution.assignments.size());
  }
  m_solution.complete = true;
  if (m_log_level > 0) {
    std::cerr << "\n#### NEW FULL soln: sc.prod "
              << m_solution.total_scalar_product_weight << "; p-edges "
              << m_solution.total_p_edges_weight;
    if (m_log_level > 1) {
      print_assignments(m_solution);
    }
  }
  return true;
}

bool SolutionStorage::add_partial_solution(const SearchBranch& branch) {
  const auto& node = branch.get_current_node_wrapper().get();
  const auto& new_p_edges_weight = node.total_p_edge_weights;
  if (m_solution.complete) {
    TKET_ASSERT(m_solution.total_p_edges_weight >= new_p_edges_weight);
    return false;
  }
  const auto& new_scalar_prod = node.current_scalar_product;
  const auto weight_optional = get_acceptable_scalar_product();
  if (weight_optional && new_scalar_prod > weight_optional.value()) {
    return false;
  }
  // EITHER we're better than the required weight,
  // or there IS no required weight.
  // But, are we ALSO better than the current partial solution?
  // Either we embed more p-edges, OR an equal amount,
  // but with smaller scalar product.
  if (new_p_edges_weight <= m_solution.total_p_edges_weight ||
      (new_p_edges_weight == m_solution.total_p_edges_weight &&
       new_scalar_prod >= m_solution.total_scalar_product_weight)) {
    return false;
  }
  copy_solution_into_storage(
      branch, new_p_edges_weight, new_scalar_prod, m_solution);
  if (m_log_level > 0) {
    std::cerr << "\n\n#### NEW part soln: sc.prod "
              << m_solution.total_scalar_product_weight << "; p-edges "
              << m_solution.total_p_edges_weight;
    if (m_log_level > 1) {
      print_assignments(m_solution);
    }
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
