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

#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

#include <stdexcept>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SearchNodeWrapper::SearchNodeWrapper() {}

SearchNodeWrapper::SearchNodeWrapper(SearchNode node)
    : m_node(std::move(node)) {}

const SearchNode& SearchNodeWrapper::get() const { return m_node; }

SearchNode& SearchNodeWrapper::get_mutable() { return m_node; }

SearchNodeWrapper& SearchNodeWrapper::add_p_edge_weights(WeightWSM dw) {
  m_node.total_p_edge_weights += dw;
  return *this;
}

SearchNodeWrapper& SearchNodeWrapper::add_scalar_product(WeightWSM dw) {
  m_node.current_scalar_product += dw;
  return *this;
}

std::size_t SearchNodeWrapper::remove_element_from_domain(
    VertexWSM pv, VertexWSM target_vertex, Assignments& assignments) {
  const auto citer = m_node.pattern_v_to_possible_target_v.find(pv);
  if (citer == m_node.pattern_v_to_possible_target_v.cend() ||
      citer->second.empty()) {
    return 0;
  }
  return remove_element_from_domain(
      pv, target_vertex, citer->second, assignments);
}

std::size_t SearchNodeWrapper::remove_element_from_domain(
    VertexWSM pv, VertexWSM target_vertex, std::set<VertexWSM>& domain,
    Assignments& assignments) {
  if (domain.erase(target_vertex) == 1) {
    if (domain.empty()) {
      return 0;
    }
    if (domain.size() == 1) {
      const auto single_tv_opt = get_optional_value(assignments, pv);
      const VertexWSM& this_tv = *domain.cbegin();
      if (single_tv_opt) {
        if (this_tv != single_tv_opt.value()) {
          return 0;
        }
      } else {
        // Definitely a new assignment.
        m_node.chosen_assignments.emplace_back(pv, this_tv);
        assignments[pv] = this_tv;
      }
    }
  }
  return domain.size();
}

void SearchNodeWrapper::overwrite_domain(
    const std::vector<VertexWSM>& new_domain, VertexWSM pv,
    Assignments& assignments) {
  if (new_domain.size() <= 1) {
    m_node.pattern_v_to_possible_target_v.erase(pv);
    if (new_domain.empty()) {
      throw std::runtime_error(
          "SearchNodeWrapper::overwrite_domain : empty new domain");
    }
    const auto new_tv = *new_domain.cbegin();
    const auto existing_tv_opt = get_optional_value(assignments, pv);
    if (existing_tv_opt) {
      if (existing_tv_opt.value() != new_tv) {
        throw std::runtime_error(
            "SearchNodeWrapper::overwrite_domain : TV mismatch");
      }
    } else {
      assignments[pv] = new_tv;
      m_node.chosen_assignments.emplace_back(pv, new_tv);
    }
    return;
  }
  auto& domain = m_node.pattern_v_to_possible_target_v.at(pv);
  domain = {new_domain.cbegin(), new_domain.cend()};
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
