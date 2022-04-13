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

#include "WeightSubgrMono/Searching/NodeWSM.hpp"

#include <sstream>

#include "Utils/Assert.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NodeWSM::NodeWSM() : m_scalar_product(0), m_total_p_edge_weights(0) {}

const std::vector<std::pair<VertexWSM, VertexWSM>>&
NodeWSM::get_new_assignments() const {
  return m_new_assignments;
}

void NodeWSM::set_scalar_product(WeightWSM scalar_product) {
  m_scalar_product = scalar_product;
}
WeightWSM NodeWSM::get_scalar_product() const { return m_scalar_product; }
void NodeWSM::set_total_pattern_edge_weights(WeightWSM new_weight) {
  m_total_p_edge_weights = new_weight;
}
WeightWSM NodeWSM::get_total_pattern_edge_weights() const {
  return m_total_p_edge_weights;
}
void NodeWSM::clear_new_assignments() { m_new_assignments.clear(); }

void NodeWSM::overwrite_domain(
    VertexWSM pattern_v, std::set<VertexWSM> new_domain) {
  TKET_ASSERT(!new_domain.empty());
  auto& domain = m_pattern_v_to_possible_target_v.at(pattern_v);
  if (domain.size() == 1) {
    // We can't "undo" an assignment.
    TKET_ASSERT(new_domain.size() == 1);
    TKET_ASSERT(*domain.cbegin() == *new_domain.cbegin());
    return;
  }
  domain = std::move(new_domain);
  if (domain.size() == 1) {
    m_new_assignments.emplace_back(pattern_v, *domain.cbegin());
  }
}

void NodeWSM::overwrite_domain(
    VertexWSM pattern_v, const std::vector<VertexWSM>& new_domain) {
  TKET_ASSERT(!new_domain.empty());
  auto& domain = m_pattern_v_to_possible_target_v.at(pattern_v);
  std::optional<VertexWSM> existing_tv_opt;
  if (domain.size() == 1) {
    // We can't "undo" an assignment.
    existing_tv_opt = *domain.cbegin();
  }
  domain = {new_domain.cbegin(), new_domain.cend()};
  if (domain.size() != 1) {
    return;
  }
  if (existing_tv_opt) {
    TKET_ASSERT(existing_tv_opt.value() == *domain.cbegin());
  } else {
    m_new_assignments.emplace_back(pattern_v, *domain.cbegin());
  }
}

void NodeWSM::force_assignment(
    const std::pair<VertexWSM, VertexWSM>& assignment) {
  auto& domain = m_pattern_v_to_possible_target_v.at(assignment.first);
  if (domain.size() == 1 && *domain.cbegin() == assignment.second) {
    return;
  }
  domain.clear();
  domain.insert(assignment.second);
  m_new_assignments.emplace_back(assignment);
}

NodeWSM::ErasureResult NodeWSM::erase_assignment(
    const std::pair<VertexWSM, VertexWSM>& assignment) {
  auto& domain = m_pattern_v_to_possible_target_v.at(assignment.first);
  auto iter = domain.find(assignment.second);

  ErasureResult result;
  if (iter == domain.end()) {
    result.assignment_was_possible = false;
  } else {
    result.assignment_was_possible = true;
    domain.erase(iter);
    if (domain.size() == 1) {
      m_new_assignments.emplace_back(assignment.first, *domain.cbegin());
    }
  }
  result.valid = !domain.empty();
  return result;
}

NodeWSM::SoftErasureResult NodeWSM::attempt_to_erase_assignment(
    const std::pair<VertexWSM, VertexWSM>& assignment) {
  auto& domain = m_pattern_v_to_possible_target_v.at(assignment.first);
  auto iter = domain.find(assignment.second);
  if (iter == domain.end()) {
    return SoftErasureResult::TV_WAS_NOT_PRESENT;
  }
  if (domain.size() <= 2) {
    return SoftErasureResult::TV_REMAINS;
  }
  domain.erase(iter);
  return SoftErasureResult::TV_REMAINS;
}

bool NodeWSM::alldiff_reduce(std::size_t n_assignments_processed) {
  for (auto index = n_assignments_processed; index < m_new_assignments.size();
       ++index) {
    const auto& assignment = m_new_assignments[index];
    const VertexWSM& pv = assignment.first;
    const VertexWSM& tv = assignment.second;

    for (auto& entry : m_pattern_v_to_possible_target_v) {
      auto& domain = entry.second;

      if (entry.first == pv) {
        TKET_ASSERT(domain.size() == 1);
        TKET_ASSERT(*domain.cbegin() == tv);
        continue;
      }
      if (domain.erase(tv) == 1) {
        switch (domain.size()) {
          case 0:
            return false;
          case 1:
            m_new_assignments.emplace_back(entry.first, *domain.cbegin());
            break;
          default:
            break;
        }
      }
    }
  }
  return true;
}

const PossibleAssignments& NodeWSM::get_possible_assignments() const {
  return m_pattern_v_to_possible_target_v;
}

void NodeWSM::set_possible_assignments(
    PossibleAssignments possible_assignments) {
  m_pattern_v_to_possible_target_v = std::move(possible_assignments);
  m_new_assignments.clear();
  for (const auto& entry : m_pattern_v_to_possible_target_v) {
    const auto& domain = entry.second;
    if (domain.size() == 1) {
      m_new_assignments.emplace_back(entry.first, *domain.cbegin());
    }
  }
}

std::string NodeWSM::str() const {
  std::stringstream ss;
  ss << "\n"
     << m_pattern_v_to_possible_target_v.size() << " p-vertices. Domains: ";
  for (const auto& entry : m_pattern_v_to_possible_target_v) {
    ss << "\nDom(" << entry.first << ") = {";
    for (const auto& tv : entry.second) {
      ss << " " << tv << " ";
    }
    ss << "}";
  }
  ss << "\nAssigned p-edges weights " << m_total_p_edge_weights << "; sc.prod "
     << m_scalar_product << "; new assignments:\n[";
  for (const auto& entry : m_new_assignments) {
    ss << " " << entry.first << ":" << entry.second << " ";
  }
  ss << "]\n";
  return ss.str();
}

/*
bool NodeWSM::create_initial_node(PossibleAssignments
pattern_v_to_possible_target_v) {
  //m_superficially_valid = true;
  m_number_of_assignments_processed_by_alldiff_propagator = 0;
  m_scalar_product = 0;
  m_total_p_edge_weights = 0;
  m_pattern_v_to_possible_target_v = std::move(pattern_v_to_possible_target_v);
  m_initial_assigned_p_vertex = std::nullopt;
  m_newly_assigned_p_vertices.clear();
  for(const auto& entry : m_pattern_v_to_possible_target_v) {
    switch(entry.second.size()) {
      case 0: return false;
      case 1: m_newly_assigned_p_vertices.emplace_back(entry.first);
      default: break;
    }
  }
  return reduce_with_alldiff_propagation();
}
*/

/*
NodeWSM::InitialisationResult NodeWSM::initialise_from_assignment(
      VertexWSM pattern_v, VertexWSM target_v, NodeWSM& previous_node) {

  InitialisationResult result;

  // Partially reduce the previous node.
  {
    auto iter = previous_node.m_pattern_v_to_possible_target_v.find(pattern_v);
    TKET_ASSERT(iter != previous_node.m_pattern_v_to_possible_target_v.end());
    auto& domain = iter->second;
    TKET_ASSERT(domain.size() >= 2);
    TKET_ASSERT(domain.erase(target_v) == 1);
    switch (domain.size()) {
    case 0: result.previous_node_is_valid = false; break;
    case 1:
previous_node.m_newly_assigned_p_vertices.push_back(*domain.cbegin());
      // fall through
    default:
      result.previous_node_is_valid = true;
      break;
    }
  }
  // Copy across the data selectively.
  m_number_of_assignments_processed_by_alldiff_propagator = 0;
  m_scalar_product = previous_node.m_scalar_product;
  m_total_p_edge_weights = previous_node.m_total_p_edge_weights;
  m_initial_assigned_p_vertex = pattern_v;

  // Probably slightly quicker than a loop over all keys except pv.
  m_pattern_v_to_possible_target_v =
previous_node.m_pattern_v_to_possible_target_v; auto& domain =
m_pattern_v_to_possible_target_v.at(pattern_v); domain.clear();
  domain.insert(target_v);
  m_newly_assigned_p_vertices.resize(1);
  m_newly_assigned_p_vertices[0] = pattern_v;
  result.valid = reduce_with_alldiff_propagation();
  return result;
}
*/

/*
bool NodeWSM::reduce_with_alldiff_propagation() {
  TKET_ASSERT(m_number_of_assignments_processed_by_alldiff_propagator <=
m_newly_assigned_p_vertices.size());

  while(m_number_of_assignments_processed_by_alldiff_propagator <
m_newly_assigned_p_vertices.size()) { const VertexWSM pv =
m_newly_assigned_p_vertices[m_number_of_assignments_processed_by_alldiff_propagator];
    ++m_number_of_assignments_processed_by_alldiff_propagator;
    const auto& domain = m_pattern_v_to_possible_target_v.at(pv);
    TKET_ASSERT(domain.size() == 1);
    const VertexWSM tv = *domain.cbegin();

    // Now erase tv from all other domains.
    for(auto& entry : m_pattern_v_to_possible_target_v) {
      if(entry.first == pv) {
        continue;
      }
      if(entry.second.erase(tv) == 1) {
        // It's changed.
        switch(entry.second.size()) {
          case 0: return false;
          case 1: m_newly_assigned_p_vertices.push_back(pv);
          default: break;
        }
      }
    }
  }
  return true;
}
*/

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
