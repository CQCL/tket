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

#include "WeightSubgrMono/Searching/DomainsAccessor.hpp"

#include <sstream>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/Searching/NodesRawData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DomainsAccessor::DomainsAccessor(NodesRawDataWrapper& raw_data_wrapper)
    : m_raw_data(raw_data_wrapper.m_raw_data) {}

const std::vector<VertexWSM>& DomainsAccessor::get_pattern_vertices() const {
  return m_raw_data.pattern_vertices;
}

const std::set<VertexWSM>& DomainsAccessor::get_domain(VertexWSM pv) const {
  const auto& data = m_raw_data.domains_data.at(pv);
  return data.entries[data.entries_back_index].domain;
}

const std::set<VertexWSM>&
DomainsAccessor::get_unassigned_pattern_vertices_superset() const {
  for (unsigned index = m_raw_data.current_node_level; index != 0; --index) {
    if (m_raw_data.nodes_data[index].nogood) {
      continue;
    }
    const auto& vertices =
        m_raw_data.nodes_data[index].unassigned_vertices_superset;
    if (!vertices.empty()) {
      return vertices;
    }
  }
  return m_raw_data.nodes_data[0].unassigned_vertices_superset;
}

bool DomainsAccessor::domain_created_in_current_node(VertexWSM pv) const {
  const auto& data = m_raw_data.domains_data.at(pv);
  return data.entries[data.entries_back_index].node_level ==
         m_raw_data.current_node_level;
}

std::set<VertexWSM>& DomainsAccessor::get_domain_nonconst(VertexWSM pv) {
  auto& data = m_raw_data.domains_data.at(pv);
  return data.entries[data.entries_back_index].domain;
}

std::set<VertexWSM>&
DomainsAccessor::get_candidate_vertices_for_assignment_nonconst() {
  return m_raw_data.get_current_node_nonconst()
      .pvs_adjacent_to_newly_assigned_vertices;
}
const std::set<VertexWSM>&
DomainsAccessor::get_candidate_vertices_for_assignment() const {
  return m_raw_data.get_current_node().pvs_adjacent_to_newly_assigned_vertices;
}

const std::vector<std::pair<VertexWSM, VertexWSM>>&
DomainsAccessor::get_new_assignments() const {
  return m_raw_data.get_current_node().new_assignments;
}

void DomainsAccessor::clear_new_assignments() {
  m_raw_data.get_current_node_nonconst().new_assignments.clear();
}

WeightWSM DomainsAccessor::get_scalar_product() const {
  return m_raw_data.get_current_node().scalar_product;
}
DomainsAccessor& DomainsAccessor::set_scalar_product(WeightWSM scalar_product) {
  m_raw_data.get_current_node_nonconst().scalar_product = scalar_product;
  return *this;
}

WeightWSM DomainsAccessor::get_total_p_edge_weights() const {
  return m_raw_data.get_current_node().total_p_edge_weights;
}
DomainsAccessor& DomainsAccessor::set_total_p_edge_weights(
    WeightWSM total_weight) {
  m_raw_data.get_current_node_nonconst().total_p_edge_weights = total_weight;
  return *this;
}

bool DomainsAccessor::alldiff_reduce_current_node(
    std::size_t n_assignments_already_processed) {
  auto& node = m_raw_data.get_current_node_nonconst();
  TKET_ASSERT(!node.nogood);
  auto& new_assignments = node.new_assignments;

  while (n_assignments_already_processed < new_assignments.size()) {
    const auto assignment = new_assignments[n_assignments_already_processed];
    ++n_assignments_already_processed;
    node.unassigned_vertices_superset.erase(assignment.first);

    for (auto& domains_information : m_raw_data.domains_data) {
      const VertexWSM& pv = domains_information.first;
      if (pv == assignment.first) {
        continue;
      }
      auto& data_for_this_pv = m_raw_data.domains_data.at(pv);
      auto& existing_domain =
          data_for_this_pv.entries[data_for_this_pv.entries_back_index].domain;
      auto iter = existing_domain.find(assignment.second);
      if (iter == existing_domain.end()) {
        continue;
      }
      if (existing_domain.size() == 1) {
        // Erasing TV would make a nogood.
        return false;
      }
      if (existing_domain.size() == 2) {
        // Erasing TV would make a new assignment.
        VertexWSM tv_other = *existing_domain.cbegin();
        if (tv_other == assignment.second) {
          tv_other = *existing_domain.crbegin();
        }
        new_assignments.emplace_back(pv, tv_other);
      }
      // Now, we've taken care of everything EXCEPT
      // erasing TV from the domain.
      if (data_for_this_pv.entries[data_for_this_pv.entries_back_index]
              .node_level == m_raw_data.current_node_level) {
        // Erase in-place.
        existing_domain.erase(iter);
      } else {
        // We must make a new domain and copy the data across.
        auto& new_entry = get_element_with_resize(
            data_for_this_pv.entries, data_for_this_pv.entries_back_index + 1);
        new_entry.domain =
            data_for_this_pv.entries[data_for_this_pv.entries_back_index]
                .domain;
        TKET_ASSERT(new_entry.domain.erase(assignment.second) == 1);
        new_entry.node_level = m_raw_data.current_node_level;
        ++data_for_this_pv.entries_back_index;
      }
    }
  }
  return true;
}

namespace {

// We need this for the "overwrite_domain" functions.
// Handles everything EXCEPT the final domain overwrite,
// which the caller should do (there are different container types).
struct DomainOverwriteData {
  ReductionResult reduction_result;

  // If true, the caller should just overwrite it in place;
  // all indices, vectors, etc. have been checked and resized if necessary.
  bool caller_should_overwrite_domain;

  template <class DomainObject>
  DomainOverwriteData(
      VertexWSM pattern_v, const DomainObject& new_domain,
      NodesRawData& raw_data) {
    caller_should_overwrite_domain = false;
    if (new_domain.empty()) {
      reduction_result = ReductionResult::NOGOOD;
      return;
    }
    auto& data_for_this_pv = raw_data.domains_data.at(pattern_v);
    auto& existing_domain_data =
        data_for_this_pv.entries[data_for_this_pv.entries_back_index];
    auto& existing_domain = existing_domain_data.domain;
    TKET_ASSERT(!existing_domain.empty());
    TKET_ASSERT(new_domain.size() <= existing_domain.size());

    // We should check that the new domain is a subset;
    // but this is expensive, so just check one element.
    const VertexWSM first_new_tv = *new_domain.cbegin();
    TKET_ASSERT(existing_domain.count(first_new_tv) != 0);
    if (new_domain.size() == existing_domain.size()) {
      reduction_result = ReductionResult::SUCCESS;
      return;
    }
    caller_should_overwrite_domain = true;

    // Now, the new domain is nonempty and smaller than the existing one.
    if (new_domain.size() == 1) {
      reduction_result = ReductionResult::NEW_ASSIGNMENTS;
      raw_data.get_current_node_nonconst().new_assignments.emplace_back(
          pattern_v, first_new_tv);
      raw_data.get_current_node_nonconst().unassigned_vertices_superset.erase(
          pattern_v);
    } else {
      reduction_result = ReductionResult::SUCCESS;
    }

    if (existing_domain_data.node_level == raw_data.current_node_level) {
      // We can just overwrite the domain in-place; the caller will do that.
      return;
    }
    // We need to make a new domain object.
    ++data_for_this_pv.entries_back_index;
    auto& new_entry = get_element_with_resize(
        data_for_this_pv.entries, data_for_this_pv.entries_back_index);
    new_entry.node_level = raw_data.current_node_level;
  }
};
}  // namespace

ReductionResult DomainsAccessor::overwrite_domain(
    VertexWSM pv, const std::set<VertexWSM>& new_domain) {
  const DomainOverwriteData overwrite_data(pv, new_domain, m_raw_data);
  if (overwrite_data.caller_should_overwrite_domain) {
    // Simple copy.
    get_domain_nonconst(pv) = new_domain;
  }
  return overwrite_data.reduction_result;
}

ReductionResult DomainsAccessor::overwrite_domain_with_set_swap(
    VertexWSM pv, std::set<VertexWSM>& new_domain) {
  const DomainOverwriteData overwrite_data(pv, new_domain, m_raw_data);
  if (overwrite_data.caller_should_overwrite_domain) {
    get_domain_nonconst(pv).swap(new_domain);
  }
  return overwrite_data.reduction_result;
}

ReductionResult DomainsAccessor::overwrite_domain(
    VertexWSM pv, const std::vector<VertexWSM>& new_domain) {
  const DomainOverwriteData overwrite_data(pv, new_domain, m_raw_data);
  if (overwrite_data.caller_should_overwrite_domain) {
    get_domain_nonconst(pv) = {new_domain.cbegin(), new_domain.cend()};
  }
  return overwrite_data.reduction_result;
}

DomainsAccessor::IntersectionResult
DomainsAccessor::intersect_domain_with_complement_set(
    VertexWSM pattern_v, const std::set<VertexWSM>& forbidden_target_vertices) {
  IntersectionResult result;
  {
    const auto& current_domain = get_domain(pattern_v);
    if (disjoint(current_domain, forbidden_target_vertices)) {
      result.changed = false;
      result.new_domain_size = current_domain.size();
      result.reduction_result = ReductionResult::SUCCESS;
      return result;
    }
  }
  // The domain definitely is changing.
  result.changed = true;
  auto& data_for_this_pv = m_raw_data.domains_data.at(pattern_v);
  if (data_for_this_pv.entries[data_for_this_pv.entries_back_index]
          .node_level == m_raw_data.current_node_level) {
    // We can simply overwrite.
    auto& domain =
        data_for_this_pv.entries[data_for_this_pv.entries_back_index].domain;
    for (VertexWSM tv : forbidden_target_vertices) {
      if (domain.erase(tv) == 1 && domain.empty()) {
        result.reduction_result = ReductionResult::NOGOOD;
        return result;
      }
    }
  } else {
    ++data_for_this_pv.entries_back_index;
    resize_if_index_is_invalid(
        data_for_this_pv.entries, data_for_this_pv.entries_back_index);
    auto& entry = data_for_this_pv.entries[data_for_this_pv.entries_back_index];
    entry.node_level = m_raw_data.current_node_level;

    // Is it quicker to build it up as here, or copy then grind it down?
    // Depends on all the sizes, and how many end up being erased...
    entry.domain.clear();
    for (VertexWSM old_tv :
         data_for_this_pv.entries[data_for_this_pv.entries_back_index - 1]
             .domain) {
      if (forbidden_target_vertices.count(old_tv) == 0) {
        entry.domain.insert(old_tv);
      }
    }
    if (entry.domain.empty()) {
      // Erase it again!
      --data_for_this_pv.entries_back_index;
      result.reduction_result = ReductionResult::NOGOOD;
      return result;
    }
  }
  // At this stage, we've created the new domain
  // (necessarily different from the old one).
  const auto& domain =
      data_for_this_pv.entries[data_for_this_pv.entries_back_index].domain;
  result.new_domain_size = domain.size();
  TKET_ASSERT(!domain.empty());
  if (domain.size() > 1) {
    result.reduction_result = ReductionResult::SUCCESS;
  } else {
    m_raw_data.get_current_node_nonconst().new_assignments.emplace_back(
        pattern_v, *domain.cbegin());
    result.reduction_result = ReductionResult::NEW_ASSIGNMENTS;
  }
  return result;
}

std::set<VertexWSM>& DomainsAccessor::
    get_current_node_unassigned_pattern_vertices_superset_to_overwrite() {
  return m_raw_data.get_current_node_nonconst().unassigned_vertices_superset;
}

std::string DomainsAccessor::str(bool full) const {
  std::stringstream ss;
  if (full) {
    ss << "\n@@@@@ ALL NODES: (curr.lev=" << m_raw_data.current_node_level
       << ")";
    for (unsigned level = 0; level <= m_raw_data.current_node_level; ++level) {
      const auto& node = m_raw_data.nodes_data[level];
      ss << "\n+++ node " << level << ":" << node.str();
    }
    ss << "\nDOMAINS: ";
    for (const auto& entry : m_raw_data.domains_data) {
      ss << "\nDOM(" << entry.first << "):" << entry.second.str();
    }
    ss << "\n";
    return ss.str();
  }
  ss << "\ncurr.node lev=" << m_raw_data.current_node_level << "; DOMAINS: ";
  for (const auto& entry : m_raw_data.domains_data) {
    ss << "\n  DOM(" << entry.first << "): ";
    const auto& dom_data =
        entry.second.entries[entry.second.entries_back_index];
    ss << "(since lev " << dom_data.node_level
       << "): " << tket::WeightedSubgraphMonomorphism::str(dom_data.domain);
  }
  ss << "\n";
  return ss.str();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
