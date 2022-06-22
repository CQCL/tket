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
#include <tkassert/Assert.hpp>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/Searching/NodesRawData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DomainsAccessor::DomainsAccessor(NodesRawDataWrapper& raw_data_wrapper)
    : m_raw_data(raw_data_wrapper.m_raw_data) {}

unsigned DomainsAccessor::get_number_of_pattern_vertices() const {
  return m_raw_data.domains_data.size();
}

bool DomainsAccessor::current_node_is_valid() const {
  return !m_raw_data.nodes_data.top().nogood;
}

const std::set<VertexWSM>& DomainsAccessor::get_domain(VertexWSM pv) const {
  const auto& data = m_raw_data.domains_data.at(pv);
  return data.entries.top().domain;
}

bool DomainsAccessor::domain_created_in_current_node(VertexWSM pv) const {
  return m_raw_data.domains_data.at(pv).entries.top().node_index ==
         m_raw_data.current_node_index();
}

std::set<VertexWSM>& DomainsAccessor::get_domain_nonconst(VertexWSM pv) {
  return m_raw_data.domains_data.at(pv).entries.top().domain;
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
  NodesRawData::NodeData& node = m_raw_data.get_current_node_nonconst();
  TKET_ASSERT(!node.nogood);
  std::vector<std::pair<VertexWSM, VertexWSM>>& new_assignments =
      node.new_assignments;

  while (n_assignments_already_processed < new_assignments.size()) {
    const std::pair<VertexWSM, VertexWSM> assignment =
        new_assignments[n_assignments_already_processed];
    ++n_assignments_already_processed;

    for (unsigned pv = 0; pv < m_raw_data.domains_data.size(); ++pv) {
      if (pv == assignment.first) {
        continue;
      }
      NodesRawData::DomainData& data_for_this_pv = m_raw_data.domains_data[pv];
      std::set<VertexWSM>& existing_domain =
          data_for_this_pv.entries.top().domain;
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
      if (data_for_this_pv.entries.top().node_index ==
          m_raw_data.current_node_index()) {
        // Erase in-place.
        existing_domain.erase(iter);
      } else {
        // We must make a new domain and copy the data across.
        data_for_this_pv.entries.push();
        data_for_this_pv.entries.top().node_index =
            m_raw_data.current_node_index();
        data_for_this_pv.entries.top().domain =
            data_for_this_pv.entries.one_below_top().domain;
        TKET_ASSERT(
            data_for_this_pv.entries.top().domain.erase(assignment.second) ==
            1);
      }
    }
  }
  return true;
}

ReductionResult DomainsAccessor::overwrite_domain_with_set_swap(
    VertexWSM pv, std::set<VertexWSM>& new_domain) {
  if (new_domain.empty()) {
    return ReductionResult::NOGOOD;
  }
  auto& data_for_this_pv = m_raw_data.domains_data.at(pv);
  auto& existing_domain_data = data_for_this_pv.entries.top();
  std::set<VertexWSM>& existing_domain = existing_domain_data.domain;

  // We should do a detailed assert that "new_domain" is a subset
  // of "existing_domain", but this is expensive; thus we do a very
  // cheap assert.
  TKET_ASSERT(existing_domain.size() >= new_domain.size());
  TKET_ASSERT(existing_domain.count(*new_domain.cbegin()) != 0);

  if (existing_domain.size() == new_domain.size()) {
    // No change, no action required (even if the new domain has size 1;
    // it's not a NEW assignment).
    return ReductionResult::SUCCESS;
  }
  auto result = ReductionResult::SUCCESS;

  if (new_domain.size() == 1) {
    result = ReductionResult::NEW_ASSIGNMENTS;
    auto& current_node = m_raw_data.get_current_node_nonconst();
    current_node.new_assignments.emplace_back(pv, *new_domain.cbegin());
  }
  // We must be careful; "existing_domain" might be shared with previous nodes.
  if (existing_domain_data.node_index == m_raw_data.current_node_index()) {
    // We can just overwrite the domain in-place; the caller will do that.
    existing_domain.swap(new_domain);
  } else {
    // We need to make a new domain object.
    data_for_this_pv.entries.push();
    data_for_this_pv.entries.top().node_index = m_raw_data.current_node_index();
    data_for_this_pv.entries.top().domain.swap(new_domain);
  }
  return result;
}

DomainsAccessor::IntersectionResult
DomainsAccessor::intersect_domain_with_complement_set(
    VertexWSM pattern_v, const std::set<VertexWSM>& forbidden_target_vertices) {
  IntersectionResult result;
  {
    const std::set<VertexWSM>& current_domain = get_domain(pattern_v);

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
  if (data_for_this_pv.entries.top().node_index ==
      m_raw_data.current_node_index()) {
    // Domain is not shared; we can simply overwrite.
    std::set<VertexWSM>& domain = data_for_this_pv.entries.top().domain;

    for (VertexWSM tv : forbidden_target_vertices) {
      if (domain.erase(tv) == 1 && domain.empty()) {
        result.reduction_result = ReductionResult::NOGOOD;
        return result;
      }
    }
  } else {
    // Shared domain; copy on write:
    data_for_this_pv.entries.push();
    data_for_this_pv.entries.top().node_index = m_raw_data.current_node_index();

    const std::set<VertexWSM>& old_domain =
        data_for_this_pv.entries.one_below_top().domain;
    std::set<VertexWSM>& new_domain = data_for_this_pv.entries.top().domain;

    // Is it quicker to build it up as here, or copy then grind it down?
    // Depends on all the sizes, and how many end up being erased...
    new_domain.clear();

    for (VertexWSM old_tv : old_domain) {
      if (forbidden_target_vertices.count(old_tv) == 0) {
        new_domain.insert(old_tv);
      }
    }

    if (new_domain.empty()) {
      // The caller will move up, no further action here.
      result.reduction_result = ReductionResult::NOGOOD;
      return result;
    }
  }

  // At this stage, we've created and filled the new domain object
  // (necessarily different from the old one).
  const std::set<VertexWSM>& domain = data_for_this_pv.entries.top().domain;
  result.new_domain_size = domain.size();
  TKET_ASSERT(!domain.empty());
  if (domain.size() > 1) {
    result.reduction_result = ReductionResult::SUCCESS;
  } else {
    // domain size 1, so a new assignment (as we KNOW that the domain has
    // changed).
    m_raw_data.get_current_node_nonconst().new_assignments.emplace_back(
        pattern_v, *domain.cbegin());
    result.reduction_result = ReductionResult::NEW_ASSIGNMENTS;
  }
  return result;
}

std::string DomainsAccessor::str(bool full) const {
  std::stringstream ss;
  if (full) {
    ss << "\n@@@@@ ALL NODES: (curr.node index="
       << m_raw_data.current_node_index() << ")";
    for (unsigned node_index = 0; node_index < m_raw_data.nodes_data.size();
         ++node_index) {
      const auto& node = m_raw_data.nodes_data[node_index];
      ss << "\n+++ node " << node_index << ":" << node.str();
    }
    ss << "\nDOMAINS: ";
    for (unsigned pv = 0; pv < m_raw_data.domains_data.size(); ++pv) {
      ss << "\nDOM(" << pv << "):" << m_raw_data.domains_data[pv].str();
    }
    ss << "\n";
    return ss.str();
  }
  ss << "\ncurr.node index=" << m_raw_data.current_node_index()
     << "; DOMAINS: ";
  for (unsigned pv = 0; pv < m_raw_data.domains_data.size(); ++pv) {
    ss << "\n  DOM(" << pv << "): ";
    const auto& dom_data = m_raw_data.domains_data[pv].entries.top();
    ss << "(since index " << dom_data.node_index
       << "): " << tket::WeightedSubgraphMonomorphism::str(dom_data.domain);
  }
  ss << "\n";
  return ss.str();
}

const std::vector<VertexWSM>&
DomainsAccessor::get_unassigned_pattern_vertices_superset() const {
  const auto& candidate =
      m_raw_data.get_current_node().unassigned_vertices_superset;
  if (!candidate.empty()) {
    return candidate;
  }
  TKET_ASSERT(m_raw_data.nodes_data.size() > 1);
  return m_raw_data.nodes_data.one_below_top().unassigned_vertices_superset;
}

std::vector<VertexWSM>&
DomainsAccessor::get_unassigned_pattern_vertices_superset_to_overwrite() {
  return m_raw_data.get_current_node_nonconst().unassigned_vertices_superset;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
