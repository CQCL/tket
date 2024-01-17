// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/Searching/DomainsAccessor.hpp"

#include <sstream>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Common/TemporaryRefactorCode.hpp"
#include "tkwsm/Searching/NodesRawData.hpp"

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

const boost::dynamic_bitset<>& DomainsAccessor::get_domain(VertexWSM pv) const {
  return m_raw_data.domains_data.at(pv).entries.top().domain;
}

std::size_t DomainsAccessor::get_domain_size(VertexWSM pv) const {
  return m_raw_data.domains_data.at(pv).entries.top().domain.count();
}

bool DomainsAccessor::domain_created_in_current_node(VertexWSM pv) const {
  return m_raw_data.domains_data.at(pv).entries.top().node_index ==
         m_raw_data.current_node_index();
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
      auto& existing_domain_bitset = data_for_this_pv.entries.top().domain;

      // Check if erasing TV would make a nogood or assignment.
      // So we need to know a few vertices tv1, tv2, ...
      {
        const auto tv1 = existing_domain_bitset.find_first();
        // It must currently be nonempty.
        TKET_ASSERT(tv1 < existing_domain_bitset.size());

        if (!existing_domain_bitset.test(assignment.second)) {
          // TV is not present in the domain; nothing to change.
          continue;
        }
        // Does it have 1 or 2 elements currently?
        const auto tv2 = existing_domain_bitset.find_next(tv1);
        if (tv2 < existing_domain_bitset.size()) {
          // It has at least 2 vertices. Does it have another?
          const auto tv3 = existing_domain_bitset.find_next(tv2);
          if (tv3 >= existing_domain_bitset.size()) {
            // It has EXACTLY 2 vertices: tv1, tv2.
            // One of them must be the TV we're erasing.
            // The other will form a new assignment.
            VertexWSM tv_other = tv1;
            if (tv_other == assignment.second) {
              tv_other = tv2;
            }
            TKET_ASSERT(tv_other != assignment.second);
            new_assignments.emplace_back(pv, tv_other);
          }
        } else {
          // It has EXACTLY one vertex: tv1.
          // It MUST equal TV, and then erasing it would make a nogood.
          TKET_ASSERT(tv1 == assignment.second);
          return false;
        }
      }

      // Now, we've taken care of everything EXCEPT
      // erasing TV from the domain (which we KNOW is present).
      // But we cannot immediately erase, since the domain might be shared
      // across several nodes.
      if (data_for_this_pv.entries.top().node_index ==
          m_raw_data.current_node_index()) {
        // Erase in-place.
        TKET_ASSERT(existing_domain_bitset.test_set(assignment.second, false));
      } else {
        // We must make a new domain and copy the data across.
        data_for_this_pv.entries.push();
        data_for_this_pv.entries.top().node_index =
            m_raw_data.current_node_index();
        data_for_this_pv.entries.top().domain =
            data_for_this_pv.entries.one_below_top().domain;
        TKET_ASSERT(data_for_this_pv.entries.top().domain.test_set(
            assignment.second, false));
      }
    }
  }
  return true;
}

// TODO: make another version NOT using a swap!
DomainsAccessor::IntersectionResult DomainsAccessor::intersect_domain_with_swap(
    VertexWSM pv, boost::dynamic_bitset<>& domain_mask) {
  auto& data_for_this_pv = m_raw_data.domains_data.at(pv);
  domain_mask &= get_domain(pv);

  IntersectionResult result;
  result.new_domain_size = domain_mask.count();
  result.changed = get_domain_size(pv) != result.new_domain_size;
  if (!result.changed) {
    result.reduction_result = ReductionResult::SUCCESS;
    return result;
  }
  if (result.new_domain_size == 0) {
    result.reduction_result = ReductionResult::NOGOOD;
    return result;
  }
  if (result.new_domain_size == 1) {
    result.reduction_result = ReductionResult::NEW_ASSIGNMENTS;
    auto& current_node = m_raw_data.get_current_node_nonconst();
    current_node.new_assignments.emplace_back(
        pv, VertexWSM(domain_mask.find_first()));
  } else {
    result.reduction_result = ReductionResult::SUCCESS;
  }

  if (data_for_this_pv.entries.top().node_index !=
      m_raw_data.current_node_index()) {
    // We need to make a new domain object; it's shared.
    data_for_this_pv.entries.push();
    data_for_this_pv.entries.top().node_index = m_raw_data.current_node_index();
  }
  data_for_this_pv.entries.top().domain.swap(domain_mask);
  return result;
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
