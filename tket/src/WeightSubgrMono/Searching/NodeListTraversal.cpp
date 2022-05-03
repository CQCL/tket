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

#include "WeightSubgrMono/Searching/NodeListTraversal.hpp"

#include <sstream>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Searching/NodesRawData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NodeListTraversal::NodeListTraversal(NodesRawDataWrapper& raw_data_wrapper)
    : m_raw_data(raw_data_wrapper.m_raw_data) {}

std::set<VertexWSM> NodeListTraversal::get_used_target_vertices() const {
  std::set<VertexWSM> target_vertices;

  // Examine all PV.
  for (const auto& entry_just_for_one_pv : m_raw_data.domains_data) {
    const VertexWSM& pv = entry_just_for_one_pv.first;
    const auto& domain_data = entry_just_for_one_pv.second;

    // Examine Dom(PV) at all levels.
    for (unsigned level = 0; level <= domain_data.entries_back_index; ++level) {
      const auto& node_index = domain_data.entries[level].node_level;
      // Nogood nodes are ignored.
      // So, check that at least one node with this domain
      // is not a nogood; otherwise skip.
      if (m_raw_data.nodes_data[node_index].nogood &&
          level < domain_data.entries_back_index) {
        const auto& next_node_index = domain_data.entries[level + 1].node_level;
        bool some_node_is_good = false;
        for (auto ii = node_index; ii < next_node_index; ++ii) {
          if (!m_raw_data.nodes_data[ii].nogood) {
            some_node_is_good = true;
            break;
          }
        }
        if (!some_node_is_good) {
          // Skip.
          continue;
        }
      }
      const auto& domain = domain_data.entries[level].domain;
      if (target_vertices.empty()) {
        target_vertices = domain;
      } else {
        for (VertexWSM tv : domain) {
          target_vertices.insert(tv);
        }
      }
    }
  }
  return target_vertices;
}

bool NodeListTraversal::move_up() {
  while (m_raw_data.current_node_level > 0) {
    --m_raw_data.current_node_level;
    if (m_raw_data.nodes_data[m_raw_data.current_node_level].nogood) {
      continue;
    }
    for (auto& entry_for_single_pv : m_raw_data.domains_data) {
      auto& domain_data = entry_for_single_pv.second;
      while (domain_data.entries[domain_data.entries_back_index].node_level >
             m_raw_data.current_node_level) {
        // We've moved above the level of "junk data",
        // so shrink the back index.
        TKET_ASSERT(domain_data.entries_back_index > 0);
        --domain_data.entries_back_index;
      }
    }
    return true;
  }
  // We've hit the top!
  return false;
}

// Returns false if the current node becomes invalid.
// Also handles the new assignment PV->y, if necessary,
// for the CURRENT node (not new node).
static bool when_moving_down_check_current_domain_size(
    const std::set<VertexWSM>& existing_domain, NodesRawData& raw_data,
    VertexWSM p_vertex, VertexWSM t_vertex) {
  switch (existing_domain.size()) {
    case 0:
      TKET_ASSERT(false);
      return false;
    case 1:
      // Already Dom(PV) = {TV}, so it will become a nogood.
      raw_data.get_current_node_nonconst().nogood = true;
      return false;
    case 2: {
      // Dom(PV) = {TV, y}, so it will become Dom(PV) = {y}
      // at the current level, so we need a new assignment PV->y.
      VertexWSM tv_other = *existing_domain.cbegin();
      if (tv_other == t_vertex) {
        tv_other = *existing_domain.crbegin();
      }
      raw_data.get_current_node_nonconst().new_assignments.emplace_back(
          p_vertex, tv_other);
      break;
    }
    default:
      break;
  }
  return true;
}

// Only used when the old domain is shared with previous nodes,
// and hence must be copied; but we've also guaranteed valid vector sizes.
static void when_moving_down_copy_old_shared_domain_and_erase_tv(
    NodesRawData& raw_data, NodesRawData::DomainData& data_for_this_pv,
    VertexWSM t_vertex) {
  data_for_this_pv.entries[data_for_this_pv.entries_back_index + 1].node_level =
      raw_data.current_node_level;

  auto& new_domain =
      data_for_this_pv.entries[data_for_this_pv.entries_back_index + 1].domain;
  new_domain =
      data_for_this_pv.entries[data_for_this_pv.entries_back_index].domain;
  TKET_ASSERT(new_domain.erase(t_vertex) == 1);
}

// data_for_this_pv.entries_back_index should now be for
// the new domain we'll create.
static void complete_move_down_with_resized_vectors_and_indices(
    NodesRawData& raw_data, NodesRawData::DomainData& data_for_this_pv,
    VertexWSM p_vertex, VertexWSM t_vertex) {
  // Remember, the domain data and node data
  // at the current level is now "junk".
  auto& new_domain =
      data_for_this_pv.entries[data_for_this_pv.entries_back_index].domain;
  new_domain.clear();
  new_domain.insert(t_vertex);

  data_for_this_pv.entries[data_for_this_pv.entries_back_index].node_level =
      raw_data.current_node_level;

  auto& new_node_data = raw_data.nodes_data[raw_data.current_node_level];
  new_node_data.new_assignments.resize(1);
  new_node_data.new_assignments[0] = {p_vertex, t_vertex};
  new_node_data.nogood = false;
  new_node_data.pvs_adjacent_to_newly_assigned_vertices.clear();
  new_node_data.scalar_product =
      raw_data.nodes_data[raw_data.current_node_level - 1].scalar_product;
  new_node_data.total_p_edge_weights =
      raw_data.nodes_data[raw_data.current_node_level - 1].total_p_edge_weights;
}

void NodeListTraversal::move_down(VertexWSM p_vertex, VertexWSM t_vertex) {
  // We can only move down from a valid, fully reduced node.
  bool existing_node_valid = !m_raw_data.get_current_node().nogood;
  TKET_ASSERT(existing_node_valid);
  TKET_ASSERT(m_raw_data.get_current_node().new_assignments.empty());

  auto& data_for_this_pv = m_raw_data.domains_data.at(p_vertex);
  auto& existing_domain =
      data_for_this_pv.entries[data_for_this_pv.entries_back_index].domain;
  auto iter = existing_domain.find(t_vertex);

  // TV must be present!
  TKET_ASSERT(iter != existing_domain.end());
  existing_node_valid = when_moving_down_check_current_domain_size(
      existing_domain, m_raw_data, p_vertex, t_vertex);

  // Next, we need to erase TV from the existing domain.
  if (existing_node_valid) {
    if (data_for_this_pv.entries[data_for_this_pv.entries_back_index]
            .node_level == m_raw_data.current_node_level) {
      // The data is not shared by previous nodes, so we can just overwrite it
      // in place.
      existing_domain.erase(iter);

      // Don't forget the new singleton domain   Dom(PV) = {TV}
      // we'll create shortly!
      resize_if_index_is_invalid(
          data_for_this_pv.entries, data_for_this_pv.entries_back_index + 1);

      ++data_for_this_pv.entries_back_index;
    } else {
      resize_if_index_is_invalid(
          data_for_this_pv.entries, data_for_this_pv.entries_back_index + 2);

      when_moving_down_copy_old_shared_domain_and_erase_tv(
          m_raw_data, data_for_this_pv, t_vertex);

      data_for_this_pv.entries_back_index += 2;
    }
  } else {
    // The current node becomes invalid, so don't waste any time with domains.
    // (Note that the existing shared data is NOT changed - whether it points
    // to the current node or not is now irrelevant, as the current node
    // will simply be ignored in future).
    // HOWEVER we still need to resize, ready for the new domain.
    resize_if_index_is_invalid(
        data_for_this_pv.entries, data_for_this_pv.entries_back_index + 1);

    ++data_for_this_pv.entries_back_index;
  }
  ++m_raw_data.current_node_level;
  resize_if_index_is_invalid(
      m_raw_data.nodes_data, m_raw_data.current_node_level);
  m_raw_data.get_current_node_nonconst().unassigned_vertices_superset.clear();
  complete_move_down_with_resized_vectors_and_indices(
      m_raw_data, data_for_this_pv, p_vertex, t_vertex);
}

bool NodeListTraversal::erase_impossible_assignment(
    std::pair<VertexWSM, VertexWSM> impossible_assignment,
    ImpossibleAssignmentAction action) {
  auto max_level = m_raw_data.current_node_level;
  if (action == ImpossibleAssignmentAction::PROCESS_CURRENT_NODE) {
    if (max_level == 0) {
      return true;
    }
    --max_level;
  }
  auto& data_for_this_pv =
      m_raw_data.domains_data.at(impossible_assignment.first);

  for (unsigned ii = 0; ii <= data_for_this_pv.entries_back_index; ++ii) {
    TKET_ASSERT(
        data_for_this_pv.entries[ii].node_level <=
        m_raw_data.current_node_level);
    if (data_for_this_pv.entries[ii].node_level > max_level) {
      break;
    }
    auto& domain = data_for_this_pv.entries[ii].domain;
    if (domain.erase(impossible_assignment.second) == 0 || domain.size() >= 2) {
      // Nothing changed, OR it did change, but with no significant effect.
      continue;
    }
    // Now, it's EITHER a nogood, OR a new assignment is created.
    // Anyway, we have to go through all the nodes which use this domain.
    unsigned node_level_end = max_level + 1;
    if (ii < data_for_this_pv.entries_back_index) {
      node_level_end =
          std::min(node_level_end, data_for_this_pv.entries[ii + 1].node_level);
    }
    if (domain.size() == 0) {
      for (unsigned jj = data_for_this_pv.entries[ii].node_level;
           jj < node_level_end; ++jj) {
        m_raw_data.nodes_data[jj].nogood = true;
      }
    } else {
      // A new assignment.
      const auto new_assignment =
          std::make_pair(impossible_assignment.first, *domain.cbegin());

      for (unsigned jj = data_for_this_pv.entries[ii].node_level;
           jj < node_level_end; ++jj) {
        m_raw_data.nodes_data[jj].new_assignments.emplace_back(new_assignment);
      }
    }
  }
  if (action == ImpossibleAssignmentAction::PROCESS_CURRENT_NODE) {
    return true;
  }
  return !m_raw_data.get_current_node().nogood;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
