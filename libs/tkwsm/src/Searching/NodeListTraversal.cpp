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

#include "tkwsm/Searching/NodeListTraversal.hpp"

#include <sstream>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Searching/NodesRawData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NodeListTraversal::NodeListTraversal(NodesRawDataWrapper& raw_data_wrapper)
    : m_raw_data(raw_data_wrapper.m_raw_data) {}

// Given the index of a valid entry in domain_data_for_pv.entries,
// find the last node index which we need to check,
// for all the nodes that share this domain.
static unsigned get_final_node_index_for_shared_domain(
    const NodesRawData::DomainData& domain_data_for_pv, unsigned entries_index,
    const NodesRawData& raw_data) {
  TKET_ASSERT(entries_index < domain_data_for_pv.entries.size());
  if (entries_index + 1 == domain_data_for_pv.entries.size()) {
    // It's the final domain in the list.
    return raw_data.current_node_index();
  }
  return domain_data_for_pv.entries[entries_index + 1].node_index - 1;
}

std::set<VertexWSM> NodeListTraversal::get_used_target_vertices() const {
  std::set<VertexWSM> target_vertices;

  // Examine all PV.
  for (unsigned pv = 0; pv < m_raw_data.domains_data.size(); ++pv) {
    const NodesRawData::DomainData& domain_data = m_raw_data.domains_data[pv];

    // Add Dom(PV) at all valid nodes.
    const unsigned domain_data_entries_size = domain_data.entries.size();
    for (unsigned entries_index = 0; entries_index < domain_data_entries_size;
         ++entries_index) {
      const unsigned& node_index =
          domain_data.entries[entries_index].node_index;

      // Nogood nodes are ignored.
      // So, check that at least one node sharing this domain
      // is not a nogood; otherwise skip.
      if (m_raw_data.nodes_data[node_index].nogood) {
        bool some_good_node_shares_this_domain = false;
        const unsigned final_node_index =
            get_final_node_index_for_shared_domain(
                domain_data, entries_index, m_raw_data);
        for (unsigned ii = node_index + 1; ii <= final_node_index; ++ii) {
          if (!m_raw_data.nodes_data[ii].nogood) {
            some_good_node_shares_this_domain = true;
            break;
          }
        }
        if (!some_good_node_shares_this_domain) {
          // Ignore this domain: all nodes sharing it are nogood.
          continue;
        }
      }
      const std::set<VertexWSM>& domain =
          domain_data.entries[entries_index].domain;
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
  if (m_raw_data.nodes_data.size() <= 1) {
    // We will run out of nodes when we pop.
    // So, this is the end of the search!
    return false;
  }
  m_raw_data.nodes_data.pop();

  // Keep popping more nodes off, until we reach a node which is NOT a nogood
  // (i.e., a "good" node; although not standard terminology!)
  while (m_raw_data.nodes_data.top().nogood) {
    m_raw_data.nodes_data.pop();
    if (m_raw_data.nodes_data.empty()) {
      return false;
    }
  }

  // Now we've reached a "good" (not nogood) node.
  const unsigned current_node_index = m_raw_data.current_node_index();

  for (unsigned pv = 0; pv < m_raw_data.domains_data.size(); ++pv) {
    NodesRawData::DomainData& domain_data = m_raw_data.domains_data[pv];
    while (domain_data.entries.top().node_index > current_node_index) {
      // Data for any higher node_index is simply junk.
      domain_data.entries.pop();
      TKET_ASSERT(!domain_data.entries.empty());
    }
  }
  return true;
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
      // at this node, so we need a new assignment PV->y.
      // The elements of a set of size 2 can be found using
      // cbegin and crbegin.
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

// We've taken care of erasing the TV from the existing domain,
// now we must create the new node with the new assignment PV->TV.
static void complete_move_down(
    NodesRawData& raw_data, NodesRawData::DomainData& data_for_this_pv,
    VertexWSM p_vertex, VertexWSM t_vertex) {
  // The newly created node is "junk".
  raw_data.nodes_data.push();

  // Recall that "empty" actually means "check a previous vector".
  raw_data.get_current_node_nonconst().unassigned_vertices_superset.clear();

  auto& new_node_data = raw_data.nodes_data.top();
  new_node_data.new_assignments.resize(1);
  new_node_data.new_assignments[0] = {p_vertex, t_vertex};
  new_node_data.nogood = false;

  new_node_data.scalar_product =
      raw_data.nodes_data.one_below_top().scalar_product;

  new_node_data.total_p_edge_weights =
      raw_data.nodes_data.one_below_top().total_p_edge_weights;

  // We must also create the domain Dom(PV) = {TV}.
  data_for_this_pv.entries.push();
  data_for_this_pv.entries.top().node_index = raw_data.current_node_index();

  // The domain data for the current node is now "junk".
  std::set<VertexWSM>& new_domain = data_for_this_pv.entries.top().domain;
  new_domain.clear();
  new_domain.insert(t_vertex);
}

void NodeListTraversal::move_down(VertexWSM p_vertex, VertexWSM t_vertex) {
  // We can only move down from a valid, fully reduced node.
  bool existing_node_valid = !m_raw_data.get_current_node().nogood;
  TKET_ASSERT(existing_node_valid);
  TKET_ASSERT(m_raw_data.get_current_node().new_assignments.empty());

  NodesRawData::DomainData& data_for_this_pv =
      m_raw_data.domains_data.at(p_vertex);
  std::set<VertexWSM>& existing_domain = data_for_this_pv.entries.top().domain;
  auto iter = existing_domain.find(t_vertex);

  // TV must be present!
  TKET_ASSERT(iter != existing_domain.end());
  existing_node_valid = when_moving_down_check_current_domain_size(
      existing_domain, m_raw_data, p_vertex, t_vertex);

  // Next, we need to erase TV from the existing domain.
  if (existing_node_valid) {
    if (data_for_this_pv.entries.top().node_index ==
        m_raw_data.current_node_index()) {
      // The data is not shared by previous nodes, so we can just overwrite it
      // in place.
      existing_domain.erase(iter);
    } else {
      data_for_this_pv.entries.push();
      data_for_this_pv.entries.top().node_index =
          m_raw_data.current_node_index();
      data_for_this_pv.entries.top().domain =
          data_for_this_pv.entries.one_below_top().domain;
      TKET_ASSERT(data_for_this_pv.entries.top().domain.erase(t_vertex) == 1);
    }
  }

  complete_move_down(m_raw_data, data_for_this_pv, p_vertex, t_vertex);
}

static void fill_nogood_or_new_assignment_in_all_shared_nodes(
    VertexWSM pv, const std::set<VertexWSM>& new_domain,
    unsigned data_for_this_pv_entry_index,
    const NodesRawData::DomainData& data_for_this_pv, NodesRawData& raw_data) {
  // Which nodes share the domain?
  const unsigned final_node_index = get_final_node_index_for_shared_domain(
      data_for_this_pv, data_for_this_pv_entry_index, raw_data);

  if (new_domain.empty()) {
    // The nodes are all nogoods.
    for (unsigned jj =
             data_for_this_pv.entries[data_for_this_pv_entry_index].node_index;
         jj <= final_node_index; ++jj) {
      raw_data.nodes_data[jj].nogood = true;
    }
    return;
  }

  // It's the same new assignment PV->TV to be added to all nodes.
  const auto new_assignment = std::make_pair(pv, *new_domain.cbegin());

  for (unsigned jj =
           data_for_this_pv.entries[data_for_this_pv_entry_index].node_index;
       jj <= final_node_index; ++jj) {
    // Of course, some shared nodes might be nogoods for some OTHER reason,
    // related to other domains we don't know about.
    if (!raw_data.nodes_data[jj].nogood) {
      raw_data.nodes_data[jj].new_assignments.emplace_back(new_assignment);
    }
  }
}

void NodeListTraversal::erase_impossible_assignment(
    std::pair<VertexWSM, VertexWSM> impossible_assignment) {
  auto& data_for_this_pv =
      m_raw_data.domains_data.at(impossible_assignment.first);

  const unsigned size = data_for_this_pv.entries.size();
  for (unsigned ii = 0; ii < size; ++ii) {
    TKET_ASSERT(
        data_for_this_pv.entries[ii].node_index <=
        m_raw_data.current_node_index());

    if (m_raw_data.nodes_data[data_for_this_pv.entries[ii].node_index].nogood) {
      // Don't waste time with nogood nodes.
      continue;
    }
    std::set<VertexWSM>& domain = data_for_this_pv.entries[ii].domain;
    if (domain.erase(impossible_assignment.second) == 0 || domain.size() >= 2) {
      // EITHER nothing changed, OR the change had no significant effect.
      continue;
    }
    // Now, it's EITHER a nogood, OR a new assignment is created.
    // Anyway, we go through all the nodes which share this domain and
    // process them (the same Domain(pv) may be shared across several nodes).
    fill_nogood_or_new_assignment_in_all_shared_nodes(
        impossible_assignment.first, domain, ii, data_for_this_pv, m_raw_data);
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
