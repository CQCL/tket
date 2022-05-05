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

#include "WeightSubgrMono/Searching/NodesRawData.hpp"

#include <sstream>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static std::vector<VertexWSM> get_pattern_vertices(
    const PossibleAssignments& possible_assignments) {
  std::vector<VertexWSM> pattern_vertices;
  pattern_vertices.reserve(possible_assignments.size());
  for (const auto& entry : possible_assignments) {
    pattern_vertices.push_back(entry.first);
  }
  return pattern_vertices;
}

NodesRawData::NodesRawData(const PossibleAssignments& possible_assignments)
    : pattern_vertices(get_pattern_vertices(possible_assignments)),
      current_node_level(0) {
  nodes_data.resize(1);

  auto& node = nodes_data[0];
  node.nogood = false;
  node.scalar_product = 0;
  node.total_p_edge_weights = 0;

  for (const auto& entry : possible_assignments) {
    const VertexWSM& pv = entry.first;
    const auto& domain = entry.second;

    auto& domain_data = domains_data[pv];
    domain_data.entries_back_index = 0;
    domain_data.entries.resize(1);
    domain_data.entries[0].domain = domain;
    domain_data.entries[0].node_level = 0;

    TKET_ASSERT(!domain.empty());
    if (domain.size() == 1) {
      node.new_assignments.emplace_back(pv, *domain.cbegin());
    } else {
      node.unassigned_vertices_superset.push_back(pv);
    }
  }
}

const NodesRawData::NodeData& NodesRawData::get_current_node() const {
  return nodes_data[current_node_level];
}

NodesRawData::NodeData& NodesRawData::get_current_node_nonconst() {
  return nodes_data[current_node_level];
}

std::string NodesRawData::NodeData::str() const {
  std::stringstream ss;
  if (nogood) {
    ss << "##NOGOOD!## ";
  }
  ss << "Has " << new_assignments.size() << " ass.: [ ";
  for (const auto& entry : new_assignments) {
    ss << entry.first << ":" << entry.second << " ";
  }
  ss << "];  sc.prod " << scalar_product << "; p-edge weight "
     << total_p_edge_weights;
  return ss.str();
}

std::string NodesRawData::DomainData::str() const {
  std::stringstream ss;
  // ss << "\nDomains: (curr. node level=" << current_node_level << ")";
  for (unsigned ii = 0; ii <= entries_back_index; ++ii) {
    ss << "\n  lev=" << entries[ii].node_level << ", Dom: "
       << tket::WeightedSubgraphMonomorphism::str(entries[ii].domain);
  }
  ss << "\n";
  return ss.str();
}

NodesRawDataWrapper::NodesRawDataWrapper(
    const PossibleAssignments& possible_assignments)
    : m_raw_data(possible_assignments) {}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
