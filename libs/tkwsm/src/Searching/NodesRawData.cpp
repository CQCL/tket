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

#include "tkwsm/Searching/NodesRawData.hpp"

#include <sstream>
#include <stdexcept>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Common/TemporaryRefactorCode.hpp"

#include "WeightSubgrMono/Common/TemporaryRefactorCode.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {


NodesRawData::NodesRawData(
    const DomainInitialiser::InitialDomains& initial_domains,
    std::size_t num_tv) 
    : number_of_tv(num_tv) 
      {
  nodes_data.push();
  auto& node = nodes_data.top();
  node.nogood = false;
  node.scalar_product = 0;
  node.total_p_edge_weights = 0;

  TKET_ASSERT(initial_domains.size() >= 2);
  domains_data.resize(initial_domains.size());

  for (unsigned pv = 0; pv < initial_domains.size(); ++pv) {
    const std::set<VertexWSM>& domain_set = initial_domains[pv];
    auto& domain_data = domains_data[pv];
    domain_data.entries.push();
    domain_data.entries[0].node_index = 0;

    TemporaryRefactorCode::set_bitset(domain_set, domain_data.entries[0].domain, number_of_tv);
    
    switch (domain_set.size()) {
      case 0: {
        std::stringstream ss;
        ss << "NodesRawData: Domain(" << pv << ") is empty!";
        throw std::runtime_error(ss.str());
      }
      case 1:
        node.new_assignments.emplace_back(pv, *domain_set.cbegin());
        break;
      default:
        node.unassigned_vertices_superset.push_back(pv);
    }
    TKET_ASSERT(*domain_set.crbegin() < number_of_tv);
  }
}

unsigned NodesRawData::current_node_index() const {
  return nodes_data.size() - 1;
}

const NodesRawData::NodeData& NodesRawData::get_current_node() const {
  return nodes_data.top();
}

NodesRawData::NodeData& NodesRawData::get_current_node_nonconst() {
  return nodes_data.top();
}

std::string NodesRawData::NodeData::str() const {
  std::stringstream ss;
  if (nogood) {
    ss << "##NOGOOD!## ";
  }
  ss << "Has " << new_assignments.size() << " ass.: [ ";
  for (const std::pair<VertexWSM, VertexWSM>& entry : new_assignments) {
    ss << entry.first << ":" << entry.second << " ";
  }
  ss << "];  sc.prod " << scalar_product << "; p-edge weight "
     << total_p_edge_weights;
  return ss.str();
}

std::string NodesRawData::DomainData::str() const {
  std::stringstream ss;
  const unsigned size = entries.size();

  std::set<VertexWSM> dom_temp;

  for (unsigned ii = 0; ii < size; ++ii) {
    TemporaryRefactorCode::set_domain_from_bitset(dom_temp, entries[ii].domain);
    ss << "\n  node_index=" << entries[ii].node_index << ", Dom: "
       << tket::WeightedSubgraphMonomorphism::str(dom_temp);
  }
  ss << "\n";
  return ss.str();
}

NodesRawDataWrapper::NodesRawDataWrapper(
    const DomainInitialiser::InitialDomains& initial_domains,
    std::size_t number_of_tv)
    : m_raw_data(initial_domains, number_of_tv) {}

const NodesRawData& NodesRawDataWrapper::get_raw_data_for_debug() const {
  return m_raw_data;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
