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

#include "ArchitectureMapping.hpp"

#include <sstream>
#include <stdexcept>

#include "Utils/Assert.hpp"

namespace tket {

ArchitectureMapping::ArchitectureMapping(const Architecture& arch)
    : m_arch(arch) {
  const auto uids = arch.nodes();
  m_vertex_to_node_mapping.reserve(uids.size());
  for (const UnitID& uid : uids) {
    m_vertex_to_node_mapping.emplace_back(Node(uid));
  }

  for (size_t ii = 0; ii < m_vertex_to_node_mapping.size(); ++ii) {
    const auto& node = m_vertex_to_node_mapping[ii];
    {
      const auto citer = m_node_to_vertex_mapping.find(node);
      // GCOVR_EXCL_START
      TKET_ASSERT(
          citer == m_node_to_vertex_mapping.cend() ||
          AssertMessage() << "Duplicate node " << node.repr() << " at vertices "
                          << citer->second << ", " << ii);
      // GCOVR_EXCL_STOP
    }
    m_node_to_vertex_mapping[node] = ii;
  }
}

ArchitectureMapping::ArchitectureMapping(
    const Architecture& arch,
    const std::vector<std::pair<unsigned, unsigned>>& edges)
    : m_arch(arch) {
  auto& node_to_vertex_mapping = m_node_to_vertex_mapping;
  auto& vertex_to_node_mapping = m_vertex_to_node_mapping;

  const auto add_node = [&node_to_vertex_mapping,
                         &vertex_to_node_mapping](unsigned nn) {
    const Node node(nn);
    if (node_to_vertex_mapping.count(node) == 0) {
      node_to_vertex_mapping[node] = vertex_to_node_mapping.size();
      vertex_to_node_mapping.push_back(node);
    }
  };

  // The nodes are labelled 0,1,2,... in order of appearance.
  // Nothing special about this ordering, just for backwards compatibility.
  for (const auto& entry : edges) {
    add_node(entry.first);
    add_node(entry.second);
  }

  // Check that the nodes agree with the architecture object.
  const auto uids = arch.nodes();
  // GCOVR_EXCL_START
  TKET_ASSERT(
      uids.size() == m_vertex_to_node_mapping.size() ||
      AssertMessage() << "passed in " << edges.size() << " edges, giving "
                      << m_vertex_to_node_mapping.size()
                      << " vertices; but the architecture object has "
                      << uids.size() << " vertices");
  // GCOVR_EXCL_STOP

  for (const UnitID& uid : uids) {
    const Node node(uid);
    // GCOVR_EXCL_START
    TKET_ASSERT(
        m_node_to_vertex_mapping.count(node) != 0 ||
        AssertMessage()
            << "passed in " << edges.size() << " edges, giving "
            << m_vertex_to_node_mapping.size()
            << " vertices; but the architecture object has an unknown node "
            << node.repr());
    // GCOVR_EXCL_STOP
  }
}

size_t ArchitectureMapping::number_of_vertices() const {
  return m_vertex_to_node_mapping.size();
}

const Node& ArchitectureMapping::get_node(size_t vertex) const {
  const auto num_vertices = number_of_vertices();
  // GCOVR_EXCL_START
  TKET_ASSERT(
      vertex < num_vertices || AssertMessage() << "invalid vertex " << vertex
                                               << " (architecture only has "
                                               << num_vertices << " vertices)");
  // GCOVR_EXCL_STOP

  return m_vertex_to_node_mapping[vertex];
}

size_t ArchitectureMapping::get_vertex(const Node& node) const {
  const auto citer = m_node_to_vertex_mapping.find(node);
  // GCOVR_EXCL_START
  TKET_ASSERT(
      citer != m_node_to_vertex_mapping.cend() ||
      AssertMessage() << "node " << node.repr() << " has no vertex number");
  // GCOVR_EXCL_STOP
  return citer->second;
}

const Architecture& ArchitectureMapping::get_architecture() const {
  return m_arch;
}

std::vector<Swap> ArchitectureMapping::get_edges() const {
  std::vector<Swap> edges;
  for (auto [node1, node2] : m_arch.get_all_edges_vec()) {
    edges.emplace_back(get_swap(get_vertex(node1), get_vertex(node2)));
  }
  return edges;
}

}  // namespace tket
