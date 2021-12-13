#include "ArchitectureMapping.hpp"

#include <sstream>
#include <stdexcept>

namespace tket {
namespace tsa_internal {

ArchitectureMapping::ArchitectureMapping(const Architecture& arch)
    : m_arch(arch) {
  const auto uids = arch.get_all_nodes();
  m_vertex_to_node_mapping.reserve(uids.size());
  for (const UnitID& uid : uids) {
    m_vertex_to_node_mapping.emplace_back(Node(uid));
  }

  for (size_t ii = 0; ii < m_vertex_to_node_mapping.size(); ++ii) {
    const auto& node = m_vertex_to_node_mapping[ii];
    {
      const auto citer = m_node_to_vertex_mapping.find(node);
      if (citer != m_node_to_vertex_mapping.cend()) {
        std::stringstream ss;
        ss << "Duplicate node " << node.repr() << " at vertices "
           << citer->second << ", " << ii;
        throw std::runtime_error(ss.str());
      }
    }
    m_node_to_vertex_mapping[node] = ii;
  }
}

size_t ArchitectureMapping::number_of_vertices() const {
  return m_vertex_to_node_mapping.size();
}

const Node& ArchitectureMapping::get_node(size_t vertex) const {
  const auto num_vertices = number_of_vertices();
  if (vertex >= num_vertices) {
    std::stringstream ss;
    ss << "get_node: invalid vertex " << vertex << " (architecture only has "
       << num_vertices << " vertices)";
    throw std::runtime_error(ss.str());
  }
  return m_vertex_to_node_mapping[vertex];
}

size_t ArchitectureMapping::get_vertex(const Node& node) const {
  const auto citer = m_node_to_vertex_mapping.find(node);
  if (citer == m_node_to_vertex_mapping.cend()) {
    std::stringstream ss;
    ss << "get_vertex: node " << node.repr() << " has no vertex number";
    throw std::runtime_error(ss.str());
  }
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

}  // namespace tsa_internal
}  // namespace tket
