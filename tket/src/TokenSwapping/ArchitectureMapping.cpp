#include "ArchitectureMapping.hpp"

#include <sstream>
#include <stdexcept>

namespace tket {
namespace tsa_internal {

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
  if (uids.size() != m_vertex_to_node_mapping.size()) {
    std::stringstream ss;
    ss << "ArchitectureMapping: passed in " << edges.size() << " edges, giving "
       << m_vertex_to_node_mapping.size()
       << " vertices; but the architecture object has " << uids.size()
       << " vertices";
    throw std::runtime_error(ss.str());
  }
  for (const UnitID& uid : uids) {
    const Node node(uid);
    if (m_node_to_vertex_mapping.count(node) == 0) {
      std::stringstream ss;
      ss << "ArchitectureMapping: passed in " << edges.size()
         << " edges, giving " << m_vertex_to_node_mapping.size()
         << " vertices; but the architecture object has an unknown node "
         << node.repr();
      throw std::runtime_error(ss.str());
    }
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
