#include "NeighboursFromArchitecture.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

;

namespace tket {
namespace tsa_internal {

NeighboursFromArchitecture::NeighboursFromArchitecture(
    const ArchitectureMapping& arch_mapping)
    : m_arch_mapping(arch_mapping) {}

const std::vector<size_t>& NeighboursFromArchitecture::operator()(
    size_t vertex) {
  const auto num_vertices = m_arch_mapping.number_of_vertices();
  if (vertex >= num_vertices) {
    std::stringstream ss;
    ss << "get_neighbours: invalid vertex " << vertex << " (only have "
       << num_vertices << " vertices)";
    throw std::runtime_error(ss.str());
  }
  auto& neighbours = m_cached_neighbours[vertex];
  if (!neighbours.empty()) {
    // Already cached.
    return neighbours;
  }

  // OK, if a vertex is isolated (has no neighbours) then this is wasteful;
  // however this case should almost never occur in practice.

  const auto& source_node = m_arch_mapping.get_node(vertex);
  const auto neighbour_nodes =
      m_arch_mapping.get_architecture().get_neighbour_uids(source_node);

  neighbours.reserve(neighbour_nodes.size());

  for (const Node& node : neighbour_nodes) {
    const auto neighbour_vertex = m_arch_mapping.get_vertex(node);
    if (neighbour_vertex == vertex) {
      std::stringstream ss;
      ss << "get_neighbours: vertex " << vertex << " for node " << node.repr()
         << " has " << neighbour_nodes.size()
         << " neighbours, and lists itself as a neighbour (loops not "
            "allowed)";
      throw std::runtime_error(ss.str());
    }
    neighbours.push_back(neighbour_vertex);
  }
  std::sort(neighbours.begin(), neighbours.end());
  return neighbours;
}

}  // namespace tsa_internal
}  // namespace tket
