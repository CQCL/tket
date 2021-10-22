#include "NeighboursFromEdges.hpp"

#include <algorithm>

;

namespace tket {
namespace tsa_internal {
namespace tests {

NeighboursFromEdges::NeighboursFromEdges() {}

void NeighboursFromEdges::add_edge(const Swap& edge) {
  m_cached_neighbours[edge.first].insert(edge.second);
  m_cached_neighbours[edge.second].insert(edge.first);
}

const std::vector<size_t>& NeighboursFromEdges::operator()(size_t vertex) {
  const auto& neighbours_set = m_cached_neighbours[vertex];
  m_neighbours_storage = {neighbours_set.cbegin(), neighbours_set.cend()};
  return m_neighbours_storage;
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
