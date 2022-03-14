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

#include "NeighboursFromArchitecture.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace tket {

NeighboursFromArchitecture::NeighboursFromArchitecture(
    const ArchitectureMapping& arch_mapping)
    : m_arch_mapping(arch_mapping) {}

const std::vector<size_t>& NeighboursFromArchitecture::operator()(
    size_t vertex) {
  const auto num_vertices = m_arch_mapping.number_of_vertices();
  // GCOVR_EXCL_START
  TKET_ASSERT(
      vertex < num_vertices ||
      AssertMessage() << "get_neighbours: invalid vertex " << vertex
                      << " (only have " << num_vertices << " vertices)");
  // GCOVR_EXCL_STOP
  auto& neighbours = m_cached_neighbours[vertex];
  if (!neighbours.empty()) {
    // Already cached.
    return neighbours;
  }

  // OK, if a vertex is isolated (has no neighbours) then this is wasteful;
  // however this case should almost never occur in practice.

  const auto& source_node = m_arch_mapping.get_node(vertex);
  const auto neighbour_nodes =
      m_arch_mapping.get_architecture().get_neighbour_nodes(source_node);

  neighbours.reserve(neighbour_nodes.size());

  for (const Node& node : neighbour_nodes) {
    const auto neighbour_vertex = m_arch_mapping.get_vertex(node);
    // GCOVR_EXCL_START
    TKET_ASSERT(
        neighbour_vertex != vertex ||
        AssertMessage()
            << "get_neighbours: vertex " << vertex << " for node "
            << node.repr() << " has " << neighbour_nodes.size()
            << " neighbours, and lists itself as a neighbour (loops not "
               "allowed)");
    // GCOVR_EXCL_STOP
    neighbours.push_back(neighbour_vertex);
  }
  std::sort(neighbours.begin(), neighbours.end());
  return neighbours;
}

}  // namespace tket
