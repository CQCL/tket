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

#include "DistancesFromArchitecture.hpp"

#include <sstream>
#include <stdexcept>

namespace tket {

DistancesFromArchitecture::DistancesFromArchitecture(
    const ArchitectureMapping& arch_mapping)
    : m_arch_mapping(arch_mapping) {}

void DistancesFromArchitecture::register_shortest_path(
    const std::vector<size_t>& path) {
  // To avoid quadratic growth for really long paths,
  // just do various slices.
  if (path.size() <= 5) {
    register_shortest_path_with_limits(path, 0, path.size());
    return;
  }
  const size_t middle = path.size() / 2;
  if (path.size() <= 10) {
    register_shortest_path_with_limits(path, 0, middle);
    register_shortest_path_with_limits(path, middle, path.size());
    register_edge(path[middle - 1], path[middle]);
    return;
  }
  register_shortest_path_with_limits(path, 0, 5);
  register_shortest_path_with_limits(path, path.size() - 5, path.size());
  if (path.size() >= 15) {
    register_shortest_path_with_limits(path, middle - 2, middle + 3);
  }
}

void DistancesFromArchitecture::register_shortest_path_with_limits(
    const std::vector<size_t>& path, size_t begin, size_t end) {
  for (size_t ii = begin; ii < end; ++ii) {
    for (size_t jj = ii + 1; jj < end; ++jj) {
      m_cached_distances[get_swap(path[ii], path[jj])] = jj - ii;
    }
  }
}

void DistancesFromArchitecture::register_edge(size_t vertex1, size_t vertex2) {
  m_cached_distances[get_swap(vertex1, vertex2)] = 1;
}

size_t DistancesFromArchitecture::operator()(size_t vertex1, size_t vertex2) {
  if (vertex1 == vertex2) {
    return 0;
  }
  // Automatically set to zero if it doesn't exist yet.
  auto& distance_entry = m_cached_distances[get_swap(vertex1, vertex2)];
  if (distance_entry == 0) {
    const auto& arch = m_arch_mapping.get_architecture();
    distance_entry = arch.get_distance(
        m_arch_mapping.get_node(vertex1), m_arch_mapping.get_node(vertex2));

    // This message should no longer be triggered for disconnected
    // architectures, since get_distance now should throw if v1, v2 are in
    // different connected components. However, leave the check in, in case some
    // other bizarre error causes distance zero to be returned.
    // GCOVR_EXCL_START
    TKET_ASSERT(
        distance_entry > 0 ||
        AssertMessage() << "DistancesFromArchitecture: architecture has "
                        << arch.n_nodes() << " vertices, "
                        << arch.n_connections() << " edges; "
                        << " and d(" << vertex1 << "," << vertex2
                        << ")=0. "
                           "Is the graph connected?");
    // GCOVR_EXCL_STOP
  }
  return distance_entry;
}

}  // namespace tket
