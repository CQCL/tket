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

#include "NeighboursFromEdges.hpp"

#include <algorithm>

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
