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

#include "WeightSubgrMono/GraphTheoretic/CloseNeighboursCalculator.hpp"

#include <algorithm>

#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

void CloseNeighboursCalculator::operator()(
    unsigned num_levels, CloseVertices& close_vertices,
    const NeighboursData& ndata, VertexWSM v) {
  close_vertices.resize(num_levels);
  if (close_vertices.empty()) {
    return;
  }
  close_vertices[0] = ndata.get_neighbours_expensive(v);
  m_vertices_seen.clear();
  m_vertices_seen.insert(v);
  for (VertexWSM old_v : close_vertices[0]) {
    m_vertices_seen.insert(old_v);
  }
  for (unsigned ii = 1; ii < close_vertices.size(); ++ii) {
    close_vertices[ii].clear();
    for (VertexWSM old_v : close_vertices[ii - 1]) {
      const auto& neighbours_and_weights =
          ndata.get_neighbours_and_weights(old_v);
      for (const auto& entry : neighbours_and_weights) {
        const VertexWSM& new_v = entry.first;
        if (m_vertices_seen.count(new_v) != 0) {
          continue;
        }
        m_vertices_seen.insert(new_v);
        close_vertices[ii].push_back(new_v);
      }
    }
    std::sort(close_vertices[ii].begin(), close_vertices[ii].end());
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
