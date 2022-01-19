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

#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <vector>

namespace tket {
namespace graphs {

class AdjacencyData;

/**
 * General easily reusable routines related to graphs. These are not
 * methods of the class AdjacencyData because they can be done efficiently
 * just with public methods of AdjacencyData.
 */
struct GraphRoutines {
  /**
   * Splits the graph into connected components.
   * @param adjacency_data The graph to split.
   * @return A collection of vertex sets. Each set gives all vertices in a
   * single component.
   */
  static std::vector<std::set<std::size_t>> get_connected_components(
      const AdjacencyData& adjacency_data);
};

}  // namespace graphs
}  // namespace tket
