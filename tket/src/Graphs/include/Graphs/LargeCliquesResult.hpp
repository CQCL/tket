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
#include <set>
#include <vector>

namespace tket {
namespace graphs {

class AdjacencyData;

/** Try to find large cliques in a single connected component of a graph. */
struct LargeCliquesResult {
  /**
   * These should all be valid cliques, of the same size.
   * However, might not be ALL the cliques.
   */
  std::vector<std::set<std::size_t>> cliques;

  /**
   * Are the cliques guaranteed to be of maximum possible size,
   * rather than just quite big?
   */
  bool cliques_are_definitely_max_size = false;

  /**
   * Given a connected component of a graph, try to find all the largest cliques
   * (unless the internal_size_limit is hit).
   *
   * @param adjacency_data The full graph.
   * @param vertices_in_component The vertices in the single connected component
   * we consider here. The caller is responsible for finding the components
   * correctly, although it's not checked.
   * @param internal_size_limit If this limit is hit, then the returned cliques
   * might not be of maximum size (although they should all be of the SAME
   * size), and we might be missing some of them. For colouring we don't need
   * the actual largest clique, so it may be acceptable. This limit might be hit
   * even if the actual final number of max size cliques is much smaller.
   */
  LargeCliquesResult(
      const AdjacencyData& adjacency_data,
      const std::set<std::size_t>& vertices_in_component,
      std::size_t internal_size_limit = 100);
};

}  // namespace graphs
}  // namespace tket
