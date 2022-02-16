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
#include <utility>
#include <vector>

#include "Utils/RNG.hpp"

namespace tket {
namespace graphs {

class AdjacencyData;

namespace tests {

/**
 * For having a whole sequence of checked edges
 * to add to a graph in a specific order,
 * and thus construct an increasing sequence of graphs to test
 * (not "subgraphs" because we are adding edges, not vertices).
 */
struct EdgeSequence {
  /** Will be used to check the edges for internal consistency. */
  AdjacencyData& adjacency_data;

  /** For convenience, have an RNG available to the caller. */
  RNG& rng;

  /**
   * Upon construction, for convenience store AdjacencyData and RNG objects.
   * The caller must ensure that the objects remain valid.
   */
  EdgeSequence(AdjacencyData& adj, RNG& rng);

  /** Erase all stored data in "edges" and the AdjacencyData object. */
  void clear();

  /**
   * Check if the edge already existed. If it didn't, add it to "edges".
   *
   * @param v1 The first vertex in the edge
   * @param v2 The second vertex in the edge
   * @return True if a new edge v1-v2 was added to "edges".
   * False if the edge already existed, so no action was taken.
   */
  bool add_edge(std::size_t v1, std::size_t v2);

  /** These are the added edges, in order. */
  std::vector<std::pair<std::size_t, std::size_t>> edges;
};

}  // namespace tests
}  // namespace graphs
}  // namespace tket
