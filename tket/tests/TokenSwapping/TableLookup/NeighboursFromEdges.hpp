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

#include <set>

#include "TokenSwapping/NeighboursInterface.hpp"
#include "TokenSwapping/SwapFunctions.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

/** Simply take a collection of swaps (or edges) and construct the neighbours
 * data. */
class NeighboursFromEdges : public NeighboursInterface {
 public:
  NeighboursFromEdges();

  template <class SwapContainer>
  explicit NeighboursFromEdges(const SwapContainer& edges);

  /** Add the edges one-by-one if desired.
   * @param edge An edge which you know is present in the graph.
   */
  void add_edge(const Swap& edge);

  /** The caller must not call this too soon, before "add_edge" calls are
   * completed.
   * @param vertex A vertex in the graph
   * @return All other vertices adjecent to the vertex (stored internally).
   */
  virtual const std::vector<size_t>& operator()(size_t vertex) override;

 private:
  /** The key is the vertex, the value is the list of neighbours. */
  std::map<size_t, std::set<size_t>> m_cached_neighbours;

  std::vector<size_t> m_neighbours_storage;
};

template <class SwapContainer>
NeighboursFromEdges::NeighboursFromEdges(const SwapContainer& edges) {
  for (const Swap& edge : edges) {
    add_edge(edge);
  }
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
