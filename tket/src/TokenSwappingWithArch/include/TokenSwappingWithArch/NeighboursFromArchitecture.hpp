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

#include "ArchitectureMapping.hpp"
#include "TokenSwapping/NeighboursInterface.hpp"

namespace tket {

/** Stores and returns upon request the adjacent vertices to a given vertex
 *  on a graph, using an underlying Architecture object.
 */
class NeighboursFromArchitecture : public NeighboursInterface {
 public:
  /** The objects must remain valid AND unchanged
   *  for the lifetime of this object.
   *  @param arch_mapping An object which contains a reference to an
   *    Architecture object internally, and handles Node -> vertex size_t
   * conversions.
   */
  explicit NeighboursFromArchitecture(const ArchitectureMapping& arch_mapping);

  /** For extra convenience, the list of neighbours is always sorted
   *  in increasing order (so you can do binary search, etc.)
   *  @param vertex A vertex.
   *  @return A sorted list of all adjacent vertices, stored internally.
   */
  virtual const std::vector<size_t>& operator()(size_t vertex) override;

 private:
  const ArchitectureMapping& m_arch_mapping;

  /** The key is the vertex, the value is the list of neighbours. */
  std::map<size_t, std::vector<size_t>> m_cached_neighbours;
};

}  // namespace tket
