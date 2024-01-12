// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <forward_list>

#include "GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct DerivedGraphStructs {
  /** This is an actual count of paths in the original graph,
   * used as an edge weight in a derived graph.
   */
  typedef std::size_t Count;

  /** This is a list of neighbours in a derived graph, where the "counts"
   * are edge weights (actual counts of vertices in the original graph).
   * Sorted by vertex number, for easy neighbours lookup.
   */
  typedef std::vector<std::pair<VertexWSM, Count>> NeighboursAndCounts;

  /** We create derived graphs lazily (since many vertices and edges
   * might never be needed).
   * For convenience, we want to pass around references or pointers,
   * but we can't just store the objects in a vector
   * because the references might be invalidated upon insertions.
   * Thus we use a linked list.
   */
  typedef std::forward_list<NeighboursAndCounts> NeighboursAndCountsStorage;
  typedef NeighboursAndCountsStorage::iterator NeighboursAndCountsIter;

  /** These are the counts alone, which are in "NeighboursAndCounts",
   * sorted (for fast compatibility checking). */
  typedef std::vector<Count> SortedCounts;
  typedef std::forward_list<SortedCounts> SortedCountsStorage;
  typedef SortedCountsStorage::iterator SortedCountsIter;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
