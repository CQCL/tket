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
#include "GeneralStructs.hpp"
#include <forward_list>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct DerivedGraphStructs {

  /** This is an actual count of paths in the original graph,
   * used as an edge weight in a derived graph.
   */
  typedef std::size_t Count;

  /** This is a list of neighbours in a dervied graph, where the "counts"
   * are edge weights (actual counts of vertices in the original graph).
   */
  typedef std::vector<std::pair<VertexWSM, Count>> NeighboursAndCounts;

  typedef std::forward_list<NeighboursAndCounts> List;
  typedef List::iterator Iter;
};


/** Put all the raw data for graphs in one place. */
class DerivedGraphsStorage {
public:

  /** The whole point is that we want references which are not invalidated
   * by lazy evaluation of other parts of the graphs.
   */
  DerivedGraphStructs::Iter get_new_neighbours_and_counts_iter();

private:
  
  // Stores the neighbours data, but without understanding the meaning.
  // (I.e., it's up to the caller to assign this data to appropriate graphs;
  // this class knows nothing about which graphs/vertices have which
  // neighbours and counts objects).
  DerivedGraphStructs::List m_stored_neighbours_and_counts;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
