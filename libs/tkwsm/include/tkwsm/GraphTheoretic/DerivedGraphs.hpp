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
#include "tkwsm/GraphTheoretic/DerivedGraphStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DerivedGraphsCalculator;
class NeighboursData;

/** The whole point is that the references remain valid,
 * taking care of dependencies on other data,
 * even though everything is lazily evaluated.
 */
class DerivedGraphs {
 public:
  /** The calculator is shared between several objects,
   * to reduce memory consumption.
   * @param ndata Data about the original graph we're deriving from.
   * @param calculator An object which knows how to calculate edges and weights
   * for derived graphs.
   */
  DerivedGraphs(
      const NeighboursData& ndata, DerivedGraphsCalculator& calculator);

  /** The point is, this is cheap to copy, AND the iters
   * provide references which remain valid even after new insertions.
   */
  struct VertexData {
    DerivedGraphStructs::Count triangle_count;

    // Of course one could consider more and more derived graphs,
    // and recursion: D2 of D2, D2 of D3 of D2, etc.
    // but quick tests seem to show that it isn't worthwhile
    // (takes more time to compute than the time saved).

    DerivedGraphStructs::NeighboursAndCountsIter d2_neighbours;
    DerivedGraphStructs::SortedCountsIter d2_sorted_counts_iter;

    DerivedGraphStructs::NeighboursAndCountsIter d3_neighbours;
    DerivedGraphStructs::SortedCountsIter d3_sorted_counts_iter;
  };

  /** Get data for a vertex. Notice that the data is cheap to copy,
   * and the pointers contained within it remain valid even after
   * new insertions.
   * @param v A vertex in the graph.
   * @return Data about that vertex, in the derived graphs.
   */
  VertexData get_data(VertexWSM v);

 private:
  const NeighboursData& m_neighbours_data;
  DerivedGraphsCalculator& m_calculator;
  DerivedGraphStructs::NeighboursAndCountsStorage m_storage;
  DerivedGraphStructs::SortedCountsStorage m_counts_storage;

  // KEY: a vertex
  // VALUE: data for that vertex; note that VertexData is cheap to copy.
  std::map<VertexWSM, VertexData> m_data_for_vertices;

  void fill(VertexWSM v, VertexData& vertex_data);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
