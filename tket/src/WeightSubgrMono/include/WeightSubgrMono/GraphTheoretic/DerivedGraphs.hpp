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
#include "DerivedGraphStructs.hpp"

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
  DerivedGraphs(
      const NeighboursData& ndata, DerivedGraphsCalculator& calculator,
      DerivedGraphStructs::NeighboursAndCountsStorage& storage,
      DerivedGraphStructs::SortedCountsStorage& counts_storage);

  /** The point is, this is cheap to copy, AND the iters
   * provide references which remain valid.
   */
  struct VertexData {
    DerivedGraphStructs::Count triangle_count;
    DerivedGraphStructs::NeighboursAndCountsStorage::Iter d2_neighbours;
    DerivedGraphStructs::SortedCountsStorage::Iter d2_sorted_counts_iter;

    DerivedGraphStructs::NeighboursAndCountsStorage::Iter d3_neighbours;
    DerivedGraphStructs::SortedCountsStorage::Iter d3_sorted_counts_iter;
  };

  VertexData get_data(VertexWSM v);

 private:
  const NeighboursData& m_neighbours_data;
  DerivedGraphsCalculator& m_calculator;
  DerivedGraphStructs::NeighboursAndCountsStorage& m_storage;
  DerivedGraphStructs::SortedCountsStorage& m_counts_storage;

  std::map<VertexWSM, VertexData> m_data;

  void fill(VertexWSM v, VertexData& vertex_data);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
