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

class DerivedGraphsUpdater;

class DerivedGraph {
public:

  explicit DerivedGraph(DerivedGraphsUpdater& updater);

  /** If the data does not already exist, will evaluate it lazily.
   * Crucially, the reference remains valid (as long as the
   * internal DerivedGraphsUpdater object remains valid),
   * even though other data may be created.
   */
  const DerivedGraphStructs::NeighboursAndCounts& get_neighbours(VertexWSM v);

private:
  DerivedGraphsUpdater& m_updater;
  std::map<VertexWSM, DerivedGraphStructs::Iter> m_data;

  /** Requires the entry for v NOT to exist before this call. */
  void add_neighbours(VertexWSM v, DerivedGraphStructs::Iter iter);
  
  friend class DerivedGraphsUpdater;
};


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
