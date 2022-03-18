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
#include <memory>

#include "DerivedGraphsCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;
struct DerivedGraphs;

class DerivedGraphsUpdater {
 public:
  /** The updater needs to have a reference to derived graph objects
   * (so it knows what to update); but then, within this constructor,
   * the derived graph objects are given a pointer to this updater object
   * (so, each graph knows who to inform, when data needs to be updated).
   */
  DerivedGraphsUpdater(
      const NeighboursData& ndata, DerivedGraphsCalculator& calculator,
      DerivedGraphsStorage& storage);

  ~DerivedGraphsUpdater();

  void fill_data_in_container(VertexWSM v);

  DerivedGraphs& get_derived_graphs();

  unsigned get_n_verts() const;
  unsigned get_n_edges() const;

 private:
  const NeighboursData& m_original_graph;
  DerivedGraphsCalculator& m_calculator;
  DerivedGraphsStorage& m_storage;
  std::unique_ptr<DerivedGraphs> m_derived_graphs_ptr;
};

struct DerivedGraphsUpdaterPair {
  DerivedGraphsUpdater patterns_updater;
  DerivedGraphsUpdater targets_updater;

  DerivedGraphsUpdaterPair(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
      DerivedGraphsCalculator& calculator, DerivedGraphsStorage& storage);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
