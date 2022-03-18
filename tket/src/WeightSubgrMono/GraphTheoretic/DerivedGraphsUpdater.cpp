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

#include "WeightSubgrMono/GraphTheoretic/DerivedGraph.hpp"
#include "WeightSubgrMono/GraphTheoretic/DerivedGraphs.hpp"
#include "WeightSubgrMono/GraphTheoretic/DerivedGraphsUpdater.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "Utils/Assert.hpp"
#include <utility>

namespace tket {
namespace WeightedSubgraphMonomorphism {

unsigned DerivedGraphsUpdater::get_n_verts() const {
  return m_original_graph.get_number_of_nonisolated_vertices();
}

unsigned DerivedGraphsUpdater::get_n_edges() const {
  return m_original_graph.get_number_of_edges();
}


DerivedGraphsUpdater::DerivedGraphsUpdater(const NeighboursData& ndata, DerivedGraphsCalculator& calculator,
        DerivedGraphsStorage& storage)
    : m_original_graph(ndata), m_calculator(calculator), m_storage(storage) {

  m_derived_graphs_ptr = std::make_unique<DerivedGraphs>(*this);
  TKET_ASSERT(m_derived_graphs_ptr);
}

DerivedGraphsUpdater::~DerivedGraphsUpdater() {}


DerivedGraphs& DerivedGraphsUpdater::get_derived_graphs() {
  return *m_derived_graphs_ptr;
}

void DerivedGraphsUpdater::fill_data_in_container(VertexWSM v) {
  auto d2_iter = m_storage.get_new_neighbours_and_counts_iter();
  auto d3_iter = m_storage.get_new_neighbours_and_counts_iter();

  m_derived_graphs_ptr->triangle_counts.fill_count(v,
    m_calculator.fill_neighbours_and_weights(
        m_original_graph, v, *d2_iter, *d3_iter));
  
  m_derived_graphs_ptr->d2_graph.add_neighbours(v, d2_iter);
  m_derived_graphs_ptr->d3_graph.add_neighbours(v, d3_iter);
}


DerivedGraphsUpdaterPair::DerivedGraphsUpdaterPair(
          const NeighboursData& pattern_ndata,
          const NeighboursData& target_ndata,
          DerivedGraphsCalculator& calculator,
          DerivedGraphsStorage& storage)
    : patterns_updater(pattern_ndata, calculator, storage),
      targets_updater(target_ndata, calculator, storage)
{}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
