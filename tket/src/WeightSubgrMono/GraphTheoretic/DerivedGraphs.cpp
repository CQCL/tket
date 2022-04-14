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

#include "WeightSubgrMono/GraphTheoretic/DerivedGraphs.hpp"

#include "WeightSubgrMono/GraphTheoretic/DerivedGraphsCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphs::DerivedGraphs(
    const NeighboursData& ndata, DerivedGraphsCalculator& calculator,
    DerivedGraphStructs::NeighboursAndCountsStorage& storage,
    DerivedGraphStructs::SortedCountsStorage& counts_storage)
    : m_neighbours_data(ndata),
      m_calculator(calculator),
      m_storage(storage),
      m_counts_storage(counts_storage) {}

DerivedGraphs::VertexData DerivedGraphs::get_data(VertexWSM v) {
  auto iter = m_data.find(v);
  if (iter == m_data.end()) {
    auto& entry = m_data[v];
    fill(v, entry);
    return entry;
  }
  return iter->second;
}

static void fill_with_sorted_counts(
    DerivedGraphStructs::SortedCounts& counts,
    const DerivedGraphStructs::NeighboursAndCounts& neighbours_and_counts) {
  counts.reserve(neighbours_and_counts.size());
  for (const auto& entry : neighbours_and_counts) {
    counts.push_back(entry.second);
  }
  std::sort(counts.begin(), counts.end());
}

void DerivedGraphs::fill(VertexWSM v, VertexData& vertex_data) {
  vertex_data.d2_neighbours = m_storage.get_new_iter();
  vertex_data.d3_neighbours = m_storage.get_new_iter();
  m_calculator.fill(
      m_neighbours_data, v, vertex_data.triangle_count,
      *vertex_data.d2_neighbours, *vertex_data.d3_neighbours);

  vertex_data.d2_sorted_counts_iter = m_counts_storage.get_new_iter();
  vertex_data.d3_sorted_counts_iter = m_counts_storage.get_new_iter();

  fill_with_sorted_counts(
      *vertex_data.d2_sorted_counts_iter, *vertex_data.d2_neighbours);
  fill_with_sorted_counts(
      *vertex_data.d3_sorted_counts_iter, *vertex_data.d3_neighbours);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
