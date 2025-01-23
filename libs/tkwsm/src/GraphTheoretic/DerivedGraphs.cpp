// Copyright Quantinuum
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

#include "tkwsm/GraphTheoretic/DerivedGraphs.hpp"

#include "tkwsm/GraphTheoretic/DerivedGraphsCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphs::DerivedGraphs(
    const NeighboursData& ndata, DerivedGraphsCalculator& calculator)
    : m_neighbours_data(ndata), m_calculator(calculator) {}

DerivedGraphs::VertexData DerivedGraphs::get_data(VertexWSM v) {
  auto iter = m_data_for_vertices.find(v);
  if (iter != m_data_for_vertices.end()) {
    return iter->second;
  }
  auto& entry = m_data_for_vertices[v];
  fill(v, entry);
  return entry;
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

template <class T>
static typename std::forward_list<T>::iterator get_new_iter(
    std::forward_list<T>& storage) {
  storage.emplace_front();
  return storage.begin();
}

void DerivedGraphs::fill(VertexWSM v, VertexData& vertex_data) {
  vertex_data.d2_neighbours = get_new_iter(m_storage);
  vertex_data.d3_neighbours = get_new_iter(m_storage);
  m_calculator.fill(
      m_neighbours_data, v, vertex_data.triangle_count,
      *vertex_data.d2_neighbours, *vertex_data.d3_neighbours);

  vertex_data.d2_sorted_counts_iter = get_new_iter(m_counts_storage);
  vertex_data.d3_sorted_counts_iter = get_new_iter(m_counts_storage);

  fill_with_sorted_counts(
      *vertex_data.d2_sorted_counts_iter, *vertex_data.d2_neighbours);
  fill_with_sorted_counts(
      *vertex_data.d3_sorted_counts_iter, *vertex_data.d3_neighbours);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
