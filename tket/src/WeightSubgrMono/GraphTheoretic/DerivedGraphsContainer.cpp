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

#include "WeightSubgrMono/GraphTheoretic/DerivedGraphsContainer.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphsContainer::VertexData& DerivedGraphsContainer::get_pattern_v_data_permanent_reference(VertexWSM v, const NeighboursData& pattern_ndata) {
  return get_vertex_data_permanent_reference(v, pattern_ndata, m_pattern_iters);
}

DerivedGraphsContainer::VertexData& DerivedGraphsContainer::get_target_v_data_permanent_reference(VertexWSM v, const NeighboursData& target_ndata) {
  return get_vertex_data_permanent_reference(v, target_ndata, m_target_iters);
}

DerivedGraphsContainer::VertexData& DerivedGraphsContainer::get_vertex_data_permanent_reference(
          VertexWSM v, const NeighboursData& ndata, std::map<VertexWSM, Iter>& map) {
  auto map_iter = map.find(v);
  if(map_iter != map.end()) {
    return *map_iter->second;
  }
  auto& element = map[v];
  // Make and store a new vertex data object.
  m_list_container.emplace_front();
  element = m_list_container.begin();

  // Initialise it!
  element->triangle_count = m_calculator.fill_neighbours_and_weights(ndata, v,
      element->depth_2_neighbours,
      element->depth_3_neighbours);
  return *element;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
