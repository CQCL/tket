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

#include "TokenSwapping/BestFullTsa.hpp"

#include "TokenSwapping/RiverFlowPathFinder.hpp"
#include "TokenSwapping/VertexMapResizing.hpp"

namespace tket {

using namespace tsa_internal;

BestFullTsa::BestFullTsa() { m_name = "BestFullTsa"; }

void BestFullTsa::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, NeighboursInterface& neighbours,
    RiverFlowPathFinder& path_finder) {
  auto vm_copy = vertex_mapping;

  m_hybrid_tsa.append_partial_solution(
      swaps, vm_copy, distances, neighbours, path_finder);

  // Still subject to experimentation, but this seems the best
  m_swap_list_optimiser.optimise_pass_with_zero_travel(swaps);
  m_swap_list_optimiser.optimise_pass_with_token_tracking(swaps);
  m_swap_list_optimiser.optimise_pass_remove_empty_swaps(swaps, vertex_mapping);
  m_swap_list_optimiser.full_optimise(swaps, vertex_mapping);

  VertexMapResizing map_resizing(neighbours);
  std::set<size_t> vertices_with_tokens_at_start;
  for (const auto& entry : vertex_mapping) {
    vertices_with_tokens_at_start.insert(entry.first);
  }
  m_table_optimiser.optimise(
      vertices_with_tokens_at_start, map_resizing, swaps,
      m_swap_list_optimiser);
}

}  // namespace tket
