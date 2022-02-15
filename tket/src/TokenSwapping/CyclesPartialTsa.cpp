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

#include "CyclesPartialTsa.hpp"

#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

CyclesPartialTsa::CyclesPartialTsa() { m_name = "Cycles"; }

void CyclesPartialTsa::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, NeighboursInterface& neighbours,
    RiverFlowPathFinder& path_finder) {
  // We'll add the calculated swaps to the path finder at the end.
  // THIS is the right place to do it, not the caller, because
  // (as far as the caller knows) it's possible that PartialTSA objects
  // reduce/reorder swaps, and so it would be invalid just to go back through
  // the appended swaps. However, THIS class knows that no reordering or
  // reduction occurs.
  const size_t initial_swap_size = swaps.size();
  for (;;) {
    const auto swap_size_before = swaps.size();
    single_iteration_partial_solution(
        swaps, vertex_mapping, distances, neighbours);
    const auto swap_size_after = swaps.size();
    TKET_ASSERT(swap_size_after >= swap_size_before);
    if (swap_size_before == swap_size_after) {
      break;
    }
  }
  const size_t final_swap_size = swaps.size();
  TKET_ASSERT(initial_swap_size <= final_swap_size);
  if (initial_swap_size == final_swap_size) {
    return;
  }
  // At least one swap was added.
  const auto current_back_id_opt = swaps.back_id();
  TKET_ASSERT(current_back_id_opt);
  auto current_id = current_back_id_opt.value();
  for (size_t remaining_swaps = final_swap_size - initial_swap_size;;) {
    const auto& swap = swaps.at(current_id);
    path_finder.register_edge(swap.first, swap.second);
    --remaining_swaps;
    if (remaining_swaps == 0) {
      break;
    }
    const auto prev_id_opt = swaps.previous(current_id);
    TKET_ASSERT(prev_id_opt);
    current_id = prev_id_opt.value();
  }
}

void CyclesPartialTsa::single_iteration_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, NeighboursInterface& neighbours) {
  if (!m_growth_manager.reset(vertex_mapping, distances, neighbours)) {
    // no solutions.
    return;
  }

  for (auto infinite_loop_guard = m_growth_manager.get_options().max_cycle_size;
       infinite_loop_guard > 0; --infinite_loop_guard) {
    if (m_growth_manager.attempt_to_close_cycles(vertex_mapping, distances)) {
      // Some solutions found.
      m_candidate_manager.append_partial_solution(
          m_growth_manager, swaps, vertex_mapping);
      return;
    }
    // No solutions so far, so grow...
    const auto growth_result =
        m_growth_manager.attempt_to_grow(vertex_mapping, distances, neighbours);
    if (growth_result.empty || growth_result.hit_cycle_length_limit) {
      return;
    }
  }
  TKET_ASSERT(!"growth_manager termination");
}

}  // namespace tsa_internal
}  // namespace tket
