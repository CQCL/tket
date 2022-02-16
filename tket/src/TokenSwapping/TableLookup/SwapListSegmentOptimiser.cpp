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

#include "TokenSwapping/SwapListSegmentOptimiser.hpp"

#include <algorithm>
#include <array>
#include <limits>

#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

const SwapListSegmentOptimiser::Output&
SwapListSegmentOptimiser::optimise_segment(
    SwapID initial_id, const std::set<size_t>& vertices_with_tokens_at_start,
    VertexMapResizing& map_resizing, SwapList& swap_list) {
  m_best_optimised_swaps.clear();

  // Nonzero if and only if a valid sequence of swaps was stored.
  m_output.initial_segment_size = 0;

  // If the mapping has too many vertices, it MAY happen that
  // adding more swaps REDUCES the number of vertices
  // (since, some vertices may move back to their original positions,
  // and hence be "ignored"). Thus, we ALLOW the lookup to fail a few times
  // due to too many vertices before we give up.
  const int max_consecutive_too_many_vertices = 5;
  int too_many_vertices_count = max_consecutive_too_many_vertices;

  VertexMapping current_map;
  {
    const auto& initial_swap = swap_list.at(initial_id);
    current_map[initial_swap.first] = initial_swap.second;
    current_map[initial_swap.second] = initial_swap.first;
  }
  size_t current_number_of_swaps = 1;
  VertexMapping current_map_copy;
  for (auto next_id_opt = swap_list.next(initial_id);;) {
    bool too_many_vertices = false;

    // As we keep adding swaps to a sequence and updating the resultant
    // target->source vertex mapping, should we look up EVERY mapping in the
    // table, or is it enough to do so only when the map increases in size, etc.
    // etc.? Desperately need some theory here! We look up almost EVERYTHING, so
    // table lookup is one possible slowdown; reducing unnecessary lookups is
    // worthwhile.
    // TODO: think of some theory, and experiment!
    bool attempt_to_optimise = current_map.size() >= 3;
    if (!attempt_to_optimise && !next_id_opt) {
      // Because it's the FINAL segment, optimise it whatever we do.
      attempt_to_optimise = true;
    }
    if (attempt_to_optimise) {
      // We're going to attempt to optimise.
      current_map_copy = current_map;
      const auto& resize_result = map_resizing.resize_mapping(current_map);
      if (resize_result.success) {
        const auto& lookup_result = m_mapping_lookup(
            current_map, resize_result.edges, vertices_with_tokens_at_start,
            current_number_of_swaps);

        if (lookup_result.success) {
          // We've got a new result from the table; do we store it?
          bool should_store = m_output.initial_segment_size == 0;
          if (!should_store) {
            // Something IS stored, but is our new solution better?
            // GCOVR_EXCL_START
            TKET_ASSERT(
                m_output.initial_segment_size >= m_best_optimised_swaps.size());
            // GCOVR_EXCL_STOP
            const size_t current_decrease =
                m_output.initial_segment_size - m_best_optimised_swaps.size();
            TKET_ASSERT(current_number_of_swaps >= lookup_result.swaps.size());
            const size_t new_decrease =
                current_number_of_swaps - lookup_result.swaps.size();
            should_store = new_decrease > current_decrease;
          }
          if (should_store) {
            m_output.initial_segment_size = current_number_of_swaps;
            m_best_optimised_swaps = lookup_result.swaps;
          }
        } else {
          if (lookup_result.too_many_vertices) {
            too_many_vertices = true;
          }
        }
      } else {
        // We couldn't resize the mapping, so there must be too many vertices.
        too_many_vertices = true;
        // Also, the vertex mapping may be corrupted, so restore it
        current_map = current_map_copy;
      }
    }

    if (too_many_vertices) {
      --too_many_vertices_count;
      if (too_many_vertices_count <= 0) {
        break;
      }
    } else {
      too_many_vertices_count = max_consecutive_too_many_vertices;
    }

    // Now add a swap.
    if (next_id_opt) {
      const auto id = next_id_opt.value();
      const Swap swap = swap_list.at(id);
      add_swap(current_map, swap);
      ++current_number_of_swaps;
      next_id_opt = swap_list.next(id);
    } else {
      // We've reached the end!
      break;
    }
  }
  fill_final_output_and_swaplist(initial_id, swap_list);
  return m_output;
}

void SwapListSegmentOptimiser::fill_final_output_and_swaplist(
    SwapID initial_id, SwapList& swap_list) {
  if (m_output.initial_segment_size == 0) {
    // No improvement was found.
    m_output.final_segment_size = 0;
    m_output.new_segment_last_id = {};
    return;
  }
  m_output.final_segment_size = m_best_optimised_swaps.size();
  TKET_ASSERT(m_output.final_segment_size <= m_output.initial_segment_size);
  const auto initial_size = swap_list.size();

  if (m_best_optimised_swaps.empty()) {
    swap_list.erase_interval(initial_id, m_output.initial_segment_size);
    m_output.new_segment_last_id = {};
  } else {
    const auto overwrite_result = swap_list.overwrite_interval(
        initial_id, m_best_optimised_swaps.cbegin(),
        m_best_optimised_swaps.cend());

    // GCOVR_EXCL_START
    TKET_ASSERT(
        overwrite_result.number_of_overwritten_elements ==
        m_best_optimised_swaps.size());
    // GCOVR_EXCL_STOP
    m_output.new_segment_last_id =
        overwrite_result.final_overwritten_element_id;

    const size_t remaining_elements_to_erase =
        m_output.initial_segment_size - m_output.final_segment_size;

    const auto next_id_opt =
        swap_list.next(overwrite_result.final_overwritten_element_id);
    if (next_id_opt) {
      swap_list.erase_interval(
          next_id_opt.value(), remaining_elements_to_erase);
    }
  }
  // GCOVR_EXCL_START
  TKET_ASSERT(
      swap_list.size() + m_output.initial_segment_size ==
      initial_size + m_output.final_segment_size);
  // GCOVR_EXCL_STOP
}

}  // namespace tsa_internal
}  // namespace tket
