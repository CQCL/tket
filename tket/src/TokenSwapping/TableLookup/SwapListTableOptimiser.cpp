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

#include "TokenSwapping/SwapListTableOptimiser.hpp"

#include <algorithm>
#include <array>
#include <limits>

#include "Utils/Assert.hpp"

namespace tket {
namespace tsa_internal {

enum class EmptySwapCheckResult {
  NOT_EMPTY,
  CONTINUE_AFTER_ERASURE,
  TERMINATE_AFTER_ERASURE
};

// current_id is KNOWN to be valid.
// vertices_with_tokens is correct just BEFORE performing the swap.
// If the swap is empty, erase it and update current_id (to the next swap).
static EmptySwapCheckResult check_for_empty_swap(
    const std::set<size_t>& vertices_with_tokens, SwapID& current_id,
    SwapList& swap_list) {
  const auto swap = swap_list.at(current_id);
  if (vertices_with_tokens.count(swap.first) != 0 ||
      vertices_with_tokens.count(swap.second) != 0) {
    return EmptySwapCheckResult::NOT_EMPTY;
  }
  const auto next_id_opt = swap_list.next(current_id);
  swap_list.erase(current_id);
  if (!next_id_opt) {
    return EmptySwapCheckResult::TERMINATE_AFTER_ERASURE;
  }
  current_id = next_id_opt.value();
  return EmptySwapCheckResult::CONTINUE_AFTER_ERASURE;
}

// current_id is KNOWN to be valid.
// vertices_with_tokens is correct just BEFORE performing the swap.
// Keep erasing empty swaps and updating current_id
// until EITHER we hit a nonempty swap, OR we run out of swaps,
// and thus return false.
static bool erase_empty_swaps_interval(
    const std::set<size_t>& vertices_with_tokens, SwapID& current_id,
    SwapList& swap_list) {
  for (auto infinite_loop_guard = 1 + swap_list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    switch (check_for_empty_swap(vertices_with_tokens, current_id, swap_list)) {
      case EmptySwapCheckResult::CONTINUE_AFTER_ERASURE:
        // Maybe more to erase!
        break;
      case EmptySwapCheckResult::NOT_EMPTY:
        return true;
      case EmptySwapCheckResult::TERMINATE_AFTER_ERASURE:
        return false;
      default:
        TKET_ASSERT(!"unknown EmptySwapCheckResult enum");
        break;
    }
  }
  // Should never get here!
  TKET_ASSERT(!"erase_empty_swaps_interval failed to terminate");
  return false;
}

// current_id is KNOWN to be valid and nonempty.
// vertices_with_tokens is correct just BEFORE we perform the current swap.
// Perform the swap (i.e., updating vertices_with_tokens),
// and advance current_id to the next swap.
static bool perform_current_nonempty_swap(
    std::set<size_t>& vertices_with_tokens, SwapID& current_id,
    const SwapList& swap_list) {
  const auto swap = swap_list.at(current_id);

  if (vertices_with_tokens.count(swap.first) == 0) {
    // No empty swaps!
    TKET_ASSERT(vertices_with_tokens.count(swap.second) != 0);
    // Second has a token, first doesn't.
    TKET_ASSERT(vertices_with_tokens.insert(swap.first).second);
    TKET_ASSERT(vertices_with_tokens.erase(swap.second) == 1);
  } else {
    // First has a token.
    if (vertices_with_tokens.count(swap.second) == 0) {
      // Second has no token.
      TKET_ASSERT(vertices_with_tokens.erase(swap.first) == 1);
      TKET_ASSERT(vertices_with_tokens.insert(swap.second).second);
    }
  }

  const auto next_id_opt = swap_list.next(current_id);
  if (!next_id_opt) {
    return false;
  }
  current_id = next_id_opt.value();
  return true;
}

void SwapListTableOptimiser::optimise(
    const std::set<size_t>& vertices_with_tokens_at_start,
    VertexMapResizing& map_resizing, SwapList& swap_list,
    SwapListOptimiser& swap_list_optimiser) {
  if (vertices_with_tokens_at_start.empty()) {
    swap_list.clear();
    return;
  }
  if (swap_list.empty()) {
    return;
  }

  // Because we'll go in both directions, we need to know
  // which tokens exist at the END of the mapping.
  auto vertices_with_tokens_at_end = vertices_with_tokens_at_start;
  {
    // Already checked to be nonempty.
    auto current_id = swap_list.front_id().value();
    bool terminated_correctly = false;
    for (auto infinite_loop_guard = 1 + swap_list.size();
         infinite_loop_guard > 0; --infinite_loop_guard) {
      if (!erase_empty_swaps_interval(
              vertices_with_tokens_at_end, current_id, swap_list)) {
        terminated_correctly = true;
        break;
      }
      if (!perform_current_nonempty_swap(
              vertices_with_tokens_at_end, current_id, swap_list)) {
        terminated_correctly = true;
        break;
      }
    }
    TKET_ASSERT(terminated_correctly);
    if (swap_list.size() <= 1) {
      return;
    }
  }
  // Now begin the forward/backward loop.
  for (auto infinite_loop_guard = 1 + swap_list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    const auto old_size = swap_list.size();
    optimise_in_forward_direction(
        vertices_with_tokens_at_start, map_resizing, swap_list,
        swap_list_optimiser);

    swap_list.reverse();
    optimise_in_forward_direction(
        vertices_with_tokens_at_end, map_resizing, swap_list,
        swap_list_optimiser);

    // Must reverse again to get back to start!
    swap_list.reverse();
    const auto new_size = swap_list.size();
    TKET_ASSERT(new_size <= old_size);
    if (new_size == old_size) {
      return;
    }
  }
  TKET_ASSERT(!"SwapListTableOptimiser::optimise");
}

void SwapListTableOptimiser::optimise_in_forward_direction(
    const std::set<size_t>& vertices_with_tokens_at_start,
    VertexMapResizing& map_resizing, SwapList& swap_list,
    SwapListOptimiser& swap_list_optimiser) {
  swap_list_optimiser.optimise_pass_with_frontward_travel(swap_list);

  m_segment_optimiser.optimise_segment(
      swap_list.front_id().value(), vertices_with_tokens_at_start, map_resizing,
      swap_list);

  if (swap_list.size() <= 1) {
    return;
  }
  // Will always remain valid. We perform this swap and then optimise
  // starting from the next one.
  auto current_id = swap_list.front_id().value();
  auto vertices_with_tokens = vertices_with_tokens_at_start;

  for (size_t infinite_loop_guard = swap_list.size(); infinite_loop_guard != 0;
       --infinite_loop_guard) {
    if (!erase_empty_swaps_interval(
            vertices_with_tokens, current_id, swap_list)) {
      return;
    }
    // We now have a valid ID with nonempty swap.
    if (!perform_current_nonempty_swap(
            vertices_with_tokens, current_id, swap_list)) {
      return;
    }

    // NOW we want to optimise from this ID.
    // However, we must be careful; maybe it will be erased, so we have
    // to get the PREVIOUS and recover from there.
    const auto previous_id_opt = swap_list.previous(current_id);

    m_segment_optimiser.optimise_segment(
        current_id, vertices_with_tokens, map_resizing, swap_list);

    // We now want to set "current_id" to the first swap of
    // the newly optimised segment (if any) - which may of course
    // be unchanged, changed, or empty.

    // If there was no previous ID, we must have been at the front
    // just before we optimised.
    auto current_id_opt = swap_list.front_id();
    if (previous_id_opt) {
      // There WAS a previous ID, so we CAN move onto the next.
      current_id_opt = swap_list.next(previous_id_opt.value());
    }
    if (!current_id_opt) {
      return;
    }
    current_id = current_id_opt.value();
  }
}

SwapListSegmentOptimiser& SwapListTableOptimiser::get_segment_optimiser() {
  return m_segment_optimiser;
}

}  // namespace tsa_internal
}  // namespace tket
