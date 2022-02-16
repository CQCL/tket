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

#include "SwapListOptimiser.hpp"

#include "Utils/Assert.hpp"
#include "VertexSwapResult.hpp"

namespace tket {
namespace tsa_internal {

void SwapListOptimiser::push_back(SwapList& list, const Swap& swap) {
  if (list.empty() || list.back() != swap) {
    list.push_back(swap);
    return;
  }
  list.pop_back();
}

// It may be that using a std::set is very slightly quicker
// (to store only the vertices containing tokens, as we don't care about the
// targets). However, it's simpler just to use the copied VertexMapping; not
// worth worrying about. (Also, if a std::map is large, then copying all keys
// into a std::set might actually be SLOWER than copying the whole map;
// HOPEFULLY the compiler can copy a whole map very quickly just by copying raw
// bytes, but for a std::set it would have to insert the keys one-by-one and do
// a lot of tree rebalancing).
void SwapListOptimiser::optimise_pass_remove_empty_swaps(
    SwapList& list, VertexMapping vertex_mapping) {
  auto id_opt = list.front_id();
  while (id_opt) {
    const auto id = id_opt.value();
    id_opt = list.next(id);
    const VertexSwapResult result(list.at(id), vertex_mapping);
    if (result.tokens_moved == 0) {
      list.erase(id);
    }
  }
}

std::optional<SwapListOptimiser::ID>
SwapListOptimiser::get_id_of_previous_blocker(SwapList& list, SwapID id) {
  const auto& initial_swap = list.at(id);

  // This is the first non-disjoint swap it hits when it moves back.
  // Guaranteed to be valid if we drop out of the loop.
  SwapID current_id = id;

  bool terminated_correctly = false;
  for (auto infinite_loop_guard = 1 + list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    const auto prev_id = list.previous(current_id);
    if (!prev_id) {
      // Right at the front!
      return {};
    }
    current_id = prev_id.value();
    const auto& new_swap = list.at(current_id);
    if (!disjoint(initial_swap, new_swap)) {
      // Blocks, OR identical
      if (new_swap != initial_swap) {
        // It blocks
        return current_id;
      }
      terminated_correctly = true;
      break;
    }
  }
  TKET_ASSERT(terminated_correctly);
  // It's hit a copy of itself
  list.erase(id);
  list.erase(current_id);
  return {};
}

bool SwapListOptimiser::move_swap_towards_front(SwapList& list, SwapID id) {
  TKET_ASSERT(list.front_id());
  if (id == list.front_id().value()) {
    return false;
  }
  const auto old_size = list.size();
  const auto previous_blocker_opt = get_id_of_previous_blocker(list, id);
  if (old_size != list.size()) {
    // The swap was erased!
    return true;
  }
  if (previous_blocker_opt) {
    // It can't move all the way to the front.
    const ID blocker = previous_blocker_opt.value();

    // Must be non-null.
    const ID previous_id = list.previous(id).value();
    if (blocker != previous_id) {
      // Do the move...erase before insert to minimise possible sizes...
      const auto swap = list.at(id);
      list.erase(id);
      const auto new_id = list.insert_after(blocker);
      list.at(new_id) = swap;
    }
    return false;
  }
  // There was no blocker, so we CAN move all the way to the front
  // (and we checked before that we're not already at the front).
  const auto swap = list.at(id);
  list.erase(id);
  list.push_front(swap);
  return false;
}

void SwapListOptimiser::optimise_pass_with_zero_travel(SwapList& list) {
  if (list.size() <= 1) {
    return;
  }
  ID current_id = list.front_id().value();

  // This moves swaps to cancel with previous identical swaps,
  // if there is nothing blocking the move.
  // However, only worth doing if previous identical swaps do exist;
  // leaves them unchanged otherwise.
  // We can be sneaky: rather than storing all previous IDs
  // for each swap, we store the NUMBER of them; we don't need
  // to know the previous location, since the move back
  // will check for that anyway.
  //
  // This probably could be cleverly optimised further
  // but would require more thought.
  for (auto& entry : m_data) {
    // This is quicker than clearing and reinserting;
    // no tree rebalancing.
    entry.second = 0;
  }
  for (auto infinite_loop_guard = 1 + list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    const auto next_id_opt = list.next(current_id);

    // C++ guarantees nonexistent values will be set to 0.
    auto& swap_count = m_data[list.at(current_id)];
    if (swap_count == 0) {
      swap_count = 1;
    } else {
      // There's a possibility of cancellation.
      const auto old_size = list.size();
      (void)get_id_of_previous_blocker(list, current_id);
      if (old_size == list.size()) {
        // No cancellation.
        ++swap_count;
      } else {
        // Cancellation occurred; "get_id_of_previous_blocker" already erased
        // both vertex swaps, but didn't update the counts.
        --swap_count;
      }
    }
    if (!next_id_opt) {
      return;
    }
    current_id = next_id_opt.value();
  }
  TKET_ASSERT(!"optimise_pass_with_zero_travel termination");
}

void SwapListOptimiser::optimise_pass_with_frontward_travel(SwapList& list) {
  if (list.size() <= 1) {
    return;
  }
  // Start one past the front.
  ID current_id = list.front_id().value();
  current_id = list.next(current_id).value();

  for (auto infinite_loop_guard = 1 + list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    const auto next_id_opt = list.next(current_id);
    move_swap_towards_front(list, current_id);
    if (!next_id_opt) {
      return;
    }
    current_id = next_id_opt.value();
  }
  TKET_ASSERT(!"optimise_pass_with_frontward_travel termination");
}

void SwapListOptimiser::optimise_pass_with_token_tracking(SwapList& list) {
  if (list.size() <= 1) {
    return;
  }
  m_token_tracker.clear();
  optimise_pass_with_token_tracking_without_clearing_tracker(list);
}

void SwapListOptimiser::
    optimise_pass_with_token_tracking_without_clearing_tracker(SwapList& list) {
  if (list.size() <= 1) {
    return;
  }
  // Put a different token at each vertex, and start swapping.
  // Now, if a TOKEN swap (rather than vertex swap)
  // repeats, then removing this vertex swap together with the preceding one
  // in which those two tokens were exchanged gives the same final result.
  // This is because, if we don't actually carry out the first swap,
  // everything proceeds as before, and all tokens except those two
  // are in the same place. When we reach the time of the second swap,
  // everything is as before EXCEPT that the two tokens have changed places;
  // thus the effect of the second swap
  // was merely to interchange those two tokens again.
  //
  // Now, m_data will store the previous LOCATIONS of vertex swaps.
  //
  // The actual values of the tokens are irrelevant,
  // as long as they are distinct.

  const auto invalid_index = VectorListHybridSkeleton::get_invalid_index();

  for (auto infinite_loop_guard = 1 + list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    // Keep looping until we stop changing.
    // The size is always decreasing or unchanged;
    // we never insert, only erase.
    const auto old_size = list.size();
    if (old_size == 0) {
      return;
    }
    for (auto& entry : m_data) {
      entry.second = invalid_index;
    }
    ID current_id = list.front_id().value();
    bool terminated_correctly = false;
    for (auto infinite_loop_guard = 1 + list.size(); infinite_loop_guard > 0;
         --infinite_loop_guard) {
      const auto& vertex_swap = list.at(current_id);
      const auto token_swap = m_token_tracker.do_vertex_swap(vertex_swap);
      const auto citer = m_data.find(token_swap);
      if (citer != m_data.cend() && citer->second != invalid_index) {
        // The swap occurred before, the entry tells us the ID.
        // Erase both swaps.
        list.erase(citer->second);
        list.erase(current_id);
        // We have to start at the beginning.
        // Changing the labels for these tokens
        // messes up other entries between the two swaps.
        // A warm restart from the middle of the swap list
        // would be a lot of extra complication, not worth it for now.
        terminated_correctly = true;
        break;
      }
      // Swap hasn't occurred before, now advance.
      m_data[token_swap] = current_id;
      const auto next_id_opt = list.next(current_id);
      if (!next_id_opt) {
        terminated_correctly = true;
        break;
      }
      current_id = next_id_opt.value();
    }
    TKET_ASSERT(terminated_correctly);
    const auto new_size = list.size();
    if (old_size == new_size) {
      return;
    }
    TKET_ASSERT(new_size < old_size);
  }
  TKET_ASSERT(!"optimise_pass_with_token_tracking termination");
}

void SwapListOptimiser::full_optimise(SwapList& list) {
  // More experimentation needed to find the best combination.
  optimise_pass_with_zero_travel(list);
  m_token_tracker.reset();
  optimise_pass_with_token_tracking_without_clearing_tracker(list);
}

void SwapListOptimiser::full_optimise(
    SwapList& list, const VertexMapping& vertex_mapping) {
  for (auto infinite_loop_guard = 1 + list.size(); infinite_loop_guard > 0;
       --infinite_loop_guard) {
    const auto old_size = list.size();
    full_optimise(list);
    optimise_pass_remove_empty_swaps(list, vertex_mapping);
    if (old_size == list.size() || list.size() == 0) {
      return;
    }
    TKET_ASSERT(list.size() < old_size);
  }
  TKET_ASSERT(!"full_optimise termination");
}

}  // namespace tsa_internal
}  // namespace tket
