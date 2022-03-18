
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

#include <set>

#include "PartialMappingLookup.hpp"
#include "SwapListSegmentOptimiser.hpp"
#include "TokenSwapping/SwapListOptimiser.hpp"
#include "VertexMapResizing.hpp"

/// TODO: The swap table optimiser currently tries to optimise many segments;
/// solving ~2300 problems with Best TSA takes ~20 seconds, most of which
/// is the table optimisation part.
/// Certainly we can cut down the number of segments optimised;
/// needs experimentation.

namespace tket {
namespace tsa_internal {

/** Uses the lookup table to reduce many intervals of a swap sequence. */
class SwapListTableOptimiser {
 public:
  /** Reduce the given list of swap in-place, by using the big lookup table.
   * Swaps may be significantly reordered, and the final end-to-end
   * permutation of vertices may change; only the partial mapping of those
   * vertices with tokens is preserved. It's not actually clear what the best
   * method is; experimentation is still needed. We can optimise any segment,
   * i.e. between any two points. But then, which other segments should we
   * choose? Should we overlap trial segments? Should we then combine with the
   * simple SwapListOptimiser again? We are lacking a lot of theory to guide us.
   * This pass will erase some empty swaps, but doesn't guarantee to find all
   * (although in practice, it never does produce empty swaps, if they were
   * previously well optimised with a swap list optimiser. Is this "luck", or is
   * there a theoretical reason?)
   * @param vertices_with_tokens_at_start Before we perform any swaps, which
   * vertices have tokens on them? Other vertices are allowed to be moved around
   * arbitrarily.
   * @param map_resizing An object to take a VertexMapping and enlarge/contract
   * it to give the desired number of vertices. So, this object knows about the
   * edges in the graph.
   * @param swap_list The sequence of swaps to be shortened.
   * @param swap_list_optimiser An object to handle non-table optimisations.
   * This is used only to do the basic passes needed to make the table effective
   * (i.e., clustering interacting swaps together).
   */
  void optimise(
      const std::set<size_t>& vertices_with_tokens_at_start,
      VertexMapResizing& map_resizing, SwapList& swap_list,
      SwapListOptimiser& swap_list_optimiser);

  /** For testing, give internal access to the segment optimiser.
   * @return a reference to the internal segment optimiser object.
   */
  SwapListSegmentOptimiser& get_segment_optimiser();

 private:
  SwapListSegmentOptimiser m_segment_optimiser;

  /** The same interface as "optimise", which goes in both directions,
   * and calls this function in a loop, repeatedly reversing and re-reversing
   * the swap list to do both directions. A bit crude, but simple and not
   * actually too inefficient.
   * @param @param vertices_with_tokens_at_start Before we perform any swaps,
   * which vertices have tokens on them?
   * @param map_resizing An object to take a VertexMapping and enlarge/contract
   * it to give the desired number of vertices.
   * @param swap_list The sequence of swaps to be shortened.
   * @param swap_list_optimiser An object to handle non-table optimisations.
   */
  void optimise_in_forward_direction(
      const std::set<size_t>& vertices_with_tokens_at_start,
      VertexMapResizing& map_resizing, SwapList& swap_list,
      SwapListOptimiser& swap_list_optimiser);
};

}  // namespace tsa_internal
}  // namespace tket
