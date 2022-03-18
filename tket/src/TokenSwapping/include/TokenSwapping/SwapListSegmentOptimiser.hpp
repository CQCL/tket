
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

#include <cstdint>
#include <set>
#include <vector>

#include "PartialMappingLookup.hpp"
#include "TokenSwapping/SwapFunctions.hpp"
#include "VertexMapResizing.hpp"

namespace tket {
namespace tsa_internal {

/** Given a swap list and a start point in the list, uses the lookup table
 * to reduce an interval of swaps, replacing them in-place by a shorter sequence
 * with the same end-to-end vertex mapping (although source->target mappings may
 * change for empty source vertices, i.e. those without a token at the
 * beginning).
 */
class SwapListSegmentOptimiser {
 public:
  struct Output {
    /** The length of the segment that was replaced.
     * Of course, this will be zero if no optimisation takes place.
     */
    size_t initial_segment_size;

    /** The length of the segment after replacement. Always <=
     * initial_segment_size. */
    size_t final_segment_size;

    /** If we did replace a segment with a shorter one, give the ID of the last
     * swap of the segment. It might be null because the new segment might be
     * empty.
     */
    std::optional<SwapID> new_segment_last_id;
  };

  /** Starting at the given ID, which must be valid, move forward to examine an
   * interval of swaps, and try to replace it with a shorter sequence looked up
   * in the table. It MAY replace a segment with a different one of equal
   * length; optimisation has probably already taken place, and couldn't break
   * it up any further. If the table suggests a different but still valid
   * interval, it MAY afford further opportunities for optimisation even if it's
   * of the same length, so we might as well splice in the new segment.
   * @param initial_id The ID within the swap list of the first swap which may
   * be replaced, where we begin optimisation.
   * @param vertices_with_tokens_at_start Just before the swap at initial_id is
   * performed, which vertices have tokens on them? Extra unused vertices are
   * allowed (but are helpful, since they may be added into the new sequence to
   * reduce length).
   * @param map_resizing An object to add/remove vertices from the mapping, with
   * knowledge of edges in the graph (not just those involved in the swap list).
   * @param swap_list The sequence of swaps to be reduced, in-place.
   * @return An object stored internally, with information about the segment
   * replacement/reduction (if any).
   */
  const Output& optimise_segment(
      SwapID initial_id, const std::set<size_t>& vertices_with_tokens_at_start,
      VertexMapResizing& map_resizing, SwapList& swap_list);

 private:
  Output m_output;
  PartialMappingLookup m_mapping_lookup;

  // Naively, a greedy-type way to optimise is to
  // reduce the SHORTEST sequence possible, by the LARGEST amount.
  // This may not always be optimal, but should be OK.
  std::vector<Swap> m_best_optimised_swaps;

  /** Once m_output.initial_segment_size and m_best_optimised_swaps have been
   * filled, fill in the rest of the data in m_output and make the swap
   * replacements in swap_list.
   * @param initial_id The ID within the swap list of the first swap which may
   * be replaced, where we begin optimisation.
   * @param swap_list The sequence of swaps to be reduced, in-place.
   */
  void fill_final_output_and_swaplist(SwapID initial_id, SwapList& swap_list);
};

}  // namespace tsa_internal
}  // namespace tket
