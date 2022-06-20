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

#include "CyclesGrowthManager.hpp"
#include "PartialTsaInterface.hpp"

namespace tket {
namespace tsa_internal {

/** Concerned with filtering and selecting candidate cycles
 *  to convert into a swap sequence. Used by CyclesPartialTsa.
 *  For further explanation, please see the comments for the class
 *  CyclesPartialTsa.
 *
 *  This is used when all cycles are valid candidates to be converted
 *  into swap sequences. This class selects the ones to use.
 *  All cycle candidates are assumed to have the same length
 *  (swaps are just cycles on 2 vertices), but have different "power",
 *  i.e. different overall contribution to the decrease of L, the sum of
 *  the distances between the current vertex of each token and its target.
 *
 *  We only want to return solutions which strictly decrease L, so that
 *  we're guaranteed to make progress (or make no change).
 *  We must select a subset of disjoint cycles, since if they
 *  were not disjoint, the returned solution might not decrease L.
 *  (We based all our calculations on treating the cycles individually,
 *  so obviously non-disjoint cycles could behave very differently).
 */
class CyclesCandidateManager {
 public:
  /** These control the behaviour of filtering for candidate selection.
   *  Experimentation needed to find the best options.
   */
  struct Options {
    // In both these options, we have a whole collection of candidate
    // swap sequences.
    // We can EITHER perform just the best single candidate,
    // OR carry out multiple swap sequences simultaneously,
    // by selecting a large disjoint subset.
    // However, returning multiple sequences, although probably faster
    // to compute overall, might give a worse end-to-end solution
    // (but this needs testing). (But of course it may actually be slower.
    // All these are just guesses, need testing!)
    // The reason is that, once the tokens
    // have shifted a little bit, it may enable better solutions
    // (sequences of higher power) which the algorithm previously
    // did not detect.

    /** Setting this to "false" means that only the best single swaps
     *  will be returned, the others being discarded. (E.g., if some swaps
     *  move two tokens closer to home, i.e. have "power" two, then
     *  "power one" swaps - those which only move one token closer to home,
     *  the other token being empty, or remaining at the same distance from
     *  its target - will be discarded).
     */
    bool return_all_good_single_swaps = false;

    /** The same as "return_all_good_single_swaps", but for cycles
     *  on >= 3 vertices. Do we return ALL cycle solutions, or only those
     *  which decrease L by the largest amount?
     */
    bool return_lower_power_solutions_for_multiswap_candidates = false;

    /** The "power" of a swap sequence is (total L decrease) / (number of
     * swaps). Since a swap can change L by -2,-1,0,1,2 (since up to 2 tokens
     * are moved one step), always |power| <= 2. But let's assume that negative
     * power candidates are discarded, and rescale to be a percentage. Discard
     * all candidates with power percentage smaller than this. Note that only
     * fairly dense problems (lots of tokens, or all clustered close together)
     * are likely to give higher powers; if all tokens are far apart, or there
     * are very few of them, then swapping two nonempty tokens is rare, so
     * immediately most candidates would not expect to reach even 50% power.
     */
    unsigned min_candidate_power_percentage = 0;
  };

  /** The "CyclesGrowthManager" object stores the candidate cycles internally,
   *  then we select the set of candidates to use, convert them into swaps,
   *  and append them to the list of swaps. (All distance data has already
   *  been calculated and ctored within the cycles).
   *
   *  @param growth_manager The object containing the candidate cycles
   *  @param swaps The list of swaps we will add to, once we convert
   *    the candidates into swaps.
   *  @param vertex_mapping The current vertex->target vertex mapping,
   *    which will be updated with the added swaps.
   */
  void append_partial_solution(
      const CyclesGrowthManager& growth_manager, SwapList& swaps,
      VertexMapping& vertex_mapping);

 private:
  Options m_options;

  /** Information about the stored candidates, for filtering. */
  struct CycleData {
    Cycles::ID id;

    /** The vertices are listed in a vector.
     *  Store the index, in the vector, of the lowest valued vertex.
     *  The purpose is to detect duplicate stored cycles (starting from
     *  a different vertex) and discard all but one of them.
     *  (Unfortunately necessary because, as cycles are being built up,
     *  we don't know which final vertices will occur, so we can get many
     *  duplicate subpaths. Is there a clever data structure to improve this?)
     */
    size_t first_vertex_index;
  };

  /** Key: a hash of the vertices in the cycle
   *  Value: information about the candidate cycles of the last cycle
   *  with that hash. (Hash collisions are expected to be very rare, and they
   *  cause no actual problem, so it's probably faster not to use complete
   *  buckets to resolve hash collisions).
   *  Used to find duplicate cycles (the same vertices in the same cyclic
   *  order, but with different start vertex in the vector).
   */
  std::map<size_t, CycleData> m_cycle_with_vertex_hash;

  /** We will discard duplicate cycles. For better constness, we don't delete
   *  cycles, we just store the IDs of those ones we want to use.
   */
  std::vector<Cycles::ID> m_cycles_to_keep;

  /** Key: a cycle ID
   *  Value: how many other cycles it touches (i.e., cycles sharing a vertex
   *    with it, so not disjoint).
   *  This will be used to select a large subset of pairwise disjoint
   *  cycles, with a simple greedy algorithm.
   */
  std::map<Cycles::ID, size_t> m_touching_data;

  /** Used by should_add_swaps_for_candidate, to see whether a cycle
   *  is disjoint from those already selected.
   */
  std::set<size_t> m_vertices_used;

  /** Fills m_cycles_to_keep (so, effectively discarding unsuitable cycles),
   *  returns the common cycle length.
   *  @param cycles The complete collection of candidate cycles.
   *  @return The number of vertices in each cycle
   *    (all cycles should be the same length).
   */
  size_t fill_initial_cycle_ids(const Cycles& cycles);

  /** Updates m_cycles_to_keep. Keep only those solutions with the
   *  highest L-decrease.
   *  @param cycles The complete collection of candidate cycles,
   *    but we already have filled m_cycles_to_keep so will
   *    only consider those cycles.
   */
  void discard_lower_power_solutions(const Cycles& cycles);

  /** Sorts m_cycles_to_keep so that those which touch
   *  the fewest other cycles are listed first.
   *  @param cycles The complete collection of candidate cycles,
   *     but we only consider those cycles with IDs in m_cycles_to_keep.
   */
  void sort_candidates(const Cycles& cycles);

  /** Checks if the candidate is disjoint from all other candidates
   *  currently used (stored in m_vertices_used). If so updates
   *  m_vertices_used and returns true (but takes no other action).
   *  Otherwise, do nothing and return false.
   *  @param cycles The complete collection of candidate cycles.
   *  @param id The single cycle under consideration.
   *  @return whether this single cycle should be added to the collection
   *    of candidates.
   */
  bool should_add_swaps_for_candidate(const Cycles& cycles, Cycles::ID id);
};

}  // namespace tsa_internal
}  // namespace tket
