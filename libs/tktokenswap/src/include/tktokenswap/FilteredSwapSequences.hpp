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

#include <map>
#include <vector>

#include "SwapConversion.hpp"

namespace tket {
namespace tsa_internal {

/** Takes a raw list of integers, where each integer represents a swap sequence
 * on the vertices {0,1,2,...,5} giving the same vertex permutation.
 *
 * NOTE: the magic number 5 (or 6) arises because we originally constructed
 * the table by exhaustively constructing swap sequences on graphs with up to
 * 6 vertices, up to a certain length. [Results were also merged together,
 * e.g. the cycle C_6, or with a few extra edges added, can be searched
 * in reasonable time to a longer length than K_6].
 * This was chosen because the complete graph K_6 has 15 edges,
 * so conveniently each edge (or swap) can be represented by a number 1-15,
 * and thus by a single hexadecimal digit.
 * Thus, 4 bits are needed for each swap, so a 64-bit integer can represent
 * swap sequences of length <= 16 (with 0 denoting the end of sequence).
 * [Although, the table currently has entries only of length <= 12].
 * [Actually, it is not hard to prove - by considering "token tracking" -
 * that optimal swap sequences on <= N vertices have
 * length <= N(N-1)/2, the same as the number of edges of K_N. Thus length
 * <= 15 already suffices to represent all possible optimal sequences
 * on <= 6 vertices].
 * If we used 5 bits, we'd be able to represent sequences of length <= 12
 * (because 5*12 = 60 < 64) on graphs with <= 8 vertices (since
 * 8*7/2 = 28 < 31).
 * If we expand the table in future, we will probably design a whole new
 * format, so we don't attempt to make it more generic at this stage.
 *
 * Given such data, FilteredSwapSequences knows how to index and store it
 * somehow (exactly how is an implementation detail - it can be thought of
 * as a "database of swap sequences"),
 * so that results can be looked up again, when given the edges bitset
 * (i.e., edges existing in the graph, i.e. vertex swaps we are allowed to
 * perform). This is for data close to the raw table data; it knows nothing
 * about vertex relabelling, which of course is a crucial component.
 *
 * The main precomputed table of data is also accessed here, via the
 * SingleSequenceData constructor.
 *
 * Note that the raw table contains several lists of integers,
 * each one denoting different swap sequences enacting a single permutation, but
 * with different edges; whereas this class only stores a single list in
 * searchable form.
 */
class FilteredSwapSequences {
 public:
  /** A result which comes from the "raw" table data in SwapSequenceTable, with
   * minimal processing. */
  struct SingleSequenceData {
    /** The edges (i.e., swaps) actually used (or 0 if none are used). [This
     * could be computed from swaps_code but there is no need to recompute each
     * time]. */
    SwapConversion::EdgesBitset edges_bitset;

    /** An integer encoding a sequence of swaps. 0 means no swaps. */
    SwapConversion::SwapHash swaps_code;

    /** The number of swaps used. Set to max() if no valid sequence was found
     * (e.g., if not present in the table). */
    unsigned number_of_swaps;

    /** Initialised with "null" values automatically, i.e. number_of_swaps
     * taking value max(). */
    SingleSequenceData();

    /** This is how we access the fixed data in the large const static global
     * table. This constructor looks up the shortest sequence of swaps enacting
     * the given permutation, and fills the entries.
     *  @param permutation_hash The hash of the desired permutation of
     * {0,1,2,...,5}, as used to look up results in the table (after
     * relabelling). See CanonicalRelabelling for explanation.
     *  @param edges_bitset The collection of edges on {0,1,2,...,5} which
     * actually exist in the graph (i.e., the swaps which are allowed).
     *  @param max_number_of_swaps Do not return any solutions with more swaps
     * than this: useful speedup to allow early termination.
     */
    SingleSequenceData(
        unsigned permutation_hash, SwapConversion::EdgesBitset edges_bitset,
        unsigned max_number_of_swaps);
  };

  /** Index and process the raw data to allow later retrieval. Can only be done
   * once (a security measure to avoid accidentally reconstructing large tables
   * multiple times). The codes don't need to be sorted OR deduplicated.
   * Duplicate, redundant and suboptimal data IS tolerated, as long as it is
   * correct. Such data could lead to slowdowns from a larger table, BUT will
   * not affect the actual results (i.e., if the data contains some entries
   * inferior to others, then the inferior results will automatically never be
   * returned, because the superior ones will always be found).
   * @param codes The raw list of integers stored in the original table
   */
  void initialise(std::vector<SwapConversion::SwapHash> codes);

  /** Search for the entry with fewest swaps whose edges_bitset is a
   * subset of the given edges_bitset (so that it only uses allowed swaps).
   * If there is no suitable sequence in the table, returns a null object.
   * Stop searching early if it finds that all entries have too many swaps.
   *  @param allowed_swaps The swaps which can occur (in other words, the
   * existing edges in the graph).
   *  @param max_num_swaps Don't return any entries with more than this many
   * swaps.
   *  @return An entry with the fewest swaps, or a null entry if none exists.
   */
  SingleSequenceData get_lookup_result(
      SwapConversion::EdgesBitset allowed_swaps, unsigned max_num_swaps) const;

  /** For testing, just count how many entries we've stored.
   *  @return The total number of encoded swap sequences stored internally.
   */
  size_t get_total_number_of_entries() const;

 private:
  /** We recalculate the number of swaps each time, rather than storing.
   * We just sort by swaps_code, since this respects numbers of swaps.
   * I.e., if S1, S2 are swap sequences, and encoding(S(j)) is an integer, then
   * length(S1) < length(S2)   =>   encoding(S1) < encoding(S2).
   * Thus, minimising encoding(S) will also force minimising length(S).
   */
  struct TrimmedSingleSequenceData {
    SwapConversion::EdgesBitset edges_bitset;
    SwapConversion::SwapHash swaps_code;
  };

  /** Key: a subset of bits in edges_bitset.
   * Value: codes containing those bits in their edges bitset, sorted in
   * increasing order. No entry occurs multiple times, but the values are spread
   * out amongst the keys to balance the data better and give faster lookup.
   */
  std::map<SwapConversion::EdgesBitset, std::vector<TrimmedSingleSequenceData>>
      m_internal_data;

  /** Must be pushed back in increasing order of swaps_code. Processes and
   * stores the result for later searchability.
   * @param datum Information about a single raw entry from the table.
   */
  void push_back(TrimmedSingleSequenceData datum);
};

}  // namespace tsa_internal
}  // namespace tket
