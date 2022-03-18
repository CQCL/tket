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
#include <map>
#include <vector>

namespace tket {
namespace tsa_internal {

/** For swaps on vertices {0,1,2,...,5}, return precomputed short swap
 * sequences using given sets of edges. Should be close to optimal.
 * (Every sequence should have the joint-shortest length amongst all sequences
 * using those particular swaps, but not every possible sequence is included).
 *
 * (The only possibility of non-optimality is that some solutions
 * using many edges might be missing. It was constructed using breadth-first
 * searches of all possible sequences up to certain depths on various graphs
 * with <= 6 vertices. Due to time/space limitations some non-complete graphs
 * were searched as well as complete graphs K4, K5, K6.
 *
 * Note that, by token tracking, any swap sequence of n vertices of length
 * > n(n-1)/2 can be reduced in length, so in fact any optimal swap sequence
 * on n vertices has length <= n(n-1)/2, the number of edges in
 * the complete graph K(n).
 *
 * Of course, ideally we'd search K6 up to depth 15, but searching up to depth 9
 * already consumed ~30 mins of CPU time and most of the memory capacity of an
 * ordinary laptop. More efficient exhaustive search algorithms with clever
 * pruning might cut it down a bit, but (since each added depth increases the
 * difficulty roughly by a factor of 14) it would require significant
 * computational effort to reach even depth 12 for K6, and depth 15 probably
 * requires a supercomputer, or a very large distributed computation,
 * or significantly more intelligent methods).
 *
 * The table size is far smaller than the precomputation needed to create it.
 * The creation considered millions of sequences, but the table has only a few
 * thousand entries.
 *
 * The table currently contains ALL optimal swap sequences on <= 5 vertices,
 * and also all swap sequences of length:
 *     <= 9 on 6 vertices (K6, depth 9);
 *     <= 12 on cycles with <= 6 vertices (C5, C6);
 *     <= 12 on a few other special graphs with 6 vertices.
 *
 * Superficially redundant solutions have been removed:
 *
 * (a): If sequences S1, S2 have equal length but the edges set E1 is a subset
 * of E2, keep only S1, since every graph allowing S2 would also allow S1.
 *
 * (b): If sequences S1, S2 have len(S1) < len(S2), keep S2 exactly when E2 is
 * NOT a subset of E1 (since, then there are graphs containing E2 which do NOT
 * contain E1, so that S2 may be possible when S1 is impossible).
 *
 * Finally, to save space, every sequence was checked before insertion, and
 * inserted ONLY if its inverse was not already present in the table (since
 * inverting permutations is trivial for swaps: just reverse the order). Hence,
 * the table is only about half the size that it would otherwise be.
 *
 * But, whilst these sequences are universally valid,
 * this class knows nothing about HOW to look up results in the table
 * efficiently. The current lookup algorithms are quite crude (but actually
 * faster than fancier algorithms for this table size), but there is some
 * possibility of speedup (although not result improvements) if a really fancy
 * search/filtering algorithm can be found.
 *
 * NOTE: the format is reasonable, but still not as compressed as possible;
 * it still contains multiple isomorphic entries. A more complicated hashing
 * scheme is required to cut down on these isomorphic copies. (E.g., perm hash
 * 2, meaning the mapping 0->1, 1->0, i.e. (01), contains 0x262, 0x484, 0x737,
 * meaning swap sequences [02 12 02], [04 14 04], [13 03 13]. It is easily seen
 * that all 3 are isomorphic. The first two are of the form [ab cb ab] == [ac],
 * and the third has the form [ab cb ab] == [ca].) It seems like we'd need a
 * scheme involving integer hashing of graphs, with few isomorphic collisions,
 * but such algoritms need to be pretty simple and fast or they're not worth
 * doing except for much larger table sizes.
 */
struct SwapSequenceTable {
  /** The integer type used to encode a swap sequence on vertices {0,1,2,3,4,5}.
   */
  typedef std::uint_fast64_t Code;

  /** The KEY is a "permutation hash", i.e. a number representing a permutation
   * on {0,1,2,3,4,5}. (Not all possible permutations are represented, though;
   * suitable vertex relabelling changes many different permutations to the same
   * hash).
   *
   * See CanonicalRelabelling.hpp, SwapConversion.hpp for more explanation.
   *
   * The VALUE is a list of integers encoding a swap sequence, which all induce
   * the permutation on {0,1,2,3,4,5} with the given hash.
   * (Obviously, different sequences are allowed, because some swaps might not
   * be possible, i.e. the graph might be incomplete).
   */
  typedef std::map<unsigned, std::vector<Code>> Table;

  /** The actual large precomputed table. The entries are already sorted
   * and duplications/redundancies/suboptimality have been removed.
   * However, currently this raw data is processed by
   * FilteredSwapSequences which tolerates such imperfections.
   * Thus it is easy to add more sequences to the table without worrying
   * about them (as long as the newly added data is actually correct).
   *  @return A large precomputed raw table of data.
   */
  static Table get_table();
};

}  // namespace tsa_internal
}  // namespace tket
