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
#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

// PERMUTATION HASH EXPLANATION:
// Certain permutations on the vertices [0,1,2,3,4,5] are represented by an
// unsigned value, the "permutation hash". In fact, ANY permutation on a list of
// ANY 6 distinct objects can be reduced to one of these by a suitable vertex
// relabelling, which is what this class is for.
//
// Note that not every permutation on [0,1,2,3,4,5] corresponds to a permutation
// hash (in fact, very few do); most of them still need relabelling, just as for
// arbitrary labels.
//
// The permutation hashes are done by taking a partition of 6, with parts in
// decreasing order, e.g. 6 = 3+2+1 = 4+2 = 3+3 = 2+2+1+1, etc. We remove all
// the 1 entries, and stick the digits together into a single decimal:
//
//   3+2+1 -> 32,   4+2 -> 42,   2+2+1+1 -> 22
//
// Each digit represents the length of a slice of the 6 elements 012345:
//
// 32 -> (012)(34)(5),    42 -> (0123)(45),    22 -> (01)(23)(4)(5).
//
// The notation (abcd) represents a cyclic shift on the elements a,b,c,d.
// Thus a -> b -> c -> d -> a.
//
// EVERY permutation on 6 ARBITRARY distinct objects is equivalent to one of
// these, after suitable vertex relabelling. This follows because permutations
// can be decomposed into disjoint cycles.
//

/** Given a permutation with arbitrary vertex labels, currently size <= 6, we
 * want to relabel the vertices so that we can look up an isomorphic mapping in
 * a table. This class gives one possible way. Still some scope for research and
 * improvement here; we want to cut down the number of "isomorphic" copies as
 * much as possible (whatever "isomorphic" means in this context) to make the
 * lookup table fairly small.
 */
class CanonicalRelabelling {
 public:
  /** For looking up mappings in the table. */
  struct Result {
    /** Will be empty if there are too many vertices. (Current limit is 6,
     * although this may be updated in future). */
    VertexMapping old_to_new_vertices;

    /** Element[i], for new vertex i, is the old vertex number which corresponds
     * to i. Empty if too many vertices.
     */
    std::vector<size_t> new_to_old_vertices;

    /** Set equal to zero if too many vertices. Any permutation on <= 6 vertices
     * is assigned a number, to be looked up in the table. 0 is the identity
     * permutation. */
    unsigned permutation_hash;

    /** Were there too many vertices in the mapping to look up in the table? */
    bool too_many_vertices;

    /** Was it the identity mapping? If so, no need to relabel OR look up in a
     * table. */
    bool identity;
  };

  /** The returned Result object is stored internally.
   * @param desired_mapping A (source vertex) -> (target vertex) permutation on
   * arbitrary vertex labels.
   * @return An object withe information such as (1) how to relabel vertices;
   * (2) The permutation on NEW vertices, for looking up in a table.
   */
  const Result& operator()(const VertexMapping& desired_mapping);

  CanonicalRelabelling();

 private:
  Result m_result;

  VertexMapping m_desired_mapping;
  VertexMapping m_work_mapping;

  /** The relabelling/permutation hashing is all based upon decomposing an
   * arbitrarily labelled permutation into disjoint cycles, then relabelling the
   * vertices within the cycles in a reasonable way.
   */
  std::vector<std::vector<size_t>> m_cycles;

  /** The indices in "m_cycles" after sorting appropriately. */
  std::vector<unsigned> m_sorted_cycles_indices;
};

}  // namespace tsa_internal
}  // namespace tket
