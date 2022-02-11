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

/** Tracks which token is on which vertex;
 *  every vertex has a different token.
 *  Only intended for a specific optimisation pass in SwapListOptimiser.
 *  Does not require contiguous vertex numbers or token numbers.
 *  Does not require to be initialised with all vertices at the start.
 *  Thus, operations take time O(log N), with N being the current number
 *  of vertices seen, NOT the total number of vertices.
 *  The tokens are "artificial", i.e. nothing to do with an actual
 *  Token Swapping problem; they are there to track full vertex mappings
 *  induced by a sequence of swaps.
 */
class DynamicTokenTracker {
 public:
  /** Call before starting a new sequence of swaps. */
  void clear();

  /** Logically the same effect as clear, but doesn't actually clear.
   *  Instead, fills existing map entries.
   *  Should be a bit faster for many reuses than clearing every time,
   *  because it will need fewer tree rebalances inside the maps.
   */
  void reset();

  /** Swap the tokens at the given vertices,
   *  and return the TOKENS that were swapped.
   *  Note that every vertex is assumed initially to have a token
   *  with the same vertex value (i.e., the token equals the INITIAL
   *  vertex). Thus we don't need to know in advance which vertices
   *  exist, they will be lazily stored only when needed.
   *  @param swap The two vertices to be swapped.
   *  @return The two TOKENS on those vertices which were swapped.
   */
  Swap do_vertex_swap(const Swap& swap);

  /** Checks if the swap sequence performed on the other tracker object
   *  results in the same vertex permutation.
   *  This is NOT the same as just checking equality of data,
   *  because a vertex could be unmentioned in our sequence,
   *  and thus not appear anywhere internally; but in the other sequence
   *  it could appear, but end up back where it started.
   *  @param other Another DynamicTokenTracker object
   *  @return Whether the swaps performed on this object and the other object
   *    resulted in the same vertex permutation on the whole graph
   *    (remembering that some vertices may be mentioned in one object
   *    but not the other).
   */
  bool equal_vertex_permutation_from_swaps(
      const DynamicTokenTracker& other) const;

 private:
  VertexMapping m_vertex_to_token;

  /** Get the token, but if it doesn't already exist, create it.
   *  @param vertex The vertex
   *  @return The token at that vertex, or equal to the vertex number
   *    IF it doesn't already exist.
   */
  size_t get_token_at_vertex(size_t vertex);

  /** Does every token mentioned in this object lie at the same vertex in
   *  the other object?
   *  @param other Another DynamicTokenTracker object
   *  @return Whether all tokens mentioned by this object have
   *    the same location according to the other object (remembering
   *    that unmentioned vertices are implicitly assumed to have equal tokens
   *    lying on them initially).
   */
  bool tokens_here_have_equal_locations_in_the_other_object(
      const DynamicTokenTracker& other) const;
};

}  // namespace tsa_internal
}  // namespace tket
