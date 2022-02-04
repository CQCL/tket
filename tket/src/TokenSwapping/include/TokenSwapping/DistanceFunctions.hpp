// Copyright 2019-2021 Cambridge Quantum Computing
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

#include <optional>
#include <set>
#include <stdexcept>

#include "TokenSwapping/DistancesInterface.hpp"
#include "VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** The sum of the distances of each nonempty token to its home.
 *  (This is also referred to as "L" in various places, coming from the 2016
 *  paper "Approximation and Hardness of Token Swapping").
 *  @param vertex_mapping (current vertex where a token lies)->(target vertex)
 * mapping.
 *  @param distances An object to calculate distances between vertices.
 *  @return the sum, over all tokens, of (current vertex)->(target vertex)
 * distances.
 */
size_t get_total_home_distances(
    const VertexMapping& vertex_mapping, DistancesInterface& distances);

/** For just the abstract move v1->v2, ignoring the token on v2,
 *  by how much does L (the total distances to home) decrease?
 *  @param vertex_mapping current source->target mapping.
 *  @param v1 First vertex.
 *  @param v2 Second vertex. Not required to be adjacent to v1.
 *  @param distances An object to calculate distances between vertices.
 *  @return The amount by which L = get_total_home_distances would decrease,
 *    IF we moved the token on v1 to v2, IGNORING the token currently on v2
 *    (which of course is impossible to do in reality if there is a token on
 * v2), and leaving all other tokens unchanged. Doesn't have to be positive, of
 * course, although positive numbers are good.
 */
int get_move_decrease(
    const VertexMapping& vertex_mapping, size_t v1, size_t v2,
    DistancesInterface& distances);

/** The same as get_move_decrease, but for an abstract swap(v1,v2).
 *  @param vertex_mapping current source->target mapping.
 *  @param v1 First vertex.
 *  @param v2 Second vertex. Not required to be adjacent to v1.
 *  @param distances An object to calculate distances between vertices.
 *  @return The amount by which L = get_total_home_distances would decrease,
 *    (which does not have to be a positive number),
 *    IF the tokens currently on v1,v2 were swapped, and all other tokens
 *    left unchanged.
 */
int get_swap_decrease(
    const VertexMapping& vertex_mapping, size_t v1, size_t v2,
    DistancesInterface& distances);

/** A simple theoretical lower bound on the number of swaps necessary
 *  to achieve a given vertex mapping. (Of course it is not always possible
 *  to achieve this bound. But the algorithm in the 2016 paper
 *  "Approximation and Hardness of Token Swapping", for example, guarantees
 *  to find a solution within a factor of 4, or a factor of 2 for trees,
 *  in the case where every vertex has a token).
 *  TODO: What happens if some vertices are empty? Not considered in the 2016
 *  paper! Need to think about it. This is still a lower bound, but how close?
 *  @param vertex_mapping current source->target mapping.
 *  @param distances An object to calculate distances between vertices.
 *  @return A number S such that every possible solution has >= S swaps.
 *    However, note that the true minimum value might be larger, but finding
 *    the value seems about as hard as finding an actual solution, and thus
 *    is possibly exponentially hard (seems to be unknown, even for trees).
 */
size_t get_swaps_lower_bound(
    const VertexMapping& vertex_mapping, DistancesInterface& distances);

}  // namespace tsa_internal
}  // namespace tket
