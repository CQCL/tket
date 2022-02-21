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

#include <string>

#include "TokenSwapping/DistanceFunctions.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** Get a string representation.
 *  @param vertex_mapping A mapping, usually representing a desired
 * source->target mapping for a Token Swapping problem.
 *  @return A string representation.
 */
std::string str(const VertexMapping& vertex_mapping);

/** Get a string representation.
 *  @param swaps An ordered list of swaps, usually the solution to a Token
 * Swapping problem.
 *  @return A string representation.
 */
std::string str(const SwapList& swaps);

/** Get a string representation.
 *  @param swaps An ordered list of swaps, usually the solution to a Token
 * Swapping problem.
 *  @return A string representation.
 */
std::string str(const std::vector<Swap>& swaps);

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
