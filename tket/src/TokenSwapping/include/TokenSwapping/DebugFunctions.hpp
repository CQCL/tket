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

#include <string>

#include "VertexMappingFunctions.hpp"

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

}  // namespace tsa_internal
}  // namespace tket
