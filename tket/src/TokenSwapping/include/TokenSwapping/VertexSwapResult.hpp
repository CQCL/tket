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
#include <optional>
#include <utility>

#include "VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** For performing a vertex swap, and checking how many tokens moved. */
struct VertexSwapResult {
  /** How many tokens moved? Must be one of 0,1,2. */
  unsigned tokens_moved;

  /** Carry out the swap on the tokens and get the result.
   *  @param swap The swap to perform.
   *  @param vertex_mapping The source to target mapping,
   *    will be updated with the swap.
   */
  VertexSwapResult(const Swap& swap, VertexMapping& vertex_mapping);

  /** Pass in the two vertex size_t numbers directly.
   *  @param v1 First vertex of the swap to perform.
   *  @param v2 Second vertex of the swap to perform.
   *  @param vertex_mapping The source to target mapping,
   *    will be updated with the swap.
   */
  VertexSwapResult(size_t v1, size_t v2, VertexMapping& vertex_mapping);

  /** If the swap moves at least one nonempty token, carry out the swap.
   *  Otherwise, does nothing.
   *  @param v1 First vertex of the swap to perform.
   *  @param v2 Second vertex of the swap to perform.
   *  @param vertex_mapping The source to target mapping,
   *    will be updated with the swap.
   *  @param swap_list The list of swaps, which will be updated with the swap.
   */
  VertexSwapResult(
      size_t v1, size_t v2, VertexMapping& vertex_mapping, SwapList& swap_list);
};

}  // namespace tsa_internal
}  // namespace tket
