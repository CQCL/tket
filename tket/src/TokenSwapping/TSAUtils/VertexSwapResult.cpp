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

#include "TokenSwapping/VertexSwapResult.hpp"

namespace tket {
namespace tsa_internal {

VertexSwapResult::VertexSwapResult(
    size_t v1, size_t v2, VertexMapping& vertex_mapping, SwapList& swap_list)
    : VertexSwapResult(v1, v2, vertex_mapping) {
  if (tokens_moved != 0) {
    swap_list.push_back(get_swap(v1, v2));
  }
}

VertexSwapResult::VertexSwapResult(
    const Swap& swap, VertexMapping& vertex_mapping)
    : VertexSwapResult(swap.first, swap.second, vertex_mapping) {}

VertexSwapResult::VertexSwapResult(
    size_t v1, size_t v2, VertexMapping& vertex_mapping) {
  if (vertex_mapping.count(v1) == 0) {
    if (vertex_mapping.count(v2) == 0) {
      tokens_moved = 0;
      return;
    }
    // No token on the first.
    vertex_mapping[v1] = vertex_mapping[v2];
    vertex_mapping.erase(v2);
    tokens_moved = 1;
    return;
  }
  // A token on the first.
  if (vertex_mapping.count(v2) == 0) {
    vertex_mapping[v2] = vertex_mapping[v1];
    vertex_mapping.erase(v1);
    tokens_moved = 1;
    return;
  }
  // Tokens on both.
  std::swap(vertex_mapping[v1], vertex_mapping[v2]);
  tokens_moved = 2;
}

}  // namespace tsa_internal
}  // namespace tket
