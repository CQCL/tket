#include "VertexSwapResult.hpp"

;

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
