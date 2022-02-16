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

#include "DynamicTokenTracker.hpp"

namespace tket {
namespace tsa_internal {

void DynamicTokenTracker::clear() { m_vertex_to_token.clear(); }

void DynamicTokenTracker::reset() {
  for (auto& entry : m_vertex_to_token) {
    entry.second = entry.first;
  }
}

Swap DynamicTokenTracker::do_vertex_swap(const Swap& swap) {
  const auto v1 = swap.first;
  const auto v2 = swap.second;
  const auto t1 = get_token_at_vertex(v1);
  const auto t2 = get_token_at_vertex(v2);
  m_vertex_to_token[v1] = t2;
  m_vertex_to_token[v2] = t1;
  return get_swap(t1, t2);
}

bool DynamicTokenTracker::equal_vertex_permutation_from_swaps(
    const DynamicTokenTracker& other) const {
  return tokens_here_have_equal_locations_in_the_other_object(other) &&
         other.tokens_here_have_equal_locations_in_the_other_object(*this);
}

bool DynamicTokenTracker::tokens_here_have_equal_locations_in_the_other_object(
    const DynamicTokenTracker& other) const {
  for (const auto& vertex_token_pair : m_vertex_to_token) {
    const auto vertex = vertex_token_pair.first;
    const auto token = vertex_token_pair.second;
    const auto citer = other.m_vertex_to_token.find(vertex);

    if (citer == other.m_vertex_to_token.cend()) {
      // If it's unmentioned by the other, then the vertex MUST be fixed
      // to give the same permutation.
      // Otherwise, the other object doesn't know where the token moved to.
      if (vertex != token) {
        return false;
      }
    } else {
      if (token != citer->second) {
        return false;
      }
    }
  }
  return true;
}

size_t DynamicTokenTracker::get_token_at_vertex(size_t vertex) {
  const auto iter = m_vertex_to_token.find(vertex);
  if (iter == m_vertex_to_token.end()) {
    m_vertex_to_token[vertex] = vertex;
    return vertex;
  }
  return iter->second;
}

}  // namespace tsa_internal
}  // namespace tket
