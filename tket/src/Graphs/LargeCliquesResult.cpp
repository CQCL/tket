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

#include "LargeCliquesResult.hpp"

#include "AdjacencyData.hpp"

using std::set;
using std::size_t;
using std::vector;

namespace tket {
namespace graphs {

LargeCliquesResult::LargeCliquesResult(
    const AdjacencyData& adjacency_data,
    const set<size_t>& vertices_in_component, size_t internal_size_limit) {
  // ...and, before swapping, everything in here will have size one greater than
  // "result".
  vector<set<size_t>> extended_result;

  // At each stage, every set will have the same size, containing a clique of
  // that size.
  cliques.reserve(vertices_in_component.size());
  for (size_t i : vertices_in_component) {
    cliques.emplace_back(set{i});
  }
  bool hit_internal_limit = false;

  // A little trick: the vertex indices can always be stored in order: v1 < v2 <
  // ... Therefore, within each vertex set, if we only allow adding vertices
  // with LARGER index than the largest already stored, we will automatically
  // discard duplicates.
  for (size_t counter = 0; counter <= vertices_in_component.size(); ++counter) {
    // Extend the existing results.
    for (const auto& clique : cliques) {
      // Guaranteed to be nonempty.
      const size_t largest_index = *clique.crbegin();

      if (extended_result.size() >= internal_size_limit) {
        hit_internal_limit = true;
        break;
      }
      // We only have to check the neighbours of ONE vertex to form a larger
      // clique.
      const auto& neighbours = adjacency_data.get_neighbours(largest_index);
      for (size_t new_v : neighbours) {
        if (new_v <= largest_index) {
          continue;
        }
        // We have a new vertex; does it adjoin EVERY existing vertex?
        bool joins_every_vertex = true;
        for (size_t existing_v : clique) {
          if (adjacency_data.get_neighbours(existing_v).count(new_v) == 0) {
            joins_every_vertex = false;
            break;
          }
        }
        if (joins_every_vertex) {
          extended_result.emplace_back(clique);
          extended_result.back().insert(new_v);
          if (extended_result.size() >= internal_size_limit) {
            hit_internal_limit = true;
            break;
          }
        }
      }
    }
    if (extended_result.empty()) {
      // No extensions
      cliques_are_definitely_max_size = !hit_internal_limit;
      return;
    }
    cliques = extended_result;
    // TODO: improve memory reallocation etc. by reusing rather than clearing
    // (or, could be sneaky and do this within a single vector, jiggle with
    // indices).
    extended_result.clear();
  }
  cliques_are_definitely_max_size = false;
}

}  // namespace graphs
}  // namespace tket
