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

#include "GraphRoutines.hpp"

#include <stack>

#include "AdjacencyData.hpp"

using std::map;
using std::set;
using std::size_t;
using std::vector;

namespace tket {
namespace graphs {

vector<set<size_t>> GraphRoutines::get_connected_components(
    const AdjacencyData& adjacency_data) {
  set<size_t> vertices_seen;
  vector<set<size_t>> result;
  const size_t number_of_vertices = adjacency_data.get_number_of_vertices();

  for (size_t v = 0; v < number_of_vertices; ++v) {
    if (vertices_seen.count(v) != 0) {
      continue;
    }
    set<size_t> vertices_seen_in_this_component;
    vertices_seen_in_this_component.insert(v);

    std::stack<size_t> vertices_to_examine_next;
    vertices_to_examine_next.push(v);

    while (!vertices_to_examine_next.empty()) {
      const size_t this_vertex = vertices_to_examine_next.top();
      vertices_to_examine_next.pop();
      const auto& neighbouring_vertices =
          adjacency_data.get_neighbours(this_vertex);

      for (size_t j : neighbouring_vertices) {
        if (vertices_seen_in_this_component.count(j) == 0) {
          vertices_to_examine_next.push(j);
          vertices_seen_in_this_component.insert(j);
        }
      }
    }
    // Now, push back the vertices seen.
    result.emplace_back(vertices_seen_in_this_component);
    for (size_t k : vertices_seen_in_this_component) {
      vertices_seen.insert(k);
    }
  }
  return result;
}

}  // namespace graphs
}  // namespace tket
