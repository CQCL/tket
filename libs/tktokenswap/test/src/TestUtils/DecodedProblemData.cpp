// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "DecodedProblemData.hpp"

#include <catch2/catch_test_macros.hpp>
#include <tktokenswap/GeneralFunctions.hpp>
#include <tktokenswap/VertexSwapResult.hpp>
#include <algorithm>

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

/// TODO: move somewhere more appropriate.
// initially, "vm" has keys equal to the vertices with tokens; the values are
// ignored. Change to the desired source->target mapping, as used in all problem
// solving, induced by the swaps. Return the number of empty swaps.
static unsigned get_problem_mapping(
    VertexMapping& vm, const vector<Swap>& swaps) {
  const auto init_num_tokens = vm.size();
  for (auto& entry : vm) {
    entry.second = entry.first;
  }
  unsigned empty_swaps = 0;
  for (auto swap : swaps) {
    const VertexSwapResult result(swap, vm);
    if (result.tokens_moved == 0) {
      ++empty_swaps;
    }
  }
  // Each time we had v1->t1, v2->t2 and we swapped v1,v2, we then got v1->t2,
  // v2->t1. Thus, the KEY is a vertex, the VALUE is the token currently on that
  // vertex. So, the VALUES are the tokens, which are the vertex it originally
  // came from, i.e., it's end vertex -> original vertex. So our desired problem
  // mapping source -> target is the REVERSE!!
  vm = get_reversed_map(vm);
  REQUIRE(init_num_tokens == vm.size());
  check_mapping(vm);
  return empty_swaps;
}

static const std::string& encoding_chars() {
  static const std::string chars{
      "0123456789abcdefghijklmnopqrstuvwxyz"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
  return chars;
}

static std::map<unsigned char, std::size_t> get_char_to_vertex_map_local() {
  std::map<unsigned char, std::size_t> char_to_vertex_map;
  const auto& chars = encoding_chars();
  for (std::size_t ii = 0; ii < chars.size(); ++ii) {
    char_to_vertex_map[chars[ii]] = ii;
  }
  return char_to_vertex_map;
}

static const std::map<unsigned char, std::size_t>& char_to_vertex_map() {
  static const std::map<unsigned char, std::size_t> map(
      get_char_to_vertex_map_local());
  return map;
}

DecodedProblemData::DecodedProblemData(
    const std::string& str,
    RequireContiguousVertices require_contiguous_vertices) {
  if (str.empty()) {
    return;
  }

  unsigned index = 0;
  bool separator_found = false;
  while (index < str.size()) {
    if (str.at(index) == '_') {
      ++index;
      separator_found = true;
      break;
    }
    const auto v1 = char_to_vertex_map().at(str.at(index));
    const auto v2 = char_to_vertex_map().at(str.at(index + 1));
    swaps.emplace_back(get_swap(v1, v2));
    index += 2;
  }

  std::set<std::size_t> vertices;
  for (auto swap : swaps) {
    vertices.insert(swap.first);
    vertices.insert(swap.second);
  }
  CHECK(vertices.size() >= 4);
  number_of_vertices = vertices.size();
  if (require_contiguous_vertices == RequireContiguousVertices::YES) {
    REQUIRE(*vertices.crbegin() + 1 == vertices.size());
  }

  // Now set up the vertex mapping. Initially, all vertices with tokens
  // have a token value equal to the vertex number.
  vertex_mapping.clear();
  if (separator_found) {
    unsigned num_tokens = 0;
    for (; index < str.size(); ++index) {
      const auto vv = char_to_vertex_map().at(str.at(index));
      if (require_contiguous_vertices == RequireContiguousVertices::YES) {
        // It came from a swap sequence. Therefore, there are no extra edges,
        // so every vertex must exist on a USED edge.
        REQUIRE(vertices.count(vv) != 0);
      }
      vertex_mapping[vv];
      ++num_tokens;
    }
    REQUIRE(num_tokens == vertex_mapping.size());
  } else {
    REQUIRE(index == str.size());
    for (auto vv : vertices) {
      vertex_mapping[vv];
    }
  }
  // NOW, perform the swaps.
  get_problem_mapping(vertex_mapping, swaps);
}

DecodedArchitectureData::DecodedArchitectureData() : number_of_vertices(0) {}

DecodedArchitectureData::DecodedArchitectureData(
    const std::string& solution_edges_string) {
  vector<vector<std::size_t>> neighbours(1);
  std::set<std::size_t> vertices_seen;
  for (unsigned char ch : solution_edges_string) {
    if (ch != ':') {
      const auto new_v = char_to_vertex_map().at(ch);
      neighbours.back().push_back(new_v);
      vertices_seen.insert(new_v);
      continue;
    }
    // We move onto the next vertex.
    neighbours.emplace_back();
  }
  // The last vertex N cannot have any neighbours j with j>N,
  // so we don't bother to record it in the string,
  // so it's not stored in "neighbours".
  number_of_vertices = neighbours.size() + 1;
  CHECK(number_of_vertices >= 4);
  // But everything MUST be joined to something, if the graph is connected.
  // Vertex v won't be listed if it only joins higher-numbered vertices,
  // so many vertices might not be mentioned here.
  REQUIRE(!vertices_seen.empty());
  REQUIRE(*vertices_seen.crbegin() <= neighbours.size());

  for (std::size_t ii = 0; ii < neighbours.size(); ++ii) {
    if (neighbours[ii].empty()) {
      continue;
    }
    REQUIRE(std::is_sorted(neighbours[ii].cbegin(), neighbours[ii].cend()));
    REQUIRE(neighbours[ii][0] > ii);
    REQUIRE(
        std::adjacent_find(neighbours[ii].cbegin(), neighbours[ii].cend()) ==
        neighbours[ii].cend());
    for (auto jj : neighbours[ii]) {
      edges.insert(get_swap(ii, jj));
    }
  }
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
