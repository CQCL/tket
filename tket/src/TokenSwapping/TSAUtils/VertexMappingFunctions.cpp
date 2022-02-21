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

#include "TokenSwapping/VertexMappingFunctions.hpp"

#include <sstream>
#include <stdexcept>

#include "TokenSwapping/VertexSwapResult.hpp"
#include "Utils/Assert.hpp"

namespace tket {

using namespace tsa_internal;

bool all_tokens_home(const VertexMapping& vertex_mapping) {
  for (const auto& entry : vertex_mapping) {
    if (entry.first != entry.second) {
      return false;
    }
  }
  return true;
}

void check_mapping(
    const VertexMapping& vertex_mapping, VertexMapping& work_mapping) {
  work_mapping.clear();
  for (const auto& entry : vertex_mapping) {
    // GCOVR_EXCL_START
    TKET_ASSERT(
        work_mapping.count(entry.second) == 0 ||
        AssertMessage() << "Vertices v_" << entry.first << " and v_"
                        << work_mapping[entry.second]
                        << " both have the same target vertex v_"
                        << entry.second);
    // GCOVR_EXCL_STOP
    work_mapping[entry.second] = entry.first;
  }
}

void check_mapping(const VertexMapping& vertex_mapping) {
  VertexMapping work_mapping;
  check_mapping(vertex_mapping, work_mapping);
}

void append_swaps_to_interchange_path_ends(
    const std::vector<size_t>& path, VertexMapping& vertex_mapping,
    SwapList& swap_list) {
  if (path.size() < 2 || path.front() == path.back()) {
    return;
  }
  for (size_t ii = path.size() - 1; ii > 0; --ii) {
    VertexSwapResult(path[ii], path[ii - 1], vertex_mapping, swap_list);
  }
  for (size_t ii = 2; ii < path.size(); ++ii) {
    VertexSwapResult(path[ii], path[ii - 1], vertex_mapping, swap_list);
  }
}

size_t get_source_vertex(
    VertexMapping& source_to_target_map, size_t target_vertex) {
  if (source_to_target_map.count(target_vertex) == 0) {
    // If it IS a genuine permutation mapping (which we assume),
    // then the vertex is as yet unmentioned (and hence unmoved).
    source_to_target_map[target_vertex] = target_vertex;
    return target_vertex;
  }
  for (const auto& entry : source_to_target_map) {
    if (entry.second == target_vertex) {
      return entry.first;
    }
  }
  TKET_ASSERT(!"get_source_vertex");
  return target_vertex;
}

void add_swap(VertexMapping& source_to_target_map, const Swap& swap) {
  const auto source_v1 = get_source_vertex(source_to_target_map, swap.first);
  const auto source_v2 = get_source_vertex(source_to_target_map, swap.second);
  std::swap(source_to_target_map[source_v1], source_to_target_map[source_v2]);
}

}  // namespace tket
