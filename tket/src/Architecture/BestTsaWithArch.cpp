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

#include "BestTsaWithArch.hpp"

#include "DistancesFromArchitecture.hpp"
#include "NeighboursFromArchitecture.hpp"
#include "TokenSwapping/BestFullTsa.hpp"
#include "Utils/Assert.hpp"
#include "Utils/RNG.hpp"

namespace tket {

using namespace tsa_internal;

void BestTsaWithArch::append_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    const ArchitectureMapping& arch_mapping) {
  DistancesFromArchitecture distances(arch_mapping);
  NeighboursFromArchitecture neighbours(arch_mapping);
  RNG rng;
  RiverFlowPathFinder path_finder(distances, neighbours, rng);
  BestFullTsa().append_partial_solution(
      swaps, vertex_mapping, distances, neighbours, path_finder);
}

std::vector<std::pair<Node, Node>> BestTsaWithArch::get_swaps(
    const Architecture& architecture, const NodeMapping& node_mapping) {
  std::vector<std::pair<Node, Node>> swaps;
  // Before all the conversion and object construction,
  // doesn't take long to check if it's actually trivial
  bool trivial = true;
  for (const auto& entry : node_mapping) {
    if (entry.first != entry.second) {
      trivial = false;
      break;
    }
  }
  if (trivial) {
    return swaps;
  }
  // Now convert the Nodes into raw vertices for use in TSA objects.
  const ArchitectureMapping arch_mapping(architecture);
  VertexMapping vertex_mapping;
  for (const auto& node_entry : node_mapping) {
    vertex_mapping[arch_mapping.get_vertex(node_entry.first)] =
        arch_mapping.get_vertex(node_entry.second);
  }
  TKET_ASSERT(vertex_mapping.size() == node_mapping.size());
  check_mapping(vertex_mapping);

  SwapList raw_swap_list;
  BestTsaWithArch::append_solution(raw_swap_list, vertex_mapping, arch_mapping);

  // Finally, convert the raw swaps back to nodes.
  swaps.reserve(raw_swap_list.size());
  for (auto id_opt = raw_swap_list.front_id(); id_opt;
       id_opt = raw_swap_list.next(id_opt.value())) {
    const auto& raw_swap = raw_swap_list.at(id_opt.value());
    swaps.emplace_back(std::make_pair(
        arch_mapping.get_node(raw_swap.first),
        arch_mapping.get_node(raw_swap.second)));
  }
  return swaps;
}

}  // namespace tket
