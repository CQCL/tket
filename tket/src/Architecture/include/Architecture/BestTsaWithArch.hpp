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

#include "ArchitectureMapping.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {

/** A simple wrapper around BestFullTsa from TokenSwapping,
 * using Architecture objects directly to find distances and neighbours.
 */
struct BestTsaWithArch {
  /** Given the desired vertex mapping, a list
   * of swaps (which may or may not be empty), and information about
   * the architecture (the underlying graph), append extra swaps to it
   * to produce the desired mapping.
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping. Will be updated with
   * the new added swaps.
   *  @param arch_mapping An ArchitectureMapping object, which knows the graph,
   * and how to do Node <-> vertex size_t conversions.
   */
  static void append_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      const ArchitectureMapping& arch_mapping);

  /** This specifies desired source->target vertex mappings.
   *  Any nodes not occurring as a key might be moved by the algorithm.
   */
  typedef std::map<Node, Node> NodeMapping;

  /** Given an architecture and desired source->target node mapping,
   * compute a sequence of swaps (attempts to be as short as possible)
   * which will perform that mapping.
   * Note that it may use ALL the nodes in the architecture,
   * not just the ones occurring in the node_mapping.
   * If you wish certain nodes to be fixed, specify them in the mapping
   * (with equal source and target).
   * (However, note that they might STILL be moved, as long as by the end
   * they are back at the start. If you really don't to involve a particular
   * node, you must remove it completely from the architecture).
   * KNOWN BUG: it may give an error with disconnected architectures.
   *  @param architecture The raw object containing the graph.
   *  @param node_mapping The desired source->target node mapping.
   *  @return The required list of node pairs to swap.
   */
  static std::vector<std::pair<Node, Node>> get_swaps(
      const Architecture& architecture, const NodeMapping& node_mapping);
};

}  // namespace tket
