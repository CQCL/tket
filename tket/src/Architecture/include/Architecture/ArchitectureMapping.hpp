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

#include "Architecture/Architecture.hpp"
#include "TokenSwapping/SwapFunctions.hpp"

namespace tket {

/** Intended for use with TokenSwapping.
 *  For mapping between nodes in an architecture and size_t vertex numbers.
 *  The vertex numbers are merely the indices of each Node
 *  within the vector returned by the get_all_nodes() function.
 *
 *  For now, we don't want to use Node objects as (1) this would make
 *  TokenSwapping dependent on other parts of Tket and hence less modular,
 *  (2) it would probably slow things down significantly because Nodes
 *  contain extra data, like vectors and strings, which are relatively
 *  expensive to copy; vertices get copied and moved around many times
 *  by any TSA.
 *
 *  TODO: it would be better to use a Vertex wrapper class
 *  instead of raw size_t. (Also, might change to unsigned instead of size_t).
 */
class ArchitectureMapping {
 public:
  /** The object must remain valid and unchanged
   *  throughout the life of this object.
   *  @param arch The finished Architecture object, must remain valid
   *    for the lifetime of this object.
   */
  explicit ArchitectureMapping(const Architecture& arch);

  /** If the architecture object was initialised with explicit edges,
   * use these edges (rather than the Architecture nodes() function)
   * to create the Node <-> size_t mapping, in a fixed way not dependent
   * on Architecture (the reason being that Architecture does not guarantee
   * the mapping, but if we change the labels then we change to an isomorphic
   * but different token swapping problem, which messes up testing.
   * In practice every implementation of token swapping, except for the ultimate
   * probably exponential-time optimal algorithm, is going to depend
   * on the labels. Even if we had a fast graph isomorphism routine, the labels
   * would still not be uniquely determined, as they could be permuted).
   * @param arch The finished Architecture object, must remain valid
   *    for the lifetime of this object.
   * @param edges Edges originally used to construct the Architecture object.
   *    These will uniquely determine the internal Node <-> size_t mapping.
   */
  ArchitectureMapping(
      const Architecture& arch,
      const std::vector<std::pair<unsigned, unsigned>>& edges);

  /** Convenient reference to the Architecture object we used
   *  to construct this ArchitectureMapping.
   */
  const Architecture& get_architecture() const;

  /** The number of vertices in the Architecture.
   *  @return The number of vertices
   */
  size_t number_of_vertices() const;

  /** Get the newly created vertex assigned to the node.
   *  Throws if the node is invalid.
   *  @param node The node within the original Architecture object
   *  @return The newly created vertex representing this node
   */
  size_t get_vertex(const Node& node) const;

  /** Reverse of "get_vertex", throws if the vertex is invalid.
   *  @param vertex The vertex created by this ArchitectureMapping object.
   *  @return The node corresponding to this vertex.
   */
  const Node& get_node(size_t vertex) const;

  /** Get the edges using the vertices created by this ArchitectureMapping
   *  object. The vertex numbers, of course, do not necessarily match with
   *  the Node uids of the underlying architecture object
   *  (that's why we have a mapping).
   *  @return The vector of edges in the architecture, using the new
   *    vertex numbers.
   */
  std::vector<Swap> get_edges() const;

 private:
  /// Store a reference to the Architecture passed into the constructor.
  const Architecture& m_arch;

  /// Element i is simply the node corresponding to vertex i.
  node_vector_t m_vertex_to_node_mapping;

  /// Reverse of m_vertex_to_node_mapping; look up the index of a node.
  std::map<Node, size_t> m_node_to_vertex_mapping;
};

}  // namespace tket
