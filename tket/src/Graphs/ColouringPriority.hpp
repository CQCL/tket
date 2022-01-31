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

#include <cstddef>
#include <set>
#include <string>
#include <vector>

namespace tket {
namespace graphs {

class AdjacencyData;

// TODO: if we only allow 64 colours
// (which follows if, e.g., we only have <= 64 vertices),
// then it can be done faster with bit operations.
// 64^64 > 10^115 is ridiculously big, so it's not very likely
// we'll get a real graph, which we have any realistic
// hope of colouring well, needing more than 64 colours.
// For C colours, a potential speedup of roughly a factor C.
//
// If we are going to find the minimum possible number of colours
// to colour a graph, it is exponential time in the worst case. However,
// it still might be possible to speed it up substantially from brute force,
// if we colour the vertices in a particular order.
// As we descend along the sequence of colourings, we will hit
// colouring inconsistencies and have to backtrack.
// But if we have selected a good order, we might hit inconsistencies
// quite early on average, and hence prune a lot of the tree
// of all possibilities, so we hopefully examine only a fraction
// of all candidates, yet still guarantee to find an optimal colouring.
//
// (In fact, this clearly is happening, since we are able to colour,
// apparently optimally, some fairly dense graphs on 30 vertices in << 1 second
// with 3 colours. A trivial brute force colouring would need 3^30 ~ 2*10^14
// operations, so it is many orders of magnitude faster than that).
//

/**
 * Use some simple heuristics (quite fast) to try to find a good vertex order
 * for our attempted colouring, which will then be done by brute force (the slow
 * part). A good ordering would detect colour inconsistencies rapidly and so not
 * do too much work before backtracking.
 */
class ColouringPriority {
 public:
  /**
   * Only calculate within a single connected component.
   * @param adjacency_data Data for the whole graph.
   * @param vertices_in_component The vertices to consider;
   * these must form a single component. The caller should already
   * have calculated the connected components.
   * @param initial_clique used as a seed (although not actually checked to be a
   * clique: it could instead simply be a set of vertices with many edges). For
   * a greedy-like priority algorithm we colour them first, all different
   * colours.
   */
  ColouringPriority(
      const AdjacencyData& adjacency_data,
      const std::set<std::size_t>& vertices_in_component,
      const std::set<std::size_t>& initial_clique);

  /**
   * Contains extra data about a vertex. These node objects will be put into a
   * vector which defines the coluring order.
   */
  struct NodeData {
    /** The original vertex id (index). */
    std::size_t vertex;

    /**
     * Gives the indices in the vector of nodes (NOT the vertex IDs)
     * of all earlier nodes with a vertex joined to this vertex.
     * Thus, this node must be coloured differently from any of those nodes.
     * Of course, LATER nodes can also be joined to this,
     * so this does not give ALL neighbours.
     */
    std::vector<std::size_t> earlier_neighbour_node_indices;
  };

  /**
   * The list of nodes, giving the order in which to colour the nodes.
   */
  typedef std::vector<NodeData> Nodes;

  /**
   * Get the colouring order. Was already calculated upon construction.
   */
  const Nodes& get_nodes() const;

  /**
   * Helpful for some algorithms to store and return the original
   * "clique" which was passed in.
   */
  const std::set<std::size_t>& get_initial_clique() const;

  /**
   * For debugging, it's helpful to be able to copy/paste
   * nasty examples into tests. This prints out C++ code
   * representing the neighbour data.
   */
  std::string print_raw_data(bool relabel_to_simplify = true) const;

 private:
  // Store a copy of the initial clique passed in.
  const std::set<std::size_t> m_initial_clique;
  Nodes m_nodes;
};

}  // namespace graphs
}  // namespace tket
