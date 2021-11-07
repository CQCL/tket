// Copyright 2019-2021 Cambridge Quantum Computing
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

#ifndef _TKET_Architecture_H_
#define _TKET_Architecture_H_

#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Graphs/CompleteGraph.hpp"
#include "Graphs/DirectedGraph.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Json.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/TketLog.hpp"
#include "Utils/UnitID.hpp"
namespace tket {

using dist_vec = graphs::dist_vec;

class ArchitectureInvalidity : public std::logic_error {
 public:
  explicit ArchitectureInvalidity(const std::string &message)
      : std::logic_error(message) {}
};

static std::vector<std::pair<Node, Node>> as_nodepairs(
    const std::vector<std::pair<unsigned, unsigned>> &edges) {
  std::vector<std::pair<Node, Node>> nodepairs;
  nodepairs.reserve(edges.size());
  for (auto &[m, n] : edges) {
    nodepairs.push_back({Node(m), Node(n)});
  }
  return nodepairs;
}

/**
 * Generic architecture class
 *
 * Constraint: the underlying graph type must have node type @ref Node.
 *
 * @tparam T underlying graph type
 */
template <typename T>
class ArchitectureBase : public T {
 public:
  using T::edge_exists;
  using T::get_all_edges_vec;
  using T::nodes;
  using T::get_all_nodes_vec;
  using T::n_nodes;
  using T::node_exists;
  using T::T;

  /**
   * Construct from a vector of pairs of indices in the default register.
   */
  explicit ArchitectureBase(
      const std::vector<std::pair<unsigned, unsigned>> &edges)
      : ArchitectureBase(as_nodepairs(edges)) {}

  /** compute articulation points of Architecture
   *  if `subarc` is given, the articulation that do not affect the
   * sub-architecture are ignored */
  node_set_t get_articulation_points() const;
  node_set_t get_articulation_points(const ArchitectureBase &subarc) const;

  /* returns new Architecture that is generated by a subset of nodes of `this`
   */
  ArchitectureBase create_subarch(const std::vector<Node> &nodes);

  // Returns vectors of nodes which correspond to lines of specified length
  std::vector<node_vector_t> get_lines(
      std::vector<unsigned> required_lengths) const;

  // get architecture diameter
  unsigned get_diameter() const;

  /** Removes 'num' nodes from architecture, with 'worseness' of
   * nodes determined by a heuristic. */
  node_set_t remove_worst_nodes(unsigned num);

  /**
   * gives the connectivity matrix of the architecture
   * @return connectivity matrix
   */
  MatrixXb get_connectivity() const;

 protected:
  // Returns node with least connectivity given some distance matrix.
  std::optional<Node> find_worst_node(const ArchitectureBase &orig_g);
};

typedef ArchitectureBase<graphs::DirectedGraph<Node>> Architecture;

JSON_DECL(Architecture::Connection)
JSON_DECL(Architecture)

class FullyConnected : public ArchitectureBase<graphs::CompleteGraph<Node>> {
 public:
  FullyConnected() : ArchitectureBase<graphs::CompleteGraph<Node>>() {}
  FullyConnected(unsigned n) {
    for (unsigned i = 0; i < n; i++) {
      nodes_.insert(Node("fcNode", i));
    }
  }
};

JSON_DECL(FullyConnected)

// Subclass, constructor generates adjacency matrix corresponding to a ring
// architecture
class RingArch : public Architecture {
 public:
  explicit RingArch(unsigned numberOfNodes);
  // nodes() does not guarantee to return nodes in any order
  // this returns the canonical ordering of nodes

 private:
  static std::vector<Connection> get_edges(unsigned numberOfNodes);
};

// Subclass, constructor generates adjacency matrix corresponding to a
// SquareGrid architecture
class SquareGrid : public Architecture {
 public:
  // Converts square indexing to qubit indexing
  Vertex squind_to_qind(
      const unsigned ver, const unsigned hor, const unsigned layer = 0) const {
    return (ver * dimension_c + hor) + single_layer_nodes() * layer;
  }
  // Returns number of nodes in a single 2d layer
  unsigned single_layer_nodes() const { return dimension_c * dimension_r; }
  // Gets number of columns of square grid architecture
  unsigned get_columns() const { return dimension_c; }
  // Gets number of rows of square grid architecture
  unsigned get_rows() const { return dimension_r; }
  // Gets number of layers of square grid architecture
  unsigned get_layers() const { return layers; }
  // Converts qubit indexing to square indexing
  std::pair<unsigned, unsigned> qind_to_squind(const Vertex &qn) const {
    unsigned col = qn % dimension_c;
    unsigned row = (qn - col) / dimension_c;
    return {row, col};
  }
  // SquareGrid constructor
  // dim_c equiv 'x', dim_r equiv 'y'
  SquareGrid(
      const unsigned dim_r, const unsigned dim_c, const unsigned _layers = 1);

 private:
  static std::vector<Connection> get_edges(
      const unsigned dim_r, const unsigned dim_c, const unsigned layers = 1);

  unsigned dimension_r;
  unsigned dimension_c;
  unsigned layers;
};

int tri_lexicographical_comparison(
    const dist_vec &dist1, const dist_vec &dist2);

}  // namespace tket

#endif
