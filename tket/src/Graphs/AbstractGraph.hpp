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

#ifndef _TKET_AbstractGraph_H
#define _TKET_AbstractGraph_H

#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tket::graphs {

class NodeDoesNotExistError : public std::logic_error {
  using std::logic_error::logic_error;
};

class EdgeDoesNotExistError : public std::logic_error {
  using std::logic_error::logic_error;
};

/**
 * Abstract class for representing graphs.
 *
 * @tparam T type of nodes in the graph
 */
template <typename T>
class AbstractGraph {
 protected:
  using Edge = std::pair<T, T>;
  std::set<T> nodes_;

 public:
  /** Construct an empty graph */
  AbstractGraph() : nodes_() {}

  /** Construct from vector of nodes */
  explicit AbstractGraph(const std::vector<T> &nodes)
      : nodes_({nodes.begin(), nodes.end()}) {}

  /** Check if an edge exists */
  virtual bool edge_exists(const T &node1, const T &node2) const = 0;

  /** Check if a node exists */
  bool node_exists(const T &node) const { return nodes_.contains(node); }

  const std::set<T> &get_all_nodes() const { return nodes_; }

  std::set<T> get_all_nodes_set() const { return nodes_; }

  std::vector<T> get_all_nodes_vec() const {
    return {nodes_.begin(), nodes_.end()};
  }

  /** Number of nodes */
  unsigned n_nodes() const { return nodes_.size(); }

  virtual std::vector<Edge> get_all_edges_vec() const = 0;
};

}  // namespace tket::graphs

#endif  // _TKET_AbstractGraph_H
