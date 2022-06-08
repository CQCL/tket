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

#include <vector>

#include "AbstractGraph.hpp"

namespace tket::graphs {

template <typename T>
class CompleteGraph : public AbstractGraph<T> {
 public:
  using AbstractGraph<T>::node_exists;
  using AbstractGraph<T>::edge_exists;
  using AbstractGraph<T>::get_all_nodes_vec;
  using AbstractGraph<T>::get_all_edges_vec;
  using AbstractGraph<T>::n_nodes;
  using typename AbstractGraph<T>::Edge;

  /** Construct an empty graph. */
  CompleteGraph() {}

  /** All edges as a vector. */
  std::vector<Edge> get_all_edges_vec() const override {
    std::vector<Edge> edges;
    std::vector<T> nodes = get_all_nodes_vec();
    unsigned n = nodes.size();
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = i + 1; j < n; j++) {
        edges.push_back({nodes[i], nodes[j]});
      }
    }
    return edges;
  }

  /** Add a new node to the graph. */
  void add_node(const T& node) { nodes_.insert(node); }

  /** Test whether two nodes are connected. */
  bool edge_exists(const T& node1, const T& node2) const override {
    if (!node_exists(node1) || !node_exists(node2)) {
      throw NodeDoesNotExistError(
          "The UIDs passed to CompleteGraph::edge_exists must exist.");
    }
    return true;
  }

  unsigned get_distance(const T& node1, const T& node2) const override {
    if (!node_exists(node1) || !node_exists(node2)) {
      throw NodeDoesNotExistError(
          "The UIDs passed to CompleteGraph::edge_exists must exist.");
    }
    if (node1 == node2) return 0;
    return 1;
  }

  unsigned get_diameter() override {
    switch (n_nodes()) {
      case 0:
        throw std::logic_error("Graph is empty.");
      case 1:
        return 0;
      default:
        return 1;
    }
  }

  bool operator==(const CompleteGraph<T>& other) const {
    return nodes_ == other.nodes_;
  }

 protected:
  using AbstractGraph<T>::nodes_;
};

}  // namespace tket::graphs
