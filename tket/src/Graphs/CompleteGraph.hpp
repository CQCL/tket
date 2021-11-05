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

#ifndef _TKET_CompleteGraph_H
#define _TKET_CompleteGraph_H

#include <vector>

#include "AbstractGraph.hpp"

namespace tket::graphs {

template <typename T>
class CompleteGraph : public AbstractGraph<T> {
 public:
  using AbstractGraph<T>::node_exists;
  using AbstractGraph<T>::edge_exists;
  using AbstractGraph<T>::get_all_nodes_set;
  using AbstractGraph<T>::get_all_nodes_vec;
  using AbstractGraph<T>::get_all_edges_vec;
  using typename AbstractGraph<T>::Edge;

  CompleteGraph() {}

  CompleteGraph(unsigned n) {
    for (unsigned i = 0; i < n; i++) {
      nodes_.insert(T(i));
    }
  }

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

  void add_node(const T& node) { nodes_.insert(node); }

  bool operator==(const CompleteGraph<T>& other) const {
    return nodes_ == other.nodes_;
  }

 protected:
  using AbstractGraph<T>::nodes_;

  bool edge_exists(const T& node1, const T& node2) const override {
    if (!node_exists(node1) || !node_exists(node2)) {
      throw NodeDoesNotExistError(
          "The UIDs passed to CompleteGraph::edge_exists must exist.");
    }
    return true;
  }
};

}  // namespace tket::graphs

#endif  // _TKET_CompleteGraph_H
