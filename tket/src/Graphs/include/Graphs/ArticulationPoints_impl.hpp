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

#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/property.hpp>

#include "Graphs/DirectedGraph.hpp"
#include "Utils/UnitID.hpp"

namespace tket::graphs {

template <typename T>
using UndirectedConnGraph = typename DirectedGraph<T>::UndirectedConnGraph;

// ::detail includes all internal implementation-specific details
namespace detail {

/**
 *  Exception raised when subarchitecture is empty
 */
class NoSelectedComponent : public std::logic_error {
 public:
  using std::logic_error::logic_error;
};

/**
 *  A BicomponentGraph has biconnected components as vertices
 *  and APs as edges: two vertices are connected iff there is an AP
 *  that connects them
 *
 *  The biconnected component graph is in fact a tree: if there was a
 *  cycle, then one of the APs could be removed without losing connectivity
 */
template <typename T>
class BicomponentGraph {
 private:
  struct BicomponentGraphEdge {
    T ap;
  };
  using Graph = boost::adjacency_list<
      boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
      BicomponentGraphEdge>;
  // Map from vertices v (T) to the set of biconnected components of v
  using VertexToCompT = std::map<T, std::set<unsigned>>;

 public:
  /** Construct a biconnected component graph from the graph */
  explicit BicomponentGraph(const UndirectedConnGraph<T>& graph);

  /** Selects all components that contain vertices in `[begin, end)` */
  template <typename Range>
  void select_comps(Range nodes);

  /**
   * minimally expands the list of selected components so that the selected
   * components form a connected subgraph (this is well defined since the
   * graph is a tree)
   */
  void propagate_selected_comps();

  /** find all interior edges of the selected components subgraph
   *  i.e. edges that connect two selected components.
   *
   *  These are the subgraph APs
   */
  std::set<T> get_inner_edges();

 private:
  // the bicomponent graph
  Graph g;
  // the original graph (passed as argument to the constructor)
  const UndirectedConnGraph<T>& underlying_g;
  // binary flags of selected components
  std::vector<bool> selected_comps;
  // list of APs
  std::vector<T> aps;
  // map from vertices to set of components they belong to
  VertexToCompT vertex_to_comps;

  /** from the underlying graph, compute APs and the map from vertex to
   * components
   *  this initialises `selected_comps`, `aps`, `vertex_to_comps` and
   * `n_components` */
  void compute_components_map();

  /** using the values computed in `compute_components_map`, build the
   * BicomponentGraph this is called by the constructor and initialises `g` */
  void build_graph();

  // for testing
  friend class BicomponentGraphTester;
};

extern template class BicomponentGraph<UnitID>;
extern template class BicomponentGraph<Node>;
extern template class BicomponentGraph<Qubit>;
extern template void BicomponentGraph<UnitID>::select_comps(
    std::vector<UnitID>);
extern template void BicomponentGraph<Node>::select_comps(std::vector<Node>);
extern template void BicomponentGraph<Qubit>::select_comps(std::vector<Qubit>);

/** Class exposing private API of BicomponentGraph for testing
 */
class BicomponentGraphTester {
 private:
  using T = Node;
  detail::BicomponentGraph<T>* bicomp_g;

 public:
  explicit BicomponentGraphTester(detail::BicomponentGraph<T>* bg)
      : bicomp_g(bg) {}
  const std::vector<bool>& get_selected_comps() const {
    return bicomp_g->selected_comps;
  }
  const std::set<unsigned>& get_comps(T node) const {
    return bicomp_g->vertex_to_comps[node];
  }
  const detail::BicomponentGraph<Node>::Graph& get_graph() const {
    return bicomp_g->g;
  }
  unsigned n_components() const { return bicomp_g->selected_comps.size(); }
};

}  // namespace detail

}  // namespace tket::graphs
