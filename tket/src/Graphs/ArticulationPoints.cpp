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

#include "Graphs/ArticulationPoints.hpp"

#include <boost/property_map/property_map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "Graphs/Utils.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket::graphs {

namespace detail {

// Types used in `get_ap_to_comp`
// Maps BGL edges descriptors to their (unique) biconnected component
template <typename T>
using EdgeToCompT = std::map<utils::edge<UndirectedConnGraph<T>>, unsigned>;

template <typename T>
BicomponentGraph<T>::BicomponentGraph(const UndirectedConnGraph<T>& graph)
    : underlying_g(graph) {
  // get map from APs to their bicomponents (as well as the number of bicomps)
  compute_components_map();

  // all that remains is to build the graph
  build_graph();
}

/** Computes map from APs (identified by UnitID) to the set of
 *  biconnected components they belong to.
 *  Returns a pair with the map and the number of components
 */
template <typename T>
void BicomponentGraph<T>::compute_components_map() {
  // Map from edges to their (unique) biconnected component
  // (this is filled out by call to boost::biconnected_components)
  EdgeToCompT<T> edge_to_comp;

  std::vector<unsigned> underlying_aps;
  // call to BGL
  auto bgl_res = boost::biconnected_components(
      underlying_g, boost::make_assoc_property_map(edge_to_comp),
      std::back_inserter(underlying_aps));
  for (auto v_ap : underlying_aps) {
    aps.push_back(underlying_g[v_ap]);
  }
  unsigned n_components = bgl_res.first;

  // Now populate map from vertices to the set of biconnected components of v
  for (auto v : boost::make_iterator_range(boost::vertices(underlying_g))) {
    T node = underlying_g[v];
    vertex_to_comps.insert({node, {}});
    // we get the components `v` belongs to by looking at the components
    // of its incident edges
    if (boost::degree(v, underlying_g) > 0) {
      for (auto e :
           boost::make_iterator_range(boost::out_edges(v, underlying_g))) {
        // it's a set, we do not need to worry about duplicates
        vertex_to_comps[node].insert(edge_to_comp[e]);
      }
    } else {
      // if vertex is disconnected, create new exclusive component
      vertex_to_comps[node].insert(n_components++);
    }
  }
  selected_comps = std::vector<bool>(n_components, false);
}

template <typename T>
void BicomponentGraph<T>::build_graph() {
  // initialise bicomponent graph with `n_components` vertices
  unsigned n_components = selected_comps.size();
  g = Graph(n_components);

  // add edges according to vertex_to_comps
  for (auto ap : aps) {
    auto comps = vertex_to_comps[ap];
    for (auto it = comps.cbegin(); it != comps.cend(); ++it) {
      unsigned comp1 = *it;
      for (auto it2 = it; ++it2 != comps.cend();) {
        unsigned comp2 = *it2;
        // add edge between any 2 components that `ap` links
        boost::add_edge(comp1, comp2, {ap}, g);
      }
    }
  }
}

template <typename T>
template <typename Range>
void BicomponentGraph<T>::select_comps(Range nodes) {
  for (T node : nodes) {
    for (unsigned comp : vertex_to_comps[node]) {
      selected_comps[comp] = true;
    }
  }
}

/** Propagation of selected components (`propagate_selected_comps`)
 *
 *  Strategy to minimally expand the selected components to a connected
 * subgraph: given that the bicomp graph is a tree, we simply need to select
 * every vertex that is on a path between two selected vertices
 *
 *  We can easily achieve that using e.g. depth-first-search (DFS):
 *     i) we fix a vertex that is selected as root
 *    ii) at each vertex we discover, we check recursively if one of the
 * descendants is selected, and if so, we select the vertex
 */

/** A BGL custom visitor, in all its glory
 *  It checks if after visiting all the descendants, one of them
 *  was found to be selected.
 *  If so, then vertex must be on a path between two selected vertices
 *  and thus the vertex is selected
 */
template <typename Graph>
struct TrackUsedEdgesVisitor : public boost::default_dfs_visitor {
  using edge =
      std::pair<graphs::utils::vertex<Graph>, graphs::utils::vertex<Graph>>;

  explicit TrackUsedEdgesVisitor(std::vector<bool>& s) : selected(s) {
    tree_edges = std::make_shared<std::set<edge>>();
  }
  std::vector<bool>& selected;
  std::shared_ptr<std::set<edge>> tree_edges;

  // called as DFS leaves the edge and returns to parent
  void finish_edge(graphs::utils::edge<Graph> e, const Graph& g) {
    unsigned target = boost::target(e, g);
    unsigned source = boost::source(e, g);
    // only consider tree edges
    if (!tree_edges->count({source, target})) {
      return;
    }
    if (selected[target]) {
      selected[source] = true;
    }
  }
  // keep track of tree edges
  void tree_edge(graphs::utils::edge<Graph> e, const Graph& g) {
    unsigned target = boost::target(e, g);
    unsigned source = boost::source(e, g);
    tree_edges->insert({source, target});
  }
};

template <typename T>
void BicomponentGraph<T>::propagate_selected_comps() {
  unsigned root = 0;
  unsigned n_components = selected_comps.size();
  while (root < n_components && !selected_comps[root]) {
    ++root;
  }

  if (root == n_components) {
    throw NoSelectedComponent(
        "At least one component must be selected"
        " to be able to propagate");
  }

  // using our custom visitor, DFS will propagate selected_comps as necessary
  auto params = boost::visitor(TrackUsedEdgesVisitor<Graph>(selected_comps))
                    .root_vertex(root);
  boost::depth_first_search(g, params);
}

template <typename T>
std::set<T> BicomponentGraph<T>::get_inner_edges() {
  std::set<T> out;
  for (auto e : boost::make_iterator_range(boost::edges(g))) {
    unsigned comp1 = boost::source(e, g);
    unsigned comp2 = boost::target(e, g);
    if (selected_comps[comp1] && selected_comps[comp2]) {
      out.insert(g[e].ap);
    }
  }
  return out;
}

template class BicomponentGraph<UnitID>;
template class BicomponentGraph<Node>;
template class BicomponentGraph<Qubit>;
template void BicomponentGraph<UnitID>::select_comps(std::vector<UnitID>);
template void BicomponentGraph<Node>::select_comps(std::vector<Node>);
template void BicomponentGraph<Qubit>::select_comps(std::vector<Qubit>);

}  // namespace detail

template <typename T>
std::set<T> get_subgraph_aps(
    const UndirectedConnGraph<T>& graph,
    const UndirectedConnGraph<T>& subgraph) {
  detail::BicomponentGraph<T> bicomp_graph(graph);
  using funcT = std::function<T(unsigned)>;
  funcT to_node = [&subgraph](unsigned v) { return subgraph[v]; };
  auto selected = boost::make_iterator_range(boost::vertices(subgraph)) |
                  boost::adaptors::transformed(to_node);
  bicomp_graph.select_comps(selected);
  bicomp_graph.propagate_selected_comps();
  return bicomp_graph.get_inner_edges();
}

template std::set<UnitID> get_subgraph_aps(
    const UndirectedConnGraph<UnitID>& graph,
    const UndirectedConnGraph<UnitID>& subgraph);
template std::set<Node> get_subgraph_aps(
    const UndirectedConnGraph<Node>& graph,
    const UndirectedConnGraph<Node>& subgraph);

}  // namespace tket::graphs
