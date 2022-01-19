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

#include <algorithm>
#include <boost/graph/copy.hpp>
#include <set>
#include <type_traits>

#include "Graphs/Utils_impl.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket::graphs {

using dist_vec = std::vector<std::size_t>;

namespace utils {

// Type synonyms for many often-used BGL types

// BGL vertex descriptor
template <typename T>
using vertex = detail::vertex<T>;
// BGL edge descriptor
template <typename T>
using edge = detail::edge<T>;

/* Compile-time boolean flag to check if boost graph type is
 * directed or undirected graph
 *
 * Template Argument Graph is the boost type of the graph to be checked
 */
template <typename Graph>
constexpr bool is_directed = detail::is_directed_struct<Graph>::value;

/* These utility functions are meant to offer a more consistent
 * alternative to boost::remove_vertex, so that the vertex storage container
 * can be changed without changes in the implementation
 *
 * boost::remove_vertex behaves differently depending on the vertex container
 * type
 *  - for vecS, remove_vertex will shift vertex indices so that they remain in
 * the continuous interval [0..n-1].
 *  - for other types, remove_vertex will leave the vertex descriptors unchanged
 *
 * Instead, these wrappers behave as follows
 *  - support for explicit indices: if vertex descriptors are numbers
 *    (i.e. VertexList == vecS) or if an explicit index map is provided, indices
 * will be shifted as with the default boost vecS implementation.
 *  - support for external map: additionally, a map from vertex descriptors to
 *    arbitrary types can be provided. The map keys will then automatically be
 * changed to mirror the vertex descriptors changes.
 *
 * As a result, implementations can abstract away the vertex descriptor
 * implementation, and either rely on explicit indices throughout, or stick to
 * vertex descriptors and have the necessary translation maps be updated
 * automatically
 */
template <typename Graph>
inline void remove_vertex(vertex<Graph> v, Graph& graph) {
  detail::graph_utils_impl<Graph> impl(graph);
  impl.remove_vertex(v);
}
template <typename Graph, typename PMap>
inline void remove_vertex(vertex<Graph> v, Graph& graph, PMap& p) {
  detail::graph_utils_impl<Graph, PMap> impl(graph, p);
  impl.remove_vertex(v);
}
template <typename Graph, typename PMap, typename Map>
inline void remove_vertex(vertex<Graph> v, Graph& graph, PMap& pmap, Map& map) {
  detail::graph_utils_impl_with_map<Graph, Map, PMap> impl(graph, map, pmap);
  impl.remove_vertex(v);
}
template <typename Graph, typename Map, typename PMap>
inline void remove_vertex_with_map(
    vertex<Graph> v, Graph& graph, Map& map, PMap& pmap) {
  detail::graph_utils_impl_with_map<Graph, Map, PMap> impl(graph, map, pmap);
  impl.remove_vertex(v);
}
template <typename Graph, typename Map>
inline void remove_vertex_with_map(vertex<Graph> v, Graph& graph, Map& map) {
  detail::graph_utils_impl_with_map<Graph, Map> impl(graph, map);
  impl.remove_vertex(v);
}

// like boost::remove_edge, but optionally removes disconnected vertices
// in which case it can also update a map as above in graph_utils::remove_vertex
template <typename Graph>
inline void remove_edge(
    edge<Graph> e, Graph& graph, bool remove_unused_vertices = false) {
  detail::graph_utils_impl<Graph> impl(graph);
  impl.remove_edge(e, remove_unused_vertices);
}
template <typename Graph, typename PMap>
inline void remove_edge(
    edge<Graph> e, Graph& graph, PMap& pmap,
    bool remove_unused_vertices = false) {
  detail::graph_utils_impl<Graph, PMap> impl(graph, pmap);
  impl.remove_edge(e, remove_unused_vertices);
}
template <typename Graph, typename PMap, typename Map>
inline void remove_edge(
    edge<Graph> e, Graph& graph, PMap& pmap, Map& map,
    bool remove_unused_vertices = false) {
  detail::graph_utils_impl_with_map<Graph, Map, PMap> impl(graph, map, pmap);
  impl.remove_edge(e, remove_unused_vertices);
}
template <typename Graph, typename Map>
inline void remove_edge_with_map(edge<Graph> e, Graph& graph, Map& map) {
  detail::graph_utils_impl_with_map<Graph, Map> impl(graph, map);
  impl.remove_edge(e);
}
template <typename Graph, typename Map, typename PMap>
inline void remove_edge_with_map(
    edge<Graph> e, Graph& graph, Map& map, PMap& pmap) {
  detail::graph_utils_impl_with_map<Graph, Map, PMap> impl(graph, map, pmap);
  impl.remove_edge(e);
}

// simple wrapper around boost::copy_graph that serves
// to transform a directed into an undirected graph
template <typename GraphOut, typename GraphIn>
GraphOut symmetrise(const GraphIn& g) {
  static_assert(
      !is_directed<GraphOut>,
      "The GraphOut template type should be undirected. Use "
      "boost::copy_graph otherwise");
  GraphOut out;
  boost::copy_graph(g, out);
  return out;
}

template <typename GraphOut, typename GraphIn>
GraphOut symmetrise(
    const GraphIn& g, boost::bimap<vertex<GraphIn>, vertex<GraphOut>>& v_map) {
  static_assert(
      !is_directed<GraphOut>,
      "The GraphOut template type should be undirected. Use "
      "boost::copy_graph otherwise");
  GraphOut out;
  // when vertices are copied, do the default thing but keep track of the map
  boost::detail::vertex_copier<GraphIn, GraphOut> v_copier =
      boost::detail::make_vertex_copier(g, out);
  boost::copy_graph(
      g, out,
      boost::vertex_copy(
          [&v_map, &v_copier](vertex<GraphIn> old_v, vertex<GraphOut> new_v) {
            v_copier(old_v, new_v);
            v_map.insert({old_v, new_v});
          }));
  return out;
}

template <typename Graph>
std::size_t max_degree(const Graph& g) {
  auto degree = [&g](vertex<Graph> v) { return boost::degree(v, g); };
  auto [from, to] = boost::vertices(g);
  vertex<Graph> max_vertex =
      *std::max_element(from, to, detail::lt_with_key<vertex<Graph>>(degree));
  return degree(max_vertex);
}

template <typename Graph>
std::size_t min_degree(const Graph& g) {
  auto degree = [&g](vertex<Graph> v) { return boost::degree(v, g); };
  auto [from, to] = boost::vertices(g);
  vertex<Graph> min_vertex =
      *std::min_element(from, to, detail::lt_with_key<vertex<Graph>>(degree));
  return degree(min_vertex);
}

template <typename Graph>
std::set<vertex<Graph>> max_degree_nodes(const Graph& g) {
  const auto max = max_degree(g);
  auto [from, to] = boost::vertices(g);
  std::set<vertex<Graph>> out;
  std::copy_if(
      from, to, std::inserter(out, out.begin()),
      [&g, max](vertex<Graph> v) { return boost::degree(v, g) == max; });
  return out;
}

template <typename Graph>
std::set<vertex<Graph>> min_degree_nodes(const Graph& g) {
  const auto min = min_degree(g);
  auto [from, to] = boost::vertices(g);
  std::set<vertex<Graph>> out;
  std::copy_if(
      from, to, std::inserter(out, out.begin()),
      [&g, min](vertex<Graph> v) { return boost::degree(v, g) == min; });
  return out;
}

}  // namespace utils
}  // namespace tket::graphs
