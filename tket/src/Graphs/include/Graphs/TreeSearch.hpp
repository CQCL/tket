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

#include <utility>
#include <vector>

#include "Graphs/TreeSearch_impl.hpp"
#include "Graphs/Utils.hpp"

/* A wrapper around the BFS and DFS implementation of boost graph. The hope is
 * that this is less awkward (albeit less flexible) to use than the direct boost
 * code.
 *
 * The helpers `run_bfs` and `run_dfs`, as well as `longest_simple_path` are
 * provided (defined at the bottom) to be used in most situations.
 *
 * This file contains the public interface to be used,
 * and the file `TreeSearch_impl.hpp` contains the implementation details.
 * Note that these files are header-only as they rely all on templates.
 *
 * Examples (pseudo-code):
 *  >   // get distances of each vertex from root as vector
 *  >   auto dists = graphs::run_bfs(root, graph).get_dists();
 *  >   // get parent of each vertex of depth-first-search (DFS)
 *  >   auto parents = graphs::run_dfs(root, graph).get_parents()
 *  >
 *  >   // BFS and DFS require an integer vertex index.
 *  >   // This exists by default when the underlying edge type is vecS.
 *  >   // Otherwise, we can also pass a map-like Property Map that is used as
 * vertex index >   auto dfs = graphs::run_dfs(root, graph, vertex_index_map)
 *
 * `run_bfs` and `run_dfs` return instances of the BFS and DFS classes
 * defined in the `TreeSearch_impl.hpp` file. These store the result of the tree
 * search and allows the user to retrieve the relevant information as in the
 * examples above.
 */

namespace tket::graphs {

/* Tree search: Breadth First Search (BFS) and Depth First Search (DFS)
 *
 * Use the following helper functions to run BFS or DFS on boost graphs.
 * The returned instances are used to retrieve the computed results in the
 * desired form. See class TreeSearchBase for the list of member functions
 * available.
 */
template <typename Graph, typename PMap>
BFS<Graph, PMap> run_bfs(utils::vertex<Graph> root, Graph&& g, PMap& pmap) {
  BFS<Graph, PMap> impl(root, std::forward<Graph>(g), pmap);
  impl.run();
  return impl;
}
template <typename Graph>
BFS<Graph> run_bfs(utils::vertex<Graph> root, Graph&& g) {
  BFS<Graph> impl(root, std::forward<Graph>(g));
  impl.run();
  return impl;
}
template <typename Graph, typename PMap>
DFS<Graph, PMap> run_dfs(utils::vertex<Graph> root, Graph&& g, PMap& pmap) {
  DFS<Graph, PMap> impl(root, std::forward<Graph>(g), pmap);
  impl.run();
  return impl;
}
template <typename Graph>
DFS<Graph> run_dfs(utils::vertex<Graph> root, Graph&& g) {
  DFS<Graph> impl(root, std::forward<Graph>(g));
  impl.run();
  return impl;
}
template <typename Graph, typename Visitor>
DFS<Graph> run_dfs_with_visitor(
    utils::vertex<Graph> root, Graph&& g, Visitor vis) {
  DFS<Graph> impl(root, std::forward<Graph>(g));
  impl.run(vis);
  return impl;
}

/* Computes the longest simple path in a graph using DFS
 *
 * This is a naive implementation starting a new search from every
 * possible vertex
 */
template <typename Graph>
std::vector<utils::vertex<Graph>> longest_simple_path(
    const Graph& g, std::size_t cutoff_length = 0) {
  std::vector<utils::vertex<const Graph>> longest;
  auto root = boost::vertex(0, g);
  auto dfs = run_dfs(root, g);
  for (auto [it, end] = boost::vertices(g); it != end; ++it) {
    root = *it;
    dfs.change_root(root);
    if (dfs.max_depth() + 1 > longest.size()) {
      longest = dfs.longest_path();
      if (cutoff_length > 0 && longest.size() >= cutoff_length) {
        break;
      }
    }
  }
  return longest;
}

}  // namespace tket::graphs
