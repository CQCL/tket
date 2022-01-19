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

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include "Graphs/Utils.hpp"

// implementation details of `TreeSearch.hpp`

namespace tket::graphs {

// ::detail includes all internal implementation-specific details
namespace detail {

/* Helpers `run_bfs` and `run_dfs` return instances of TreeSearchBase
 * subclasses. These provide a host of useful member functions to obtain the
 * relevant results from the tree search.
 *
 * Note that run() must have been called before
 * using the member functions to retrieve results.
 * This is done automatically by the utility functions in `TreeSearch.hpp`
 *
 * First template argument is Graph type, second is the boost property map type
 * that should be used as a vertex index.
 */
template <typename Graph, typename BoostPMap>
class TreeSearchBase {
 protected:
  using vertex = utils::vertex<Graph>;
  using size_t = std::size_t;
  using parent_vec = std::vector<vertex>;
  using color_vec = std::vector<boost::default_color_type>;
  using index_pmap_t = BoostPMap;
  using dist_pmap_t =
      boost::iterator_property_map<typename dist_vec::iterator, index_pmap_t>;
  using parent_pmap_t =
      boost::iterator_property_map<typename parent_vec::iterator, index_pmap_t>;
  using color_pmap_t =
      boost::iterator_property_map<typename color_vec::iterator, index_pmap_t>;
  using visitor_t = std::pair<
      boost::distance_recorder<dist_pmap_t, boost::on_tree_edge>,
      boost::predecessor_recorder<parent_pmap_t, boost::on_tree_edge>>;

 public:
  /* Constructor: search through Graph g with root v, using index map pmap */
  TreeSearchBase(vertex v, Graph&& g, BoostPMap pmap)
      : root(v),
        to_index(pmap),
        graph(g),
        dists(boost::num_vertices(g)),
        parents(boost::num_vertices(g)),
        color_map(boost::num_vertices(g)),
        visitor{
            boost::record_distances(
                boost::make_iterator_property_map(dists.begin(), to_index),
                boost::on_tree_edge()),
            boost::record_predecessors(
                boost::make_iterator_property_map(parents.begin(), to_index),
                boost::on_tree_edge())} {
    // for every vertex, assign self as parent
    for (auto v : boost::make_iterator_range(boost::vertices(graph))) {
      parents[to_index[v]] = v;
    }
  }

  /* Get results from tree search */
  /* vector of vertex predecessor in search */
  const parent_vec& get_parents() const& { return parents; }
  const parent_vec&& get_parents() const&& { return std::move(parents); }
  /* vector of distances from the root */
  const dist_vec& get_dists() const& { return dists; }
  const dist_vec&& get_dists() const&& { return std::move(dists); }
  size_t get_dist(const vertex& v) const { return dists[to_index[v]]; }

  /* Rebase tree search from some new root */
  void change_root(vertex v) {
    // naive implementation
    // This might be made more efficient in the subclasses
    vertex old_root = root;
    root = v;
    if (root != old_root) {
      std::fill(dists.begin(), dists.end(), 0);
      std::iota(parents.begin(), parents.end(), 0);
      run();
    }
  }

  /* Vector from target to root (following reverse edges) */
  parent_vec path_to_root(vertex target) const {
    std::vector<vertex> path;

    vertex current = target;
    path.push_back(target);
    while (current != root) {
      vertex parent = parents[to_index[current]];
      if (parent == current) {
        // graph is disconnected and there is no path to root
        return {};
      }
      current = parent;
      path.push_back(current);
    }
    return path;
  }

  /* Vector from root to target along graph edges */
  parent_vec path_from_root(vertex target) const {
    parent_vec out = path_to_root(target);
    std::reverse(out.begin(), out.end());
    return out;
  }

  /* Depth of tree search */
  size_t max_depth() const {
    auto it = std::max_element(dists.cbegin(), dists.cend());
    if (it == dists.cend()) {
      throw std::invalid_argument(
          "TreeSearch::max_depth: There is no entry in distance "
          "vector");
    }
    return *it;
  }

  /* Vertex of max depth
   * if more than one exist, returns first found */
  vertex max_depth_vertex() const {
    const size_t target_depth = max_depth();
    for (vertex v : boost::make_iterator_range(boost::vertices(this->graph))) {
      if (get_dist(v) == target_depth) {
        return v;
      }
    }
    throw std::runtime_error("max_depth_vertex: No vertex has maximal depth");
  }

  /* Path from root to deepest vertex */
  parent_vec longest_path() const { return path_from_root(max_depth_vertex()); }

  /* Run tree search */
  virtual void run() = 0;

 protected:
  template <typename Visitor>
  auto get_default_named_params(Visitor& vis) {
    return boost::visitor(vis).color_map(
        boost::make_iterator_property_map(color_map.begin(), to_index));
  }

  vertex root;            // root of tree search
  index_pmap_t to_index;  // vertex index pmap
  Graph graph;            // graph to be searched
  dist_vec dists;         // distance from root
  parent_vec parents;     // parent towards root
  color_vec color_map;    // color used by tree search algorithm
  visitor_t visitor;      // visitor passed to boost tree search alg
};

/* A vertex index map, mapping boost vertex descriptors to unsigned
 * is required for tree searches (to index the parents and dists vectors).
 *
 * The TreeSearch class accepts a std::map-like PMap argument to define such a
 * mapping. This is converted to a boost property map using
 * boost::associative_property_map. If none is provided, the default
 * boost::vertex_index is used (by default, this is available on boost graphs
 * using the vecS underlying EdgeOutList type
 */

/* Note on the implementation:
 * The third template argument (a bool) should always be left to its
 * default value. It is used to flag whether the Graph type considered
 * has an implicit integer vertex index.
 * This flag is used in the template specialisation below to change
 * the template definition when `HasImplicitIndex == True`,
 * i.e. the first definition of `TreeSearch` is for the case `HasImplicitIndex
 * == False` while the template partial specialisation below is for the case
 * `HasImplicitIndex == True`.
 */
template <
    typename Graph, typename PMap = void,
    bool HasImplicitIndex = utils::detail::has_integral_index<Graph>>
class TreeSearch
    : public TreeSearchBase<Graph, boost::associative_property_map<PMap>> {
 private:
  using Base = TreeSearchBase<Graph, boost::associative_property_map<PMap>>;
  using vertex = typename Base::vertex;

 public:
  TreeSearch(vertex v, Graph&& g, PMap& pmap)
      : Base(v, std::forward<Graph>(g), boost::make_assoc_property_map(pmap)) {}
};

template <typename Graph>
class TreeSearch<Graph, void, true>
    : public TreeSearchBase<
          Graph,
          decltype(boost::get(boost::vertex_index, std::declval<Graph>()))> {
 private:
  using PMap = decltype(boost::get(boost::vertex_index, std::declval<Graph>()));
  using Base = TreeSearchBase<Graph, PMap>;
  using vertex = typename Base::vertex;

 public:
  TreeSearch(vertex v, Graph&& g)
      : Base(v, std::forward<Graph>(g), boost::get(boost::vertex_index, g)) {}
};

}  // namespace detail

/* Return type of `run_bfs`.
 *
 * For each vertex, computes the BFS distance to the root
 * as well as the next node on the path towards root (parent).
 * These are stored in the base class TreeSearchBase and can
 * be retrieved through one of the member functions defined therein
 */
template <typename... Args>
class BFS : public detail::TreeSearch<Args...> {
 private:
  using Base = detail::TreeSearch<Args...>;

 public:
  // inherit constructors
  using Base::Base;

  using parent_vec = typename Base::parent_vec;

  /* Run bfs */
  void run() override {
    auto vis = boost::make_bfs_visitor(this->visitor);
    boost::breadth_first_search(
        this->graph, this->root, this->get_default_named_params(vis));
  }
};

/* Return type of `run_dfs`.
 *
 * For each vertex, computes the DFS distance to the root
 * as well as the next node on the path towards root (parent).
 * These are stored in the base class TreeSearchBase and can
 * be retrieved through one of the member functions defined therein
 */
template <typename... Args>
class DFS : public detail::TreeSearch<Args...> {
 private:
  using Base = detail::TreeSearch<Args...>;

 public:
  // inherit constructors
  using Base::Base;

  using parent_vec = typename Base::parent_vec;

  /* Run dfs */
  void run() override {
    auto vis = boost::make_dfs_visitor(this->visitor);
    boost::depth_first_search(
        this->graph,
        this->get_default_named_params(vis).root_vertex(this->root));
  }
};

}  // namespace tket::graphs
