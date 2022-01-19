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
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <map>
#include <set>
#include <type_traits>
#include <utility>

#include "Utils/BiMapHeaders.hpp"
#include "Utils/GraphHeaders.hpp"

// Implemetation details of `Graphs/Utils.hpp`
namespace tket::graphs::utils {

// ::detail includes all internal implementation-specific details
namespace detail {

// Type synonyms for many often-used BGL types
template <typename T>
using vertex =
    typename boost::graph_traits<std::remove_reference_t<T>>::vertex_descriptor;
template <typename T>
using edge =
    typename boost::graph_traits<std::remove_reference_t<T>>::edge_descriptor;
template <typename T>
using vertex_index_pmap = typename boost::property_map<
    std::remove_reference_t<T>, boost::vertex_index_t>::type;
template <typename T>
using vertex_index =
    typename boost::property_traits<vertex_index_pmap<T>>::value_type;

/* Compile-time boolean flag to check if boost graph type is
 * directed or undirected graph
 *
 * Template Argument Graph is the boost type of the graph to be checked
 * the second (directedness tag) should be left to default to be inferred
 * by the compiler ("SFINAE")
 *
 * This is available in the public interface as is_directed<Graph> (see
 * `Utils.hpp`)
 */
template <
    typename Graph, typename directed_category = typename boost::graph_traits<
                        std::remove_reference_t<Graph>>::directed_category>
struct is_directed_struct {
  static const bool value = true;
};
template <typename Graph>
struct is_directed_struct<Graph, boost::undirected_tag> {
  static const bool value = false;
};

/* Compile-time boolean used in template specialisations
 * and graph_utils_helper_with_map to infer whether the Graph type supports
 * integer indexing
 *
 * has_integral_index<Graph, PMap> is true iff
 * either the vertex_descriptor of Graph are integers (i.e. VertexList = vecS)
 * or an explicit property map PMap is provided with integral values
 *
 * Template argument Graph is the boost type of the graph to be checked
 * PMap is the optional property map to provide integral index values
 * The third argument should always be left to default, used for type inference
 * ("SFINAE")
 */
template <typename Graph, typename PMap = void, typename = void>
struct has_integral_index_struct {
  static const bool value = false;
};
template <typename Graph>
struct has_integral_index_struct<
    Graph, void, std::enable_if_t<std::is_integral_v<vertex<Graph>>>> {
  static const bool value = true;
};
template <typename Graph, typename PMap>
struct has_integral_index_struct<
    Graph, PMap,
    std::enable_if_t<std::is_integral_v<typename PMap::mapped_type>>> {
  static const bool value = true;
};
template <typename Graph, typename PMap = void>
constexpr bool has_integral_index =
    has_integral_index_struct<Graph, PMap>::value;

// less than comparison using key as compared values
template <typename V, typename Func>
auto lt_with_key(const Func& key) {
  return [&key](const V& v1, const V& v2) { return key(v1) < key(v2); };
}

/* This class, along with their inherited classes and template specialisations
 * below implement the graphs::utils helper functions, defined further below.
 *
 *                                 CLASS HIERARCHY
 *                                 ---------------
 *                       +---------------------+   +---------------------------+
 * Base classes:         |   graph_utils_base  |   | graph_utils_base_with_ind |
 *                       +---------------------+   +---------------------------+
 *                             |         |                  |           |
 *                             v         v                  v           v
 *  Container-specific   +---------+ +---------+       +---------+ +---------+
 *  classes              | _helper | | _helper |       | _helper | | _helper |
 *  (pmap / indexable)   | YES/NO  | |  NO/NO  |       | YES/YES | |  NO/YES |
 *                       +---------+ +---------+       +---------+ +---------+
 *                           |  |       |  |              |  |        |  |
 *                           |  |       |  |  +-------+   |  |        |  |
 *  Child class to be used   |  |       |  +->| _impl |<--+  |        |  |
 *  in utils                 |  +------------>|       |<--------------+  |
 *                           |          |     +-------+      |           |
 *                           v          v                    v           v
 *  Container-specific   +---------+ +---------+       +---------+ +---------+
 *  classes with added   | _helper | | _helper |       | _helper | | _helper |
 *  support for maps     |_with_map| |_with_map|       |_with_map| |_with_map|
 *  (pmap / indexable)   | YES/NO  | |  NO/NO  |       | YES/YES | |  NO/YES |
 *                       +---------+ +---------+       +---------+ +---------+
 *                            |          |                  |          |
 *                            |          |  +---------+     |          |
 *  Child class to be used    |          +->|  _impl  |<----+          |
 *  in utils with map support +------------>|_with_map|<---------------+
 *                                          +---------+
 *
 * graph_utils_base and graph_utils_base_with_ind implement the functionality
 * that is shared among all vertex container types (i.e. most of the code logic)
 *
 * Impl must be a subclass of this and must implement
 *  - on_remove_vertex: event handler when vertex is removed
 */
template <typename Graph, typename Impl>
class graph_utils_base {
 public:
  explicit graph_utils_base(Graph& g) : graph(g) {}

  void remove_vertex(vertex<Graph> v) {
    remove_vertex_handler(v);
    boost::remove_vertex(v, graph);
  }

  void remove_edge(edge<Graph> e, bool remove_unused_vertices = false) {
    auto [u, v] = get_vertex_pair(e);
    boost::remove_edge(e, graph);

    if (remove_unused_vertices) {
      if (boost::degree(u, graph) == 0) {
        remove_vertex(u);
      }
      if (boost::degree(v, graph) == 0) {
        remove_vertex(v);
      }
    }
  }

 protected:
  virtual void remove_vertex_handler(vertex<Graph> v) {
    Impl* impl = static_cast<Impl*>(this);
    impl->on_remove_vertex(v);
  }

  /* returns source and target in an order in which they can be
   * safely deleted if necessary
   */
  std::pair<vertex<Graph>, vertex<Graph>> get_vertex_pair(edge<Graph> e) {
    vertex<Graph> u = boost::source(e, graph);
    vertex<Graph> v = boost::target(e, graph);

    return in_delete_order(u, v);
  }

  virtual std::pair<vertex<Graph>, vertex<Graph>> in_delete_order(
      vertex<Graph> u, vertex<Graph> v) {
    return {u, v};
  }

  Graph& graph;
};

/* small adjustments to graph_utils_helper for integral vertex indices.
 * Unlike the general case, here indices get shifted when vertices are deleted
 * so that the vertex indices remain contiguous
 *
 * Impl must be a subclass of this and must implement
 *  - to_index: return index of vertex v
 *  - typename index_type: type of index returned by to_index
 *  - on_index_change: event handler when indices are shifted
 *  - on_remove_vertex: event handler when vertex is removed, after all
 *    indices were updated
 */
template <typename Graph, typename Impl>
class graph_utils_base_with_ind : public graph_utils_base<Graph, Impl> {
 private:
  using Base = graph_utils_base<Graph, Impl>;

 public:
  explicit graph_utils_base_with_ind(Graph& g) : Base(g) {}

 protected:
  using Base::graph;
  void remove_vertex_handler(vertex<Graph> v) override {
    Impl* impl = static_cast<Impl*>(this);
    using index_type = typename Impl::index_type;

    for (vertex<Graph> u : boost::make_iterator_range(boost::vertices(graph))) {
      if (impl->to_index(u) > impl->to_index(v)) {
        index_type new_ind = impl->to_index(u) - 1;
        impl->on_index_change(u, new_ind);
      }
    }
    impl->on_remove_vertex(v);
  }

  std::pair<vertex<Graph>, vertex<Graph>> in_delete_order(
      vertex<Graph> u, vertex<Graph> v) override {
    Impl* impl = static_cast<Impl*>(this);

    // return larger index first
    if (impl->to_index(u) < impl->to_index(v)) {
      return {v, u};
    } else {
      return {u, v};
    }
  }
};

/* The classes graph_utils_helper and graph_utils_helper_with_map provide
 * the container-specific implementation bits.
 * The Derived template type is the child class that inherits from this class
 *  -> see Curiously Recurring Template Pattern (CRTP)
 *  The Indexable boolean template is true iff Graph support integer vertex
 *  indexing (see has_integral_index). It should be left to default and infered
 *  by compiler
 *
 * We distinguish four cases for template specialisations, given by the four
 * combinations of two variables:
 *
 *      a) whether an explicit property map is provided or not
 *      b) whether the graph is indexable, i.e. either it has an intrinsic
 * vertex_index property (case VertexList == vecS) or the property map has
 * integral values
 *
 * Case 1: explicit property map, not indexable
 */
template <
    typename Derived, typename Graph, typename PMap = void,
    bool Indexable = has_integral_index<Graph, PMap>>
class graph_utils_helper : public graph_utils_base<Graph, Derived> {
 private:
  using Base = graph_utils_base<Graph, Derived>;

 public:
  graph_utils_helper(Graph& g, PMap& p) : Base(g), pmap(p) {}

  void on_remove_vertex(vertex<Graph> v) { pmap.erase(v); }

 private:
  PMap& pmap;
};

// Case 2: no explicit property, not indexable
template <typename Derived, typename Graph>
class graph_utils_helper<Derived, Graph, void, false>
    : public graph_utils_base<Graph, Derived> {
 private:
  using Base = graph_utils_base<Graph, Derived>;

 public:
  explicit graph_utils_helper(Graph& g) : Base(g) {}
  void on_remove_vertex(vertex<Graph>) {}
};

// Case 3: explicit property, indexable
// in this case we must ensure that we update the property map
// when vertices are removed and indices get shifted
template <typename Derived, typename Graph, typename PMap>
class graph_utils_helper<Derived, Graph, PMap, true>
    : public graph_utils_base_with_ind<Graph, Derived> {
 private:
  using Base = graph_utils_base_with_ind<Graph, Derived>;

 public:
  using index_type = typename PMap::mapped_type;
  graph_utils_helper(Graph& g, PMap& p) : Base(g), pmap(p), new_pmap(p) {
    static_assert(
        !std::is_integral_v<vertex<Graph>>,
        "Ambiguous vertex index: you cannot provide an explicit "
        "index property map if vertex descriptors are integrals.");
  }

  index_type to_index(vertex<Graph> v) { return pmap[v]; }

  void on_index_change(vertex<Graph> v, index_type new_ind) {
    new_pmap[v] = new_ind;
  }
  void on_remove_vertex(vertex<Graph> v) {
    new_pmap.erase(v);
    pmap = new_pmap;
  }

 private:
  PMap& pmap;
  PMap new_pmap;
};

// Case 4: no explicit property, indexable
// boost automatically shifts indices, so there is not much to do
template <typename Derived, typename Graph>
class graph_utils_helper<Derived, Graph, void, true>
    : public graph_utils_base_with_ind<Graph, Derived> {
 private:
  using Base = graph_utils_base_with_ind<Graph, Derived>;

 public:
  using index_type = vertex_index<Graph>;
  explicit graph_utils_helper(Graph& g) : Base(g) {}

  index_type to_index(vertex<Graph> v) { return v; }
  void on_index_change(vertex<Graph>, index_type) {}
  void on_remove_vertex(vertex<Graph>) {}
};

/* Now come the same 4 cases with the additional support of an external map
 * with vertrex indices as keys
 */
// Case 1-2: not indexable
template <
    typename Derived, typename Graph, typename Map, typename PMap = void,
    bool Indexable = has_integral_index<Graph, PMap>>
class graph_utils_helper_with_map
    : public graph_utils_helper<Derived, Graph, PMap> {
 private:
  using Base = graph_utils_helper<Derived, Graph, PMap>;

 public:
  graph_utils_helper_with_map(Graph& g, Map& m) : Base(g), map(m) {}

  void on_remove_vertex(vertex<Graph> v) {
    map.erase(v);
    Base::on_remove_vertex(v);
  }

 private:
  Map& map;
};

// Case 3: explicit property, indexable
template <typename Derived, typename Graph, typename Map, typename PMap>
class graph_utils_helper_with_map<Derived, Graph, Map, PMap, true>
    : public graph_utils_helper<Derived, Graph, PMap> {
 private:
  using Base = graph_utils_helper<Derived, Graph, PMap>;

 public:
  using index_type = typename Base::index_type;
  graph_utils_helper_with_map(Graph& g, Map& m, PMap& p)
      : Base(g, p), map(m), new_map(m) {}

  void on_index_change(vertex<Graph> v, index_type new_ind) {
    if (new_map.count(new_ind)) {
      new_map.erase(new_ind);
    }
    new_map.insert({new_ind, new_map[Base::to_index(v)]});
    Base::on_index_change(v, new_ind);
  }
  void on_remove_vertex(vertex<Graph> v) {
    auto max_ind = std::max_element(std::begin(new_map), std::end(new_map));
    new_map.erase(max_ind);
    map.clear();
    for (auto& p : new_map) {
      map.insert(p);
    }

    Base::on_remove_vertex(v);
  }

 private:
  Map& map;
  Map new_map;
};

// Case 4: no explicit property, indexable
template <typename Derived, typename Graph, typename Map>
class graph_utils_helper_with_map<Derived, Graph, Map, void, true>
    : public graph_utils_helper<Derived, Graph, void> {
 private:
  using Base = graph_utils_helper<Derived, Graph, void>;
  using StdMap = std::map<typename Map::key_type, typename Map::mapped_type>;

  static typename StdMap::value_type transformer_it(
      typename Map::value_type& elem) {
    return std::make_pair(elem.first, elem.second);
  }
  auto map_begin_it(Map& m) {
    return boost::make_transform_iterator(m.begin(), &transformer_it);
  }
  auto map_end_it(Map& m) {
    return boost::make_transform_iterator(m.end(), &transformer_it);
  }

 public:
  using index_type = typename Base::index_type;
  graph_utils_helper_with_map(Graph& g, Map& m)
      : Base(g), map(m), new_map(map_begin_it(m), map_end_it(m)) {}
  void on_index_change(vertex<Graph> v, index_type new_ind) {
    if (new_map.count(new_ind)) {
      new_map.erase(new_ind);
    }
    new_map.insert({new_ind, new_map[Base::to_index(v)]});
    Base::on_index_change(v, new_ind);
  }
  void on_remove_vertex(vertex<Graph> v) {
    auto max_ind = std::max_element(std::begin(new_map), std::end(new_map));
    new_map.erase(max_ind);
    map.clear();
    for (auto& p : new_map) {
      map.insert(p);
    }

    Base::on_remove_vertex(v);
  }

 protected:
  Map& map;
  StdMap new_map;
};

// The above implementions are accessed through the following helper classes,
// making template args nicer
template <typename Graph, typename PMap = void>
class graph_utils_impl
    : public graph_utils_helper<graph_utils_impl<Graph, PMap>, Graph, PMap> {
 private:
  using Base = graph_utils_helper<graph_utils_impl<Graph, PMap>, Graph, PMap>;

 public:
  using Base::Base;
};
template <typename Graph, typename Map, typename PMap = void>
class graph_utils_impl_with_map
    : public graph_utils_helper_with_map<
          graph_utils_impl_with_map<Graph, Map, PMap>, Graph, Map, PMap> {
 private:
  using Base = graph_utils_helper_with_map<
      graph_utils_impl_with_map<Graph, Map, PMap>, Graph, Map, PMap>;

 public:
  using Base::Base;
};

}  // namespace detail
}  // namespace tket::graphs::utils
