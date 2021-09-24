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

#ifndef _TKET_UIDConnectivity_H
#define _TKET_UIDConnectivity_H

#include <algorithm>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <type_traits>

#include "Graphs/TreeSearch.hpp"
#include "Graphs/Utils.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/UnitID.hpp"

namespace tket::graph {

namespace detail {

class UIDDoesNotExistError : public std::logic_error {
  using std::logic_error::logic_error;
};
class EdgeDoesNotExistError : public std::logic_error {
  using std::logic_error::logic_error;
};

/** vertices of connectivity graph */
template <typename UID_t>
struct UIDVertex {
  UIDVertex() : uid() {}
  explicit UIDVertex(const UID_t _uid) : uid(_uid) {}

  bool operator==(const UIDVertex<UID_t> other) { return uid == other.uid; }
  bool operator<(const UIDVertex<UID_t>& other) const {
    return uid < other.uid;
  }

  UID_t uid;
};

/** edges of connectivity graph */
struct UIDInteraction {
  explicit UIDInteraction(unsigned d = 1) : weight(d) {}

  unsigned weight;
};

/**
 * Base class of UIDConnectivity, where all the implementation resides
 *
 * UIDConnectivity instances are directed graphs where vertices are given by
 * `UnitID`s, or one of the subtypes `Qubit` or `Node`. It is a wrapper around a
 * BGL graph that provides a clean class API, taking care of mapping all BGL
 * vertices and edge pointers to UnitIDs, respectively pairs of UnitIDs.
 *
 * This class is used mainly by the Architecture and QubitGraph classes; the
 * former inherits from UIDConnnectivity<Node>, the latter from
 * UIDConnnectivity<Qubit>.
 *
 * The vertices and edges can be given weights of type double if desired, and
 * the underlying undirected graph can be computed. The container types of the
 * underlying boost graph are specified as further template parameters, but have
 * sensible defaults in the UIDConnnectivity class.
 */
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
class UIDConnectivityBase {
 protected:
  using ConnGraph = boost::adjacency_list<
      OutEdgeListS, VertexListS, boost::bidirectionalS, UIDVertex<UID_t>,
      UIDInteraction>;
  using UndirectedConnGraph = boost::adjacency_list<
      boost::setS, VertexListS, boost::undirectedS, UIDVertex<UID_t>,
      UIDInteraction>;
  using Vertex = utils::vertex<ConnGraph>;
  using UndirectedVertex = utils::vertex<UndirectedConnGraph>;
  using Connection = std::pair<UID_t, UID_t>;

 public:
  /** empty default constructor */
  UIDConnectivityBase() : graph(), uid_to_vertex() {}

  /** constructor from list of vertices */
  explicit UIDConnectivityBase(const std::vector<UID_t>& uids);

  /** constructor from list of edges */
  explicit UIDConnectivityBase(const std::vector<Connection>& edges);

  /** add vertex to interaction graph */
  void add_uid(const UID_t uid) {
    Vertex v = boost::add_vertex(UIDVertex(uid), graph);
    uid_to_vertex.left.insert({uid, v});
  }

  /** remove vertex from the interaction graph */
  void remove_uid(const UID_t uid) {
    if (!uid_exists(uid)) {
      throw UIDDoesNotExistError(
          "The UID passed to UIDConnectivity::remove_uid must "
          "exist!");
    }
    Vertex v = to_vertices(uid);
    boost::clear_vertex(v, graph);
    utils::remove_vertex_with_map(v, graph, uid_to_vertex.right);
  }

  /** add edge to the interaction graph */
  void add_connection(const UID_t uid1, const UID_t uid2, unsigned weight = 1);

  /** remove edges in the connection graph */
  void remove_connections(
      const std::vector<Connection>& edges,
      bool remove_unused_vertices = false);
  void remove_connection(
      const Connection edge, bool remove_unused_vertices = false);
  void remove_connection(
      const UID_t uid1, const UID_t uid2, bool remove_unused_vertices = false);

  /** checks if edge exists */
  bool connection_exists(const UID_t uid1, const UID_t uid2) const;

  /** check if uid is a vertex */
  bool uid_exists(const UID_t uid) const;

  /** returns connection weight between two UnitID */
  unsigned get_connection_weight(const UID_t uid1, const UID_t uid2) const;

  /** return vertex degree of UnitID */
  unsigned get_degree(const UID_t uid) const;

  /** max depth from `root` in grahp */
  std::size_t get_max_depth(const UID_t root) const;

  /** return vertex out degree of UnitID */
  unsigned get_out_degree(const UID_t uid) const;

  /** number of vertices */
  unsigned n_uids() const { return boost::num_vertices(graph); }

  /** number of edges in graph */
  inline unsigned n_connections() const { return boost::num_edges(graph); }

  /** number of vertices with deg > 0 */
  inline unsigned n_connected() const {
    auto [beg, end] = boost::vertices(graph);
    auto nonzero_deg = [this](Vertex v) { return boost::degree(v, graph) > 0; };
    return std::count_if(beg, end, nonzero_deg);
  }

  /** get all connections in a vector */
  std::set<Connection> get_connections_set() const;
  std::vector<Connection> get_connections_vec() const;

  /** returns an unweighted undirected graph
   *  with the underlying connectivity */
  UndirectedConnGraph get_undirected_connectivity() const;

  /** Run bfs on underlying undirected subgraph */
  std::vector<std::size_t> get_distances(const UID_t root) const;
  std::size_t get_distance(const UID_t uid1, const UID_t uid2) const;

  /** remove vertices with deg == 0 */
  inline void remove_stray_uids() {
    utils::remove_stray_vertices_with_map(graph, uid_to_vertex.right);
  }

  /** set of all UnitIDs in interaction graph */
  std::set<UID_t> get_all_uids_set() const;
  std::vector<UID_t> get_all_uids_vec() const;
  /** iterator of all UnitIDs in interaction graph */
  auto get_all_uids() const {
    return get_vertices_it() |
           boost::adaptors::transformed(
               [this](const auto& v) { return get_uid(v); });
  }

  /* return UIDs with greatest (undirected) degree in graph */
  std::set<UID_t> max_degree_uids() const;

  /** return UIDs with smallest (undirected) degree */
  std::set<UID_t> min_degree_uids() const;

  /** returns path from `root` to `target` */
  std::vector<UID_t> get_path(const UID_t root, const UID_t target) const;

  /** get undirected adjacent UnitIDs */
  std::set<UID_t> get_neighbour_uids(const UID_t uid) const;

  bool operator==(
      const UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>& other) const;

 protected:
  /** get edge iterator */
  auto get_edges_it() const {
    return boost::make_iterator_range(boost::edges(graph));
  }

  /** get vertex iterator */
  auto get_vertices_it() const {
    return boost::make_iterator_range(boost::vertices(graph));
  }

  using vertex_bimap = boost::bimap<UID_t, Vertex>;
  using left_map = typename vertex_bimap::left_map;
  using right_map = typename vertex_bimap::right_map;
  /* get UnitID from vertex */
  const UnitID& get_uid(Vertex v) const { return graph[v].uid; }

  left_map& to_vertices() { return uid_to_vertex.left; }
  const left_map& to_vertices() const { return uid_to_vertex.left; }
  Vertex to_vertices(const UID_t& uid) const { return to_vertices().at(uid); }

  right_map& from_vertices() { return uid_to_vertex.right; }
  const right_map& from_vertices() const { return uid_to_vertex.right; }
  UID_t from_vertices(Vertex v) const { return from_vertices().at(v); }

  ConnGraph graph;
  vertex_bimap uid_to_vertex;
};

}  // namespace detail

/**
 * UIDConnectivity instances are graphs of UnitID vertices.
 * It should be instantiated with UnitIDs, or one of its subtypes `Qubit` or
 * `Node`
 *
 * All functionality for this class is implemented in the base class
 * UIDConnectivityBase. This class only adds caching of some function calls for
 * efficiency, innvalidating cache in case of changes on the underlying graph
 */
template <
    typename UID_t, typename OutEdgeListS = boost::vecS,
    typename VertexListS = boost::vecS>
class UIDConnectivity
    : public detail::UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS> {
 private:
  using Base = detail::UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>;

 public:
  using Base::Base;
  /* useful type synonyms */
  using UndirectedConnGraph = typename Base::UndirectedConnGraph;
  using ConnGraph = typename Base::ConnGraph;
  using Connection = typename Base::Connection;
  using Vertex = typename Base::Vertex;

  // We cache distances. A value of zero in the cache implies that the nodes are
  // disconnected (unless they are equal).
  const std::vector<std::size_t>& get_distances(const UnitID& root) const&;
  std::vector<std::size_t>&& get_distances(const UnitID& root) const&&;

  /**
   * Graph distance between two nodes.
   *
   * @param uid1 first node
   * @param uid2 second node
   *
   * @return length of shortest path between the nodes
   * @throws UIDsNotConnected if there is no path between the nodes
   */
  std::size_t get_distance(const UID_t uid1, const UID_t uid2) const;

  /** Returns all nodes at a given distance from a given 'source' node */
  std::vector<UID_t> uids_at_distance(
      const UID_t& root, std::size_t distance) const;

  // we cache the undirected graph
  const UndirectedConnGraph& get_undirected_connectivity() const&;
  UndirectedConnGraph&& get_undirected_connectivity() const&&;

  // these functions invalidate caching: invalidate cache then call Base
  // function
  void remove_uid(const UID_t uid);
  void add_uid(const UID_t uid);

  void remove_stray_uids();
  void add_connection(const UID_t uid1, const UID_t uid2, unsigned weight = 1);
  void remove_connections(const std::vector<Connection>& edges);
  void remove_connection(
      const Connection edge, bool remove_unused_vertices = false);
  void remove_connection(
      const UID_t uid1, const UID_t uid2, bool remove_unused_vertices = false);

 private:
  inline void invalidate_cache() {
    distance_cache.clear();
    undir_graph = std::nullopt;
  }
  mutable std::map<UnitID, std::vector<std::size_t>> distance_cache;
  mutable std::optional<UndirectedConnGraph> undir_graph;
};

/**
 * Exception thrown because two nodes are disconnected from one another.
 *
 * @tparam UID_t node type
 */
template <typename UID_t>
class UIDsNotConnected : public std::logic_error {
 public:
  UIDsNotConnected(const UID_t& uid1, const UID_t uid2)
      : std::logic_error(
            uid1.repr() + " and " + uid2.repr() + " are not connected") {}
};

// template explicit instations, with implementations in cpp file
extern template class detail::UIDConnectivityBase<
    UnitID, boost::vecS, boost::vecS>;
extern template class detail::UIDConnectivityBase<
    Node, boost::vecS, boost::vecS>;
extern template class detail::UIDConnectivityBase<
    Qubit, boost::vecS, boost::vecS>;
extern template struct detail::UIDVertex<UnitID>;
extern template struct detail::UIDVertex<Node>;
extern template struct detail::UIDVertex<Qubit>;
extern template class UIDConnectivity<UnitID>;
extern template class UIDConnectivity<Node>;
extern template class UIDConnectivity<Qubit>;

}  // namespace tket::graph

#endif
