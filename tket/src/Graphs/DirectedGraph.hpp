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

#include "Graphs/AbstractGraph.hpp"
#include "Graphs/TreeSearch.hpp"
#include "Graphs/Utils.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket::graphs {

template <typename UID_t>
class UIDConnectivityBase : public AbstractGraph<UID_t> {
 protected:
  using ConnGraph = boost::adjacency_list<
      boost::vecS, boost::vecS, boost::bidirectionalS, UID_t, unsigned>;
  using UndirectedConnGraph = boost::adjacency_list<
      boost::setS, boost::vecS, boost::undirectedS, UID_t, unsigned>;
  using Vertex = utils::vertex<ConnGraph>;
  using UndirectedVertex = utils::vertex<UndirectedConnGraph>;
  using Connection = typename AbstractGraph<UID_t>::Edge;
  using AbstractGraph<UID_t>::nodes_;
  using AbstractGraph<UID_t>::node_exists;
  using AbstractGraph<UID_t>::edge_exists;

 public:
  /** empty default constructor */
  UIDConnectivityBase() : AbstractGraph<UID_t>(), graph(), uid_to_vertex() {}

  /** constructor from list of vertices */
  explicit UIDConnectivityBase(const std::vector<UID_t>& uids) : AbstractGraph<UID_t>(uids), graph() {
    for (const UID_t& uid : uids) {
      add_uid(uid);
    }
  }

  /** constructor from list of= edges */
  explicit UIDConnectivityBase(const std::vector<Connection>& edges) : graph() {
    for (auto [uid1, uid2] : edges) {
      if (!node_exists(uid1)) {
        add_uid(uid1);
      }
      if (!node_exists(uid2)) {
        add_uid(uid2);
      }
      boost::add_edge(to_vertices(uid1), to_vertices(uid2), graph);
    }
  }

  bool edge_exists(const UID_t& uid1, const UID_t& uid2) const override {
    if (!node_exists(uid1) || !node_exists(uid2)) {
      throw NodeDoesNotExistError(
          "The UIDs passed to UIDConnectivity::edge_exists must "
          "exist");
    }
    auto [_, exists] = boost::edge(to_vertices(uid1), to_vertices(uid2), graph);
    return exists;
  }

  /** add vertex to interaction graph */
  void add_uid(const UID_t uid) {
    nodes_.insert(uid);
    Vertex v = boost::add_vertex(uid, graph);
    uid_to_vertex.left.insert({uid, v});
  }

  /** remove vertex from the interaction graph */
  void remove_uid(const UID_t uid) {
    if (!node_exists(uid)) {
      throw NodeDoesNotExistError(
          "The UID passed to UIDConnectivity::remove_uid must "
          "exist!");
    }
    nodes_.erase(uid);
    Vertex v = to_vertices(uid);
    boost::clear_vertex(v, graph);
    utils::remove_vertex_with_map(v, graph, uid_to_vertex.right);
  }

  /** add weighted edge to the interaction graph */
  void add_connection(const UID_t uid1, const UID_t uid2, unsigned weight = 1) {
    if (!node_exists(uid1) || !node_exists(uid2)) {
      throw NodeDoesNotExistError(
          "The UIDs passed to UIDConnectivity::add_connection must "
          "exist");
    }
    boost::add_edge(to_vertices(uid1), to_vertices(uid2), weight, graph);
  }

  void remove_connection(
      const Connection edge, bool remove_unused_vertices = false) {
    if (!node_exists(edge.first) || !node_exists(edge.second)) {
      throw NodeDoesNotExistError(
          "Trying to remove an edge with non-existent vertices");
    }
    auto [e, exists] =
        boost::edge(to_vertices(edge.first), to_vertices(edge.second), graph);
    if (!exists) {
      throw EdgeDoesNotExistError(
          "The edge (" + edge.first.repr() + ", " + edge.second.repr() +
          ")"
          "cannot be removed as it does not exist");
    }
    utils::remove_edge_with_map(
        e, graph, uid_to_vertex.right, remove_unused_vertices);
  }

  void remove_connection(
      const UID_t uid1, const UID_t uid2, bool remove_unused_vertices = false) {
    remove_connection({uid1, uid2}, remove_unused_vertices);
  }

  /** remove edges in the connection graph */
  void remove_connections(
      const std::vector<Connection>& edges,
      bool remove_unused_vertices = false) {
    for (const auto& e : edges) {
      remove_connection(e, remove_unused_vertices);
    }
  }

  /** returns connection weight between two UnitID */
  unsigned get_connection_weight(const UID_t uid1, const UID_t uid2) const {
    if (!node_exists(uid1) || !node_exists(uid2)) {
      throw NodeDoesNotExistError(
          "Trying to retrieve edge weight from non-existent vertices");
    }
    auto [e, exists] = boost::edge(to_vertices(uid1), to_vertices(uid2), graph);
    if (!exists) {
      return 0.;
    }

    return graph[e];
  }

  /** return vertex degree of UnitID */
  unsigned get_degree(const UID_t uid) const {
    if (!node_exists(uid)) {
      throw NodeDoesNotExistError(
          "Trying to retrieve vertex degree from non-existent vertex");
    }
    return boost::degree(to_vertices(uid), graph);
  }

  /** max depth from `root` in grahp */
  std::size_t get_max_depth(const UID_t root) const {
    if (!node_exists(root)) {
      throw NodeDoesNotExistError("Trying to get depth from non-existent vertex");
    }
    return run_bfs(to_vertices(root), get_undirected_connectivity()).max_depth();
  }

  /** return vertex out degree of UnitID */
  unsigned get_out_degree(const UID_t uid) const {
    if (!node_exists(uid)) {
      throw NodeDoesNotExistError(
          "Trying to get outdegree from non-existent vertex");
    }
    return boost::out_degree(to_vertices(uid), graph);
  }

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
  std::set<Connection> get_connections_set() const {
    std::vector<Connection> vec = get_connections_vec();
    return std::set(vec.begin(), vec.end());
  }

  std::vector<Connection> get_connections_vec() const {
    std::vector<Connection> out;
    for (auto e : get_edges_it()) {
      UID_t source = graph[boost::source(e, graph)];
      UID_t target = graph[boost::target(e, graph)];
      out.push_back({source, target});
    }
    return out;
  }

  /** returns an unweighted undirected graph
   *  with the underlying connectivity */
  UndirectedConnGraph get_undirected_connectivity() const {
    return graphs::utils::symmetrise<UndirectedConnGraph>(graph);
  }

  /** Run bfs on underlying undirected subgraph */
  std::vector<std::size_t> get_distances(const UID_t root) const {
    if (!node_exists(root)) {
      throw NodeDoesNotExistError(
          "Trying to get distances from non-existent root vertex");
    }
    return run_bfs(to_vertices(root), get_undirected_connectivity()).get_dists();
  }

  std::size_t get_distance(const UID_t uid1, const UID_t uid2) const {
    if (uid1 == uid2) {
      return 0;
    }
    size_t d = get_distances(uid1)[to_vertices(uid2)];
    if (d == 0) {
      throw UIDsNotConnected(uid1, uid2);
    }
    return d;
  }

  /** remove vertices with deg == 0 */
  inline void remove_stray_uids() {
    utils::remove_stray_vertices_with_map(graph, uid_to_vertex.right);
  }

  /** set of all UnitIDs in interaction graph */
  std::set<UID_t> get_all_uids_set() const {
    auto uids = get_all_uids();
    std::set<UID_t> out{uids.begin(), uids.end()};
    return out;
  }

  std::vector<UID_t> get_all_uids_vec() const {
    // fix UID ordering by first collecting UIDs in a set
    auto uids = get_all_uids_set();
    std::vector<UID_t> out{uids.begin(), uids.end()};
    return out;
  }

  /** iterator of all UnitIDs in interaction graph */
  auto get_all_uids() const {
    return get_vertices_it() |
           boost::adaptors::transformed(
               [this](const auto& v) { return get_uid(v); });
  }

  /* return UIDs with greatest (undirected) degree in graph */
  std::set<UID_t> max_degree_uids() const {
    std::set<UID_t> out;
    auto max_vertices = graphs::utils::max_degree_nodes(graph);
    std::transform(
        max_vertices.begin(), max_vertices.end(), std::inserter(out, out.begin()),
        [this](Vertex v) { return UID_t(get_uid(v)); });
    return out;
  }

  /** return UIDs with smallest (undirected) degree */
  std::set<UID_t> min_degree_uids() const {
    std::set<UID_t> out;
    auto min_vertices = graphs::utils::min_degree_nodes(graph);
    std::transform(
        min_vertices.begin(), min_vertices.end(), std::inserter(out, out.begin()),
        [this](Vertex v) { return UID_t(get_uid(v)); });
    return out;
  }

  /** returns path from `root` to `target` */
  std::vector<UID_t> get_path(const UID_t root, const UID_t target) const {
    if (!node_exists(root) || !node_exists(target)) {
      throw NodeDoesNotExistError(
          "Trying to get path between non-existent vertices");
    }
    using parent_vec = typename BFS<UndirectedConnGraph>::parent_vec;

    const UndirectedConnGraph& g = get_undirected_connectivity();
    auto bfs = run_bfs(to_vertices(root), g);

    auto to_uid = [&g](UndirectedVertex v) { return g[v]; };
    parent_vec path = bfs.path_to_root(to_vertices(target));
    std::vector<UID_t> converted_path(path.size());
    std::transform(path.begin(), path.end(), converted_path.begin(), to_uid);
    return converted_path;
  }

  /** get undirected adjacent UnitIDs */
  std::set<UID_t> get_neighbour_uids(const UID_t uid) const {
    if (!node_exists(uid)) {
      throw NodeDoesNotExistError(
          "Trying to get neighbours from non-existent vertex");
    }
    std::set<UID_t> neighbours;
    for (auto [it, end] = boost::out_edges(to_vertices(uid), graph); it != end;
         ++it) {
      neighbours.insert(UID_t(get_uid(boost::target(*it, graph))));
    }
    for (auto [it, end] = boost::in_edges(to_vertices(uid), graph); it != end;
         ++it) {
      neighbours.insert(UID_t(get_uid(boost::source(*it, graph))));
    }
    return neighbours;
  }

  bool operator==(const UIDConnectivityBase<UID_t>& other) const {
    std::set<UID_t> uids = this->get_all_uids_set();
    if (uids != other.get_all_uids_set()) return false;
    for (const UID_t& u : uids) {
      for (const UID_t& v : uids) {
        if (this->edge_exists(u, v)) {
          if (!other.edge_exists(u, v))
            return false;
          else if (
              this->get_connection_weight(u, v) !=
              other.get_connection_weight(u, v))
            return false;
        } else if (other.edge_exists(u, v))
          return false;
      }
    }
    return true;
  }

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
  const UID_t& get_uid(Vertex v) const { return graph[v]; }

  left_map& to_vertices() { return uid_to_vertex.left; }
  const left_map& to_vertices() const { return uid_to_vertex.left; }
  Vertex to_vertices(const UID_t& uid) const { return to_vertices().at(uid); }

  right_map& from_vertices() { return uid_to_vertex.right; }
  const right_map& from_vertices() const { return uid_to_vertex.right; }
  UID_t from_vertices(Vertex v) const { return from_vertices().at(v); }

  ConnGraph graph;
  vertex_bimap uid_to_vertex;
};

/**
 * Exception thrown because two nodes are disconnected from one another.
 *
 * @tparam UID_t node type
 */
template <typename UID_t>
class UIDsNotConnected : public std::logic_error {
 public:
  UIDsNotConnected(const UID_t& uid1, const UID_t& uid2)
      : std::logic_error(
            uid1.repr() + " and " + uid2.repr() + " are not connected") {}
};

/**
 * UIDConnectivity instances are directed graphs. It is a wrapper around a
 * BGL graph that provides a clean class API, taking care of mapping all BGL
 * vertices and edge pointers to nodes, respectively pairs of nodes.
 *
 * The vertices and edges can be given weights of type double if desired, and
 * the underlying undirected graph can be computed.
 *
 * All functionality for this class is implemented in the base class
 * UIDConnectivityBase. This class only adds caching of some function calls for
 * efficiency, invalidating cache in case of changes on the underlying graph.
 */
template <typename UID_t>
class UIDConnectivity : public UIDConnectivityBase<UID_t> {
 private:
  using Base = UIDConnectivityBase<UID_t>;

 public:
  using Base::Base;
  using Base::node_exists;
  /* useful type synonyms */
  using UndirectedConnGraph = typename Base::UndirectedConnGraph;
  using ConnGraph = typename Base::ConnGraph;
  using Connection = typename Base::Connection;
  using Vertex = typename Base::Vertex;

  // We cache distances. A value of zero in the cache implies that the nodes are
  // disconnected (unless they are equal).
  const std::vector<std::size_t>& get_distances(const UID_t& root) const& {
    if (distance_cache.find(root) == distance_cache.end()) {
      distance_cache[root] = Base::get_distances(UID_t(root));
    }
    return distance_cache[root];
  }

  std::vector<std::size_t>&& get_distances(const UID_t& root) const&& {
    if (distance_cache.find(root) == distance_cache.end()) {
      distance_cache[root] = Base::get_distances(UID_t(root));
    }
    return std::move(distance_cache[root]);
  }

  /**
   * Graph distance between two nodes.
   *
   * @param uid1 first node
   * @param uid2 second node
   *
   * @return length of shortest path between the nodes
   * @throws UIDsNotConnected if there is no path between the nodes
   */
  std::size_t get_distance(const UID_t uid1, const UID_t uid2) const {
    if (uid1 == uid2) {
      return 0;
    }
    size_t d;
    if (distance_cache.find(uid1) != distance_cache.end()) {
      d = distance_cache[uid1][this->to_vertices(uid2)];
    } else if (distance_cache.find(uid2) != distance_cache.end()) {
      d = distance_cache[uid2][this->to_vertices(uid1)];
    } else {
      distance_cache[uid1] = Base::get_distances(uid1);
      d = distance_cache[uid1][this->to_vertices(uid2)];
    }
    if (d == 0) {
      throw UIDsNotConnected(uid1, uid2);
    }
    return d;
  }

  /** Returns all nodes at a given distance from a given 'source' node */
  std::vector<UID_t> uids_at_distance(
      const UID_t& root, std::size_t distance) const {
    auto dists = get_distances(root);
    std::vector<UID_t> out;
    for (unsigned i = 0; i < dists.size(); ++i) {
      if (dists[i] == distance) {
        out.push_back(UID_t(this->get_uid(i)));
      }
    }
    return out;
  }

  // we cache the undirected graph
  const UndirectedConnGraph& get_undirected_connectivity() const& {
    if (!undir_graph) {
      undir_graph = Base::get_undirected_connectivity();
    }
    return undir_graph.value();
  }

  UndirectedConnGraph&& get_undirected_connectivity() const&& {
    if (!undir_graph) {
      undir_graph = Base::get_undirected_connectivity();
    }
    return std::move(undir_graph.value());
  }

  // these functions invalidate caching: invalidate cache then call Base
  // function
  void remove_uid(const UID_t uid) {
    invalidate_cache();
    Base::remove_uid(uid);
  }
  void add_uid(const UID_t uid) {
    invalidate_cache();
    Base::add_uid(uid);
  }

  void remove_stray_uids() {
    invalidate_cache();
    Base::remove_stray_uids();
  }

  void add_connection(const UID_t uid1, const UID_t uid2, unsigned weight = 1) {
    invalidate_cache();
    Base::add_connection(uid1, uid2, weight);
  }

  void remove_connections(const std::vector<Connection>& edges) {
    invalidate_cache();
    Base::remove_connections(edges);
  }

  void remove_connection(
      const Connection edge, bool remove_unused_vertices = false) {
    invalidate_cache();
    Base::remove_connection(edge, remove_unused_vertices);
  }

  void remove_connection(
      const UID_t uid1, const UID_t uid2, bool remove_unused_vertices = false) {
    invalidate_cache();
    Base::remove_connection(uid1, uid2, remove_unused_vertices);
  }

 private:
  inline void invalidate_cache() {
    distance_cache.clear();
    undir_graph = std::nullopt;
  }
  mutable std::map<UID_t, std::vector<std::size_t>> distance_cache;
  mutable std::optional<UndirectedConnGraph> undir_graph;
};

}  // namespace tket::graphs

#endif
