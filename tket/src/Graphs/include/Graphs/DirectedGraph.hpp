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
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

#include "Graphs/AbstractGraph.hpp"
#include "Graphs/TreeSearch.hpp"
#include "Graphs/Utils.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket::graphs {

/**
 * Exception thrown because two nodes are disconnected from one another.
 *
 * @tparam T node type
 */
template <typename T>
class NodesNotConnected : public std::logic_error {
 public:
  NodesNotConnected(const T& node1, const T& node2)
      : std::logic_error(
            node1.repr() + " and " + node2.repr() + " are not connected") {}
};

/** Weighted edge */
struct WeightedEdge {
  explicit WeightedEdge(unsigned w = 1) : weight(w) {}
  unsigned weight;
};

template <typename T>
class DirectedGraphBase : public AbstractGraph<T> {
 protected:
  using ConnGraph = boost::adjacency_list<
      boost::vecS, boost::vecS, boost::bidirectionalS, T, WeightedEdge>;
  using UndirectedConnGraph = boost::adjacency_list<
      boost::setS, boost::vecS, boost::undirectedS, T, WeightedEdge>;
  using Vertex = utils::vertex<ConnGraph>;
  using UndirectedVertex = utils::vertex<UndirectedConnGraph>;
  using Connection = typename AbstractGraph<T>::Edge;
  using AbstractGraph<T>::nodes_;

 public:
  using AbstractGraph<T>::node_exists;
  using AbstractGraph<T>::edge_exists;
  using AbstractGraph<T>::get_all_nodes_vec;
  using AbstractGraph<T>::get_all_edges_vec;
  using AbstractGraph<T>::n_nodes;

  /** Construct an empty graph. */
  DirectedGraphBase() : AbstractGraph<T>(), graph(), node_to_vertex() {}

  /** Construct a graph with no edges. */
  explicit DirectedGraphBase(const std::vector<T>& nodes)
      : AbstractGraph<T>(nodes), graph() {
    for (const T& node : nodes) {
      add_node(node);
    }
  }

  /** Construct a graph from its edges. */
  explicit DirectedGraphBase(const std::vector<Connection>& edges) : graph() {
    for (auto [node1, node2] : edges) {
      if (!node_exists(node1)) {
        add_node(node1);
      }
      if (!node_exists(node2)) {
        add_node(node2);
      }
      add_connection(node1, node2);
    }
  }

  /** Test whether two nodes are connected. */
  bool edge_exists(const T& node1, const T& node2) const override {
    if (!node_exists(node1) || !node_exists(node2)) {
      throw NodeDoesNotExistError(
          "The nodes passed to DirectedGraph::edge_exists must exist");
    }
    auto [_, exists] =
        boost::edge(to_vertices(node1), to_vertices(node2), graph);
    return exists;
  }

  /** Add a new node to the graph. */
  void add_node(const T& node) {
    nodes_.insert(node);
    Vertex v = boost::add_vertex(node, graph);
    node_to_vertex.left.insert({node, v});
  }

  /** Remove a node from the graph. */
  void remove_node(const T& node) {
    if (!node_exists(node)) {
      throw NodeDoesNotExistError(
          "The node passed to DirectedGraph::remove_node must exist!");
    }
    nodes_.erase(node);
    Vertex v = to_vertices(node);
    boost::clear_vertex(v, graph);
    utils::remove_vertex_with_map(v, graph, node_to_vertex.right);
  }

  /** Add a weighted edge to the graph. */
  void add_connection(const T& node1, const T& node2, unsigned weight = 1) {
    if (!node_exists(node1) || !node_exists(node2)) {
      throw NodeDoesNotExistError(
          "The nodes passed to DirectedGraph::add_connection must exist");
    }
    boost::add_edge(
        to_vertices(node1), to_vertices(node2), WeightedEdge(weight), graph);
  }

  /** Remove an edge from the graph. */
  void remove_connection(const Connection& edge) {
    if (!node_exists(edge.first) || !node_exists(edge.second)) {
      throw NodeDoesNotExistError(
          "Trying to remove an edge with non-existent vertices");
    }
    auto [e, exists] =
        boost::edge(to_vertices(edge.first), to_vertices(edge.second), graph);
    if (!exists) {
      throw EdgeDoesNotExistError(
          "The edge (" + edge.first.repr() + ", " + edge.second.repr() +
          ") cannot be removed as it does not exist");
    }
    utils::remove_edge_with_map(e, graph, node_to_vertex.right);
  }

  /** Remove an edge between two nodes. */
  void remove_connection(const T& node1, const T& node2) {
    remove_connection({node1, node2});
  }

  /**
   * Get the weight of an edge between two nodes.
   *
   * Returns 0 if the edge does not exist.
   */
  unsigned get_connection_weight(const T& node1, const T& node2) const {
    if (!node_exists(node1) || !node_exists(node2)) {
      throw NodeDoesNotExistError(
          "Trying to retrieve edge weight from non-existent vertices");
    }
    auto [e, exists] =
        boost::edge(to_vertices(node1), to_vertices(node2), graph);
    if (!exists) {
      return 0;
    }

    return graph[e].weight;
  }

  /** Get the degree of a node. */
  unsigned get_degree(const T& node) const {
    if (!node_exists(node)) {
      throw NodeDoesNotExistError(
          "Trying to retrieve vertex degree from non-existent vertex");
    }
    return boost::degree(to_vertices(node), graph);
  }

  /** Get the out-degree of a node. */
  unsigned get_out_degree(const T& node) const {
    if (!node_exists(node)) {
      throw NodeDoesNotExistError(
          "Trying to get outdegree from non-existent vertex");
    }
    return boost::out_degree(to_vertices(node), graph);
  }

  /** Number of edges. */
  inline unsigned n_connections() const { return boost::num_edges(graph); }

  /** Number of nodes with degree > 0. */
  inline unsigned n_connected() const {
    auto [beg, end] = boost::vertices(graph);
    auto nonzero_deg = [this](Vertex v) { return boost::degree(v, graph) > 0; };
    return std::count_if(beg, end, nonzero_deg);
  }

  /** All edges in the graph. */
  std::set<Connection> edges() const {
    std::vector<Connection> vec = get_all_edges_vec();
    return std::set(vec.begin(), vec.end());
  }

  /** All edges as a vector. */
  std::vector<Connection> get_all_edges_vec() const override {
    std::vector<Connection> out;
    for (auto e : get_edges_it()) {
      T source = graph[boost::source(e, graph)];
      T target = graph[boost::target(e, graph)];
      out.push_back({source, target});
    }
    return out;
  }

  /** Return an unweighted undirected graph with the same connectivity. */
  UndirectedConnGraph get_undirected_connectivity() const {
    return graphs::utils::symmetrise<UndirectedConnGraph>(graph);
  }

  /** Get all distances between pairs of nodes. */
  std::vector<std::size_t> get_distances(const T& root) const {
    if (!node_exists(root)) {
      throw NodeDoesNotExistError(
          "Trying to get distances from non-existent root vertex");
    }
    return run_bfs(to_vertices(root), get_undirected_connectivity())
        .get_dists();
  }

  /** Get the distance between two nodes. */
  unsigned get_distance(const T& node1, const T& node2) const override {
    if (node1 == node2) {
      return 0;
    }
    size_t d = get_distances(node1)[to_vertices(node2)];
    if (d == 0) {
      throw NodesNotConnected(node1, node2);
    }
    return d;
  }

  /** Remove nodes with degree 0 */
  inline void remove_stray_nodes() {
    std::set<T> strays;
    for (const T& node : nodes_) {
      if (get_degree(node) == 0) {
        strays.insert(node);
      }
    }
    for (const T& node : strays) {
      remove_node(node);
    }
  }

  /* Return nodes with maximum degree. */
  std::set<T> max_degree_nodes() const {
    std::set<T> out;
    auto max_vertices = graphs::utils::max_degree_nodes(graph);
    std::transform(
        max_vertices.begin(), max_vertices.end(),
        std::inserter(out, out.begin()),
        [this](Vertex v) { return get_node(v); });
    return out;
  }

  /** Return nodes with minimum degree. */
  std::set<T> min_degree_nodes() const {
    std::set<T> out;
    auto min_vertices = graphs::utils::min_degree_nodes(graph);
    std::transform(
        min_vertices.begin(), min_vertices.end(),
        std::inserter(out, out.begin()),
        [this](Vertex v) { return get_node(v); });
    return out;
  }

  /** Returns path between two nodes. */
  std::vector<T> get_path(const T& root, const T& target) const {
    if (!node_exists(root) || !node_exists(target)) {
      throw NodeDoesNotExistError(
          "Trying to get path between non-existent vertices");
    }
    using parent_vec = typename BFS<UndirectedConnGraph>::parent_vec;

    const UndirectedConnGraph& g = get_undirected_connectivity();
    auto bfs = run_bfs(to_vertices(root), g);

    auto to_node = [&g](UndirectedVertex v) { return g[v]; };
    parent_vec path = bfs.path_to_root(to_vertices(target));
    std::vector<T> converted_path(path.size());
    std::transform(path.begin(), path.end(), converted_path.begin(), to_node);
    return converted_path;
  }

  /** Get all neighbours of a node. */
  std::set<T> get_neighbour_nodes(const T& node) const {
    if (!node_exists(node)) {
      throw NodeDoesNotExistError(
          "Trying to get neighbours from non-existent vertex");
    }
    std::set<T> neighbours;
    for (auto [it, end] = boost::out_edges(to_vertices(node), graph); it != end;
         ++it) {
      neighbours.insert(get_node(boost::target(*it, graph)));
    }
    for (auto [it, end] = boost::in_edges(to_vertices(node), graph); it != end;
         ++it) {
      neighbours.insert(get_node(boost::source(*it, graph)));
    }
    return neighbours;
  }

  bool operator==(const DirectedGraphBase<T>& other) const {
    if (nodes_ != other.nodes_) return false;
    for (const T& u : nodes_) {
      for (const T& v : nodes_) {
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

  using vertex_bimap = boost::bimap<T, Vertex>;
  using left_map = typename vertex_bimap::left_map;
  using right_map = typename vertex_bimap::right_map;
  /* get node from vertex */
  const T& get_node(Vertex v) const { return graph[v]; }

  left_map& to_vertices() { return node_to_vertex.left; }
  const left_map& to_vertices() const { return node_to_vertex.left; }
  Vertex to_vertices(const T& node) const { return to_vertices().at(node); }

  right_map& from_vertices() { return node_to_vertex.right; }
  const right_map& from_vertices() const { return node_to_vertex.right; }
  T from_vertices(Vertex v) const { return from_vertices().at(v); }

  ConnGraph graph;
  vertex_bimap node_to_vertex;
};

/**
 * DirectedGraph instances are directed graphs. It is a wrapper around a
 * BGL graph that provides a clean class API, taking care of mapping all BGL
 * vertices and edge pointers to nodes, respectively pairs of nodes.
 *
 * The vertices and edges can be given weights of type double if desired, and
 * the underlying undirected graph can be computed.
 *
 * All functionality for this class is implemented in the base class
 * DirectedGraphBase. This class only adds caching of some function calls for
 * efficiency, invalidating cache in case of changes on the underlying graph.
 */
template <typename T>
class DirectedGraph : public DirectedGraphBase<T> {
 private:
  using Base = DirectedGraphBase<T>;

 public:
  using Base::Base;
  using Base::get_all_nodes_vec;
  using Base::n_nodes;
  using Base::node_exists;
  /* useful type synonyms */
  using UndirectedConnGraph = typename Base::UndirectedConnGraph;
  using ConnGraph = typename Base::ConnGraph;
  using Connection = typename Base::Connection;
  using Vertex = typename Base::Vertex;

  /**
   * Get all distances between nodes.
   */
  const std::vector<std::size_t>& get_distances(const T& root) const& {
    // We cache distances. A value of zero in the cache implies that the nodes
    // are disconnected (unless they are equal).
    if (distance_cache.find(root) == distance_cache.end()) {
      distance_cache[root] = Base::get_distances(root);
    }
    return distance_cache[root];
  }

  /**
   * Get all distances between nodes.
   */
  std::vector<std::size_t>&& get_distances(const T& root) const&& {
    if (distance_cache.find(root) == distance_cache.end()) {
      distance_cache[root] = Base::get_distances(root);
    }
    return std::move(distance_cache[root]);
  }

  unsigned get_distance(const T& node1, const T& node2) const override {
    if (node1 == node2) {
      return 0;
    }
    size_t d;
    if (distance_cache.find(node1) != distance_cache.end()) {
      d = distance_cache[node1][this->to_vertices(node2)];
    } else if (distance_cache.find(node2) != distance_cache.end()) {
      d = distance_cache[node2][this->to_vertices(node1)];
    } else {
      distance_cache[node1] = Base::get_distances(node1);
      d = distance_cache[node1][this->to_vertices(node2)];
    }
    if (d == 0) {
      throw NodesNotConnected(node1, node2);
    }
    return d;
  }

  unsigned get_diameter() override {
    unsigned N = n_nodes();
    if (N == 0) {
      throw std::logic_error("Graph is empty.");
    }
    if (!this->diameter_) {
      this->diameter_ = 0;
      const std::vector<T> nodes = get_all_nodes_vec();
      for (unsigned i = 0; i < N; i++) {
        for (unsigned j = i + 1; j < N; j++) {
          unsigned d = get_distance(nodes[i], nodes[j]);
          if (d > *this->diameter_) this->diameter_ = d;
        }
      }
    }
    return *this->diameter_;
  }

  /** Returns all nodes at a given distance from a given 'source' node */
  std::vector<T> nodes_at_distance(const T& root, std::size_t distance) const {
    auto dists = get_distances(root);
    std::vector<T> out;
    for (unsigned i = 0; i < dists.size(); ++i) {
      if (dists[i] == distance) {
        out.push_back(this->get_node(i));
      }
    }
    return out;
  }

  /** Return an unweighted undirected graph with the same connectivity. */
  const UndirectedConnGraph& get_undirected_connectivity() const& {
    // we cache the undirected graph
    if (!undir_graph) {
      undir_graph = Base::get_undirected_connectivity();
    }
    return undir_graph.value();
  }

  /** Return an unweighted undirected graph with the same connectivity. */
  UndirectedConnGraph&& get_undirected_connectivity() const&& {
    if (!undir_graph) {
      undir_graph = Base::get_undirected_connectivity();
    }
    return std::move(undir_graph.value());
  }

  // The following functions invalidate caching.

  /** Remove a node from the graph. */
  void remove_node(const T& node) {
    invalidate_cache();
    Base::remove_node(node);
  }

  /** Add a node to the graph. */
  void add_node(const T& node) {
    invalidate_cache();
    Base::add_node(node);
  }

  /** Remove nodes with degree 0. */
  void remove_stray_nodes() {
    invalidate_cache();
    Base::remove_stray_nodes();
  }

  /** Add a (weighted) edge between two nodes. */
  void add_connection(const T& node1, const T& node2, unsigned weight = 1) {
    invalidate_cache();
    Base::add_connection(node1, node2, weight);
  }

  /** Remove an edge. */
  void remove_connection(const Connection& edge) {
    invalidate_cache();
    Base::remove_connection(edge);
  }

  /** Remove a connection between two nodes. */
  void remove_connection(const T& node1, const T& node2) {
    invalidate_cache();
    Base::remove_connection(node1, node2);
  }

 private:
  inline void invalidate_cache() {
    distance_cache.clear();
    undir_graph = std::nullopt;
  }
  mutable std::map<T, std::vector<std::size_t>> distance_cache;
  mutable std::optional<UndirectedConnGraph> undir_graph;
};

}  // namespace tket::graphs
