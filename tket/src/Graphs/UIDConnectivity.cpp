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

#include "UIDConnectivity.hpp"

#include <stdexcept>

#include "boost/iterator/transform_iterator.hpp"
#include "boost/range/adaptor/transformed.hpp"

namespace tket::graph {

namespace detail {

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::UIDConnectivityBase(
    const std::vector<UID_t>& uids)
    : graph() {
  for (const UnitID& uid : uids) {
    add_uid(UID_t(uid));
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::UIDConnectivityBase(
    const std::vector<Connection>& edges)
    : graph() {
  for (auto [uid1, uid2] : edges) {
    if (!uid_exists(uid1)) {
      add_uid(uid1);
    }
    if (!uid_exists(uid2)) {
      add_uid(uid2);
    }
    boost::add_edge(to_vertices(uid1), to_vertices(uid2), graph);
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::add_connection(
    const UID_t uid1, const UID_t uid2, unsigned weight) {
  if (!uid_exists(uid1) || !uid_exists(uid2)) {
    throw UIDDoesNotExistError(
        "The UIDs passed to UIDConnectivity::add_connection must "
        "exist");
  }
  boost::add_edge(
      to_vertices(uid1), to_vertices(uid2), UIDInteraction(weight), graph);
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::remove_connection(
    const Connection edge, bool remove_unused_vertices) {
  if (!uid_exists(edge.first) || !uid_exists(edge.second)) {
    throw UIDDoesNotExistError(
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

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::remove_connection(
    const UID_t uid1, const UID_t uid2, bool remove_unused_vertices) {
  remove_connection({uid1, uid2}, remove_unused_vertices);
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::remove_connections(
    const std::vector<Connection>& edges, bool remove_unused_vertices) {
  for (const auto& e : edges) {
    remove_connection(e, remove_unused_vertices);
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
bool UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::connection_exists(
    const UID_t uid1, const UID_t uid2) const {
  if (!uid_exists(uid1) || !uid_exists(uid2)) {
    throw UIDDoesNotExistError(
        "The UIDs passed to UIDConnectivity::connection_exists must "
        "exist");
  }
  auto [_, exists] = boost::edge(to_vertices(uid1), to_vertices(uid2), graph);
  return exists;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
bool UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::uid_exists(
    const UID_t uid) const {
  return to_vertices().find(uid) != to_vertices().end();
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
unsigned
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_connection_weight(
    const UID_t uid1, const UID_t uid2) const {
  if (!uid_exists(uid1) || !uid_exists(uid2)) {
    throw UIDDoesNotExistError(
        "Trying to retrieve edge weight from non-existent vertices");
  }
  auto [e, exists] = boost::edge(to_vertices(uid1), to_vertices(uid2), graph);
  if (!exists) {
    return 0.;
  }

  return graph[e].weight;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
unsigned UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_degree(
    const UID_t v) const {
  if (!uid_exists(v)) {
    throw UIDDoesNotExistError(
        "Trying to retrieve vertex degree from non-existent vertex");
  }
  return boost::degree(to_vertices(v), graph);
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::set<UID_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::max_degree_uids() const {
  std::set<UID_t> out;
  auto max_vertices = graph::utils::max_degree_nodes(graph);
  std::transform(
      max_vertices.begin(), max_vertices.end(), std::inserter(out, out.begin()),
      [this](Vertex v) { return UID_t(get_uid(v)); });
  return out;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::set<UID_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::min_degree_uids() const {
  std::set<UID_t> out;
  auto min_vertices = graph::utils::min_degree_nodes(graph);
  std::transform(
      min_vertices.begin(), min_vertices.end(), std::inserter(out, out.begin()),
      [this](Vertex v) { return UID_t(get_uid(v)); });
  return out;
}

// Returns the path between the nodes |root| and |target|.
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::vector<UID_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_path(
    const UID_t root, const UID_t target) const {
  if (!uid_exists(root) || !uid_exists(target)) {
    throw UIDDoesNotExistError(
        "Trying to get path between non-existent vertices");
  }
  using parent_vec = typename BFS<UndirectedConnGraph>::parent_vec;

  const UndirectedConnGraph& g = get_undirected_connectivity();
  auto bfs = run_bfs(to_vertices(root), g);

  auto to_uid = [&g](UndirectedVertex v) { return g[v].uid; };
  parent_vec path = bfs.path_to_root(to_vertices(target));
  std::vector<UID_t> converted_path(path.size());
  std::transform(path.begin(), path.end(), converted_path.begin(), to_uid);
  return converted_path;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::size_t
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_max_depth(
    const UID_t root) const {
  if (!uid_exists(root)) {
    throw UIDDoesNotExistError("Trying to get depth from non-existent vertex");
  }
  return run_bfs(to_vertices(root), get_undirected_connectivity()).max_depth();
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
unsigned UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_out_degree(
    const UID_t uid) const {
  if (!uid_exists(uid)) {
    throw UIDDoesNotExistError(
        "Trying to get outdegree from non-existent vertex");
  }
  return boost::out_degree(to_vertices(uid), graph);
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::set<UID_t> UIDConnectivityBase<
    UID_t, OutEdgeListS, VertexListS>::get_all_uids_set() const {
  auto uids = get_all_uids();
  std::set<UID_t> out{uids.begin(), uids.end()};
  return out;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::vector<UID_t> UIDConnectivityBase<
    UID_t, OutEdgeListS, VertexListS>::get_all_uids_vec() const {
  // fix UID ordering by first collecting UIDs in a set
  auto uids = get_all_uids_set();
  std::vector<UID_t> out{uids.begin(), uids.end()};
  return out;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::set<UID_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_neighbour_uids(
    const UID_t v) const {
  if (!uid_exists(v)) {
    throw UIDDoesNotExistError(
        "Trying to get neighbours from non-existent vertex");
  }
  std::set<UID_t> neighbours;
  for (auto [it, end] = boost::out_edges(to_vertices(v), graph); it != end;
       ++it) {
    neighbours.insert(UID_t(get_uid(boost::target(*it, graph))));
  }
  for (auto [it, end] = boost::in_edges(to_vertices(v), graph); it != end;
       ++it) {
    neighbours.insert(UID_t(get_uid(boost::source(*it, graph))));
  }
  return neighbours;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
auto UIDConnectivityBase<
    UID_t, OutEdgeListS, VertexListS>::get_connections_set() const
    -> std::set<Connection> {
  std::vector<Connection> vec = get_connections_vec();
  return std::set(vec.begin(), vec.end());
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
auto UIDConnectivityBase<
    UID_t, OutEdgeListS, VertexListS>::get_connections_vec() const
    -> std::vector<Connection> {
  std::vector<Connection> out;
  for (auto e : get_edges_it()) {
    UIDVertex<UID_t> source = graph[boost::source(e, graph)];
    UIDVertex<UID_t> target = graph[boost::target(e, graph)];
    out.push_back({source.uid, target.uid});
  }
  return out;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::vector<std::size_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_distances(
    const UID_t root) const {
  if (!uid_exists(root)) {
    throw UIDDoesNotExistError(
        "Trying to get distances from non-existent root vertex");
  }
  return run_bfs(to_vertices(root), get_undirected_connectivity()).get_dists();
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::size_t UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_distance(
    const UID_t uid1, const UID_t uid2) const {
  if (uid1 == uid2) {
    return 0;
  }
  size_t d = get_distances(uid1)[to_vertices(uid2)];
  if (d == 0) {
    throw UIDsNotConnected(uid1, uid2);
  }
  return d;
}

// ideally we would like to keep this private..
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
auto UIDConnectivityBase<
    UID_t, OutEdgeListS, VertexListS>::get_undirected_connectivity() const
    -> UndirectedConnGraph {
  return graph::utils::symmetrise<UndirectedConnGraph>(graph);
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
bool UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::operator==(
    const UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>& other) const {
  std::set<UID_t> uids = this->get_all_uids_set();
  if (uids != other.get_all_uids_set()) return false;
  for (const UID_t& u : uids) {
    for (const UID_t& v : uids) {
      if (this->connection_exists(u, v)) {
        if (!other.connection_exists(u, v))
          return false;
        else if (
            this->get_connection_weight(u, v) !=
            other.get_connection_weight(u, v))
          return false;
      } else if (other.connection_exists(u, v))
        return false;
    }
  }
  return true;
}

}  // namespace detail

/////////////////////////////////////////////
//             UIDConnectivity
////////////////////////////////////////////
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
const std::vector<std::size_t>&
UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::get_distances(
    const UnitID& root) const& {
  if (distance_cache.find(root) == distance_cache.end()) {
    distance_cache[root] = Base::get_distances(UID_t(root));
  }
  return distance_cache[root];
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::vector<std::size_t>&&
UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::get_distances(
    const UnitID& root) const&& {
  if (distance_cache.find(root) == distance_cache.end()) {
    distance_cache[root] = Base::get_distances(UID_t(root));
  }
  return std::move(distance_cache[root]);
}

// this function uses caching indirectly, since it calls get_distances,
// so we cannot define it in Base class (unless we made everything virtual)
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::vector<UID_t>
UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::uids_at_distance(
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

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::size_t UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::get_distance(
    const UID_t uid1, const UID_t uid2) const {
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

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
auto UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::
    get_undirected_connectivity() const& -> const UndirectedConnGraph& {
  if (!undir_graph) {
    undir_graph = Base::get_undirected_connectivity();
  }
  return undir_graph.value();
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
auto UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::
    get_undirected_connectivity() const&& -> UndirectedConnGraph&& {
  if (!undir_graph) {
    undir_graph = Base::get_undirected_connectivity();
  }
  return std::move(undir_graph.value());
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::add_uid(
    const UID_t uid) {
  invalidate_cache();
  Base::add_uid(uid);
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::remove_uid(
    const UID_t uid) {
  invalidate_cache();
  Base::remove_uid(uid);
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::remove_stray_uids() {
  invalidate_cache();
  Base::remove_stray_uids();
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::add_connection(
    const UID_t uid1, const UID_t uid2, unsigned val) {
  invalidate_cache();
  Base::add_connection(uid1, uid2, val);
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::remove_connections(
    const std::vector<Connection>& edges) {
  invalidate_cache();
  Base::remove_connections(edges);
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::remove_connection(
    const Connection edge, bool remove_unused_vertices) {
  invalidate_cache();
  Base::remove_connection(edge, remove_unused_vertices);
}
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::remove_connection(
    const UID_t uid1, const UID_t uid2, bool remove_unused_vertices) {
  invalidate_cache();
  Base::remove_connection(uid1, uid2, remove_unused_vertices);
}

template struct detail::UIDVertex<UnitID>;
template struct detail::UIDVertex<Node>;
template struct detail::UIDVertex<Qubit>;
template class detail::UIDConnectivityBase<UnitID, boost::vecS, boost::vecS>;
template class detail::UIDConnectivityBase<Node, boost::vecS, boost::vecS>;
template class detail::UIDConnectivityBase<Qubit, boost::vecS, boost::vecS>;
template class UIDConnectivity<UnitID>;
template class UIDConnectivity<Node>;
template class UIDConnectivity<Qubit>;

}  // namespace tket::graph
