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
#include <vector>

#include "boost/iterator/transform_iterator.hpp"
#include "boost/range/adaptor/transformed.hpp"

namespace tket::graphs {

namespace detail {

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::UIDConnectivityBase(
    const std::vector<UID_t>& uids)
    : graph(ConnGraph()), fc_(false) {
  for (const UnitID& uid : uids) {
    add_uid(UID_t(uid));
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::UIDConnectivityBase(
    const std::vector<Connection>& edges)
    : graph(ConnGraph()), fc_(false) {
  for (auto [uid1, uid2] : edges) {
    if (!uid_exists(uid1)) {
      add_uid(uid1);
    }
    if (!uid_exists(uid2)) {
      add_uid(uid2);
    }
    boost::add_edge(
        to_vertices(uid1), to_vertices(uid2), std::get<ConnGraph>(graph));
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::UIDConnectivityBase(
    unsigned n)
    : graph(FullConnGraph()), fc_(true) {
  for (unsigned i = 0; i < n; i++) {
    Node n("fcNode", i);
    std::get<FullConnGraph>(graph).push_back(n);
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::UIDConnectivityBase(
    const FullConnGraph& fc)
    : graph(FullConnGraph()), fc_(true) {
  std::get<FullConnGraph>(graph) = fc;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::add_connection(
    const UID_t uid1, const UID_t uid2, unsigned weight) {
  if (is_fc()) {
    throw std::logic_error("Graph is fully connected");
  }
  if (!uid_exists(uid1) || !uid_exists(uid2)) {
    throw UIDDoesNotExistError(
        "The UIDs passed to UIDConnectivity::add_connection must "
        "exist");
  }
  boost::add_edge(
      to_vertices(uid1), to_vertices(uid2), UIDInteraction(weight),
      std::get<ConnGraph>(graph));
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::remove_connection(
    const Connection edge, bool remove_unused_vertices) {
  if (is_fc()) {
    throw std::logic_error("Graph is fully connected");
  }
  if (!uid_exists(edge.first) || !uid_exists(edge.second)) {
    throw UIDDoesNotExistError(
        "Trying to remove an edge with non-existent vertices");
  }
  auto [e, exists] = boost::edge(
      to_vertices(edge.first), to_vertices(edge.second),
      std::get<ConnGraph>(graph));
  if (!exists) {
    throw EdgeDoesNotExistError(
        "The edge (" + edge.first.repr() + ", " + edge.second.repr() +
        ")"
        "cannot be removed as it does not exist");
  }
  utils::remove_edge_with_map(
      e, std::get<ConnGraph>(graph), uid_to_vertex.right,
      remove_unused_vertices);
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
  if (is_fc()) {
    return uid1 != uid2;
  } else {
    auto [_, exists] = boost::edge(
        to_vertices(uid1), to_vertices(uid2), std::get<ConnGraph>(graph));
    return exists;
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
bool UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::uid_exists(
    const UID_t uid) const {
  if (is_fc()) {
    const FullConnGraph& fc = std::get<FullConnGraph>(graph);
    return std::find(fc.begin(), fc.end(), uid) != fc.end();
  } else {
    return to_vertices().find(uid) != to_vertices().end();
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
unsigned
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_connection_weight(
    const UID_t uid1, const UID_t uid2) const {
  if (is_fc()) {
    throw std::logic_error("Graph is fully connected");
  }
  if (!uid_exists(uid1) || !uid_exists(uid2)) {
    throw UIDDoesNotExistError(
        "Trying to retrieve edge weight from non-existent vertices");
  }
  auto [e, exists] = boost::edge(
      to_vertices(uid1), to_vertices(uid2), std::get<ConnGraph>(graph));
  if (!exists) {
    return 0.;
  }

  return std::get<ConnGraph>(graph)[e].weight;
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
unsigned UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_degree(
    const UID_t v) const {
  if (!uid_exists(v)) {
    throw UIDDoesNotExistError(
        "Trying to retrieve vertex degree from non-existent vertex");
  }
  if (is_fc()) {
    return std::get<FullConnGraph>(graph).size() - 1;
  } else {
    return boost::degree(to_vertices(v), std::get<ConnGraph>(graph));
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::set<UID_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::max_degree_uids() const {
  if (is_fc()) {
    const FullConnGraph& fc = std::get<FullConnGraph>(graph);
    return {fc.begin(), fc.end()};
  } else {
    std::set<UID_t> out;
    auto max_vertices =
        graphs::utils::max_degree_nodes(std::get<ConnGraph>(graph));
    std::transform(
        max_vertices.begin(), max_vertices.end(),
        std::inserter(out, out.begin()),
        [this](Vertex v) { return UID_t(get_uid(v)); });
    return out;
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::set<UID_t>
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::min_degree_uids() const {
  if (is_fc()) {
    const FullConnGraph& fc = std::get<FullConnGraph>(graph);
    return {fc.begin(), fc.end()};
  } else {
    std::set<UID_t> out;
    auto min_vertices =
        graphs::utils::min_degree_nodes(std::get<ConnGraph>(graph));
    std::transform(
        min_vertices.begin(), min_vertices.end(),
        std::inserter(out, out.begin()),
        [this](Vertex v) { return UID_t(get_uid(v)); });
    return out;
  }
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
  if (is_fc()) {
    std::vector<UID_t> uids{root};
    if (target != root) {
      uids.push_back(target);
    }
    return uids;
  } else {
    using parent_vec = typename BFS<UndirectedConnGraph>::parent_vec;

    const UndirectedConnGraph& g = get_undirected_connectivity();
    auto bfs = run_bfs(to_vertices(root), g);

    auto to_uid = [&g](UndirectedVertex v) { return g[v].uid; };
    parent_vec path = bfs.path_to_root(to_vertices(target));
    std::vector<UID_t> converted_path(path.size());
    std::transform(path.begin(), path.end(), converted_path.begin(), to_uid);
    return converted_path;
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
std::size_t
UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_max_depth(
    const UID_t root) const {
  if (!uid_exists(root)) {
    throw UIDDoesNotExistError("Trying to get depth from non-existent vertex");
  }
  if (is_fc()) {
    unsigned n = std::get<FullConnGraph>(graph).size();
    return n == 1 ? 0 : 1;
  } else {
    return run_bfs(to_vertices(root), get_undirected_connectivity())
        .max_depth();
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
unsigned UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::get_out_degree(
    const UID_t uid) const {
  if (!uid_exists(uid)) {
    throw UIDDoesNotExistError(
        "Trying to get outdegree from non-existent vertex");
  }
  if (is_fc()) {
    return std::get<FullConnGraph>(graph).size() - 1;
  } else {
    return boost::out_degree(to_vertices(uid), std::get<ConnGraph>(graph));
  }
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
  if (is_fc()) {
    std::list<UID_t> fc = std::get<FullConnGraph>(graph);
    std::set<UID_t> neighbours = {fc.begin(), fc.end()};
    neighbours.erase(v);
    return neighbours;
  } else {
    std::set<UID_t> neighbours;
    for (auto [it, end] =
             boost::out_edges(to_vertices(v), std::get<ConnGraph>(graph));
         it != end; ++it) {
      neighbours.insert(
          UID_t(get_uid(boost::target(*it, std::get<ConnGraph>(graph)))));
    }
    for (auto [it, end] =
             boost::in_edges(to_vertices(v), std::get<ConnGraph>(graph));
         it != end; ++it) {
      neighbours.insert(
          UID_t(get_uid(boost::source(*it, std::get<ConnGraph>(graph)))));
    }
    return neighbours;
  }
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
  if (is_fc()) {
    const FullConnGraph& fc = std::get<FullConnGraph>(graph);
    for (const UID_t& v0 : fc) {
      for (const UID_t& v1 : fc) {
        if (v0 != v1) {
          out.push_back({v0, v1});
        }
      }
    }
  } else {
    for (auto e : get_edges_it()) {
      UIDVertex<UID_t> source = std::get<ConnGraph>(
          graph)[boost::source(e, std::get<ConnGraph>(graph))];
      UIDVertex<UID_t> target = std::get<ConnGraph>(
          graph)[boost::target(e, std::get<ConnGraph>(graph))];
      out.push_back({source.uid, target.uid});
    }
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
  if (is_fc()) {
    unsigned nm1 = std::get<FullConnGraph>(graph).size() - 1;
    std::vector<std::size_t> dists(nm1);
    for (unsigned i = 0; i < nm1; i++) dists[i] = 1;
    return dists;
  } else {
    return run_bfs(to_vertices(root), get_undirected_connectivity())
        .get_dists();
  }
}

// ideally we would like to keep this private..
template <typename UID_t, typename OutEdgeListS, typename VertexListS>
auto UIDConnectivityBase<
    UID_t, OutEdgeListS, VertexListS>::get_undirected_connectivity() const
    -> UndirectedConnGraph {
  if (is_fc()) {
    // TODO We could make this work by returning a std::variant. Quite lot of
    // other code would need to be adapted, so defer till we have a use case.
    throw std::logic_error("Graph is fully connected");
  } else {
    return graphs::utils::symmetrise<UndirectedConnGraph>(
        std::get<ConnGraph>(graph));
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
void UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::to_fc() {
  if (!is_fc()) {
    auto uids = get_all_uids();
    FullConnGraph fc{uids.begin(), uids.end()};
    graph = fc;
    fc_ = true;
  }
}

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
bool UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>::operator==(
    const UIDConnectivityBase<UID_t, OutEdgeListS, VertexListS>& other) const {
  if (is_fc() != other.is_fc()) return false;
  std::set<UID_t> uids = this->get_all_uids_set();
  if (uids != other.get_all_uids_set()) return false;
  if (is_fc()) return true;
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
  } else if (is_fc()) {
    return 1;
  } else {
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

template <typename UID_t, typename OutEdgeListS, typename VertexListS>
bool UIDConnectivity<UID_t, OutEdgeListS, VertexListS>::is_fc() const {
  return Base::is_fc();
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

}  // namespace tket::graphs
