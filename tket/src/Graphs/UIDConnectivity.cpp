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

namespace tket::graphs {

template <typename UID_t>
const std::vector<std::size_t>& UIDConnectivity<UID_t>::get_distances(
    const UID_t& root) const& {
  if (distance_cache.find(root) == distance_cache.end()) {
    distance_cache[root] = Base::get_distances(UID_t(root));
  }
  return distance_cache[root];
}

template <typename UID_t>
std::vector<std::size_t>&& UIDConnectivity<UID_t>::get_distances(
    const UID_t& root) const&& {
  if (distance_cache.find(root) == distance_cache.end()) {
    distance_cache[root] = Base::get_distances(UID_t(root));
  }
  return std::move(distance_cache[root]);
}

// this function uses caching indirectly, since it calls get_distances,
// so we cannot define it in Base class (unless we made everything virtual)
template <typename UID_t>
std::vector<UID_t> UIDConnectivity<UID_t>::uids_at_distance(
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

template <typename UID_t>
std::size_t UIDConnectivity<UID_t>::get_distance(
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

template <typename UID_t>
auto UIDConnectivity<UID_t>::get_undirected_connectivity() const& -> const
    UndirectedConnGraph& {
  if (!undir_graph) {
    undir_graph = Base::get_undirected_connectivity();
  }
  return undir_graph.value();
}

template <typename UID_t>
auto UIDConnectivity<UID_t>::get_undirected_connectivity()
    const&& -> UndirectedConnGraph&& {
  if (!undir_graph) {
    undir_graph = Base::get_undirected_connectivity();
  }
  return std::move(undir_graph.value());
}

template <typename UID_t>
void UIDConnectivity<UID_t>::add_uid(const UID_t uid) {
  invalidate_cache();
  Base::add_uid(uid);
}
template <typename UID_t>
void UIDConnectivity<UID_t>::remove_uid(const UID_t uid) {
  invalidate_cache();
  Base::remove_uid(uid);
}
template <typename UID_t>
void UIDConnectivity<UID_t>::remove_stray_uids() {
  invalidate_cache();
  Base::remove_stray_uids();
}
template <typename UID_t>
void UIDConnectivity<UID_t>::add_connection(
    const UID_t uid1, const UID_t uid2, unsigned val) {
  invalidate_cache();
  Base::add_connection(uid1, uid2, val);
}
template <typename UID_t>
void UIDConnectivity<UID_t>::remove_connections(
    const std::vector<Connection>& edges) {
  invalidate_cache();
  Base::remove_connections(edges);
}
template <typename UID_t>
void UIDConnectivity<UID_t>::remove_connection(
    const Connection edge, bool remove_unused_vertices) {
  invalidate_cache();
  Base::remove_connection(edge, remove_unused_vertices);
}
template <typename UID_t>
void UIDConnectivity<UID_t>::remove_connection(
    const UID_t uid1, const UID_t uid2, bool remove_unused_vertices) {
  invalidate_cache();
  Base::remove_connection(uid1, uid2, remove_unused_vertices);
}

template class UIDConnectivity<Node>;
template class UIDConnectivity<Qubit>;

}  // namespace tket::graphs
