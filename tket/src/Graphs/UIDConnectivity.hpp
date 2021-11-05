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
#include <boost/range/iterator_range_core.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <type_traits>

#include "Graphs/TreeSearch.hpp"
#include "Graphs/UIDConnectivityBase.hpp"
#include "Graphs/Utils.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/UnitID.hpp"

namespace tket::graphs {

/**
 * UIDConnectivity instances are graphs of UnitID vertices.
 * It should be instantiated with UnitIDs, or one of its subtypes `Qubit` or
 * `Node`
 *
 * All functionality for this class is implemented in the base class
 * UIDConnectivityBase. This class only adds caching of some function calls for
 * efficiency, innvalidating cache in case of changes on the underlying graph
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

// template explicit instations, with implementations in cpp file
extern template class UIDConnectivity<UnitID>;
extern template class UIDConnectivity<Node>;
extern template class UIDConnectivity<Qubit>;

}  // namespace tket::graphs

#endif
