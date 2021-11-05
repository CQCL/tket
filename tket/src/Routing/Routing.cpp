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

#include "Routing.hpp"

#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "Utils/HelperFunctions.hpp"
#include "Utils/Json.hpp"

namespace tket {

bool RoutingConfig::operator==(const RoutingConfig& other) const {
  return (this->depth_limit == other.depth_limit) &&
         (this->distrib_limit == other.distrib_limit) &&
         (this->interactions_limit == other.interactions_limit) &&
         (this->distrib_exponent == other.distrib_exponent);
}

// If unit map is same pre and both routing, then the same placement procedure
// has happened in both cases, and routing is deterministic (!!) so same
// SWAP/Bridges added assuming same config
bool Routing::circuit_modified() const {
  if (route_stats.swap_count > 0) return true;
  if (route_stats.bridge_count > 0) return true;
  if (circ_.boundary != original_boundary) return true;
  return false;
}

/* Class Constructor */
Routing::Routing(const Circuit& _circ, const Architecture& _arc)
    : circ_(_circ), slice_frontier_(circ_), original_arc_(_arc) {
  circ_.unit_bimaps_ = _circ.unit_bimaps_;
  original_boundary = circ_.boundary;

  current_arc_ = original_arc_;
  // Checks for circuit and architecture compatibility
  if (circ_.n_qubits() > current_arc_.n_nodes() || current_arc_.n_nodes() < 1) {
    throw ArchitectureMismatch(circ_.n_qubits(), current_arc_.n_nodes());
  }

  // Information for placement & running routing with subgraph of architecture
  // Initial nodes number

  // Track which nodes are actually active
  for (const UnitID& uid : current_arc_.get_all_nodes()) {
    Node n(uid);
    interaction.insert({n, n});
  }
}

void to_json(nlohmann::json& j, const RoutingConfig& config) {
  j["depth_limit"] = config.depth_limit;
  j["distrib_limit"] = config.distrib_limit;
  j["interactions_limit"] = config.interactions_limit;
  j["distrib_exponent"] = config.distrib_exponent;
}

void from_json(const nlohmann::json& j, RoutingConfig& config) {
  config.depth_limit = j.at("depth_limit").get<unsigned>();
  config.distrib_limit = j.at("distrib_limit").get<unsigned>();
  config.interactions_limit = j.at("interactions_limit").get<unsigned>();
  config.distrib_exponent = j.at("distrib_exponent").get<double>();
}

std::vector<Node> Routing::get_active_nodes() const {
  node_vector_t ret;
  ret.reserve(qmap.size());
  for (auto [qb, n] : qmap.left) {
    ret.push_back(n);
  }
  return ret;
}

qubit_mapping_t Routing::return_final_map() const {
  return bimap_to_map(final_map.left);
}

qubit_mapping_t Routing::return_initial_map() const {
  return bimap_to_map(init_map.left);
}

bool subgraph_remove_if_connected(
    Architecture& arc, const Architecture& subarc, const Node& node) {
  // do not remove if node is in subarc
  if (subarc.node_exists(node)) {
    return false;
  }
  if (subarc.n_nodes() > 0) {
    node_set_t ap = arc.get_articulation_points(subarc);

    if (ap.find(node) != ap.end()) {
      return false;
    }
  }

  arc.remove_node(node);
  return true;
}

void remove_unmapped_nodes(
    Architecture& arc, qubit_bimap_t& map, Circuit& circ) {
  std::vector<Node> unmapped_nodes;
  std::vector<Node> mapped_nodes;

  r_const_iterator_t iend = map.right.end();
  for (const UnitID& uid : arc.get_all_nodes()) {
    Node n(uid);
    r_const_iterator_t find_node = map.right.find(n);
    if (find_node == iend) {
      unmapped_nodes.push_back(n);
    } else {
      mapped_nodes.push_back(n);
    }
  }
  Architecture subarc = arc.create_subarch(mapped_nodes);

  // sort mapped nodes from least connected to most (remove least connected
  // first)
  std::sort(
      unmapped_nodes.begin(), unmapped_nodes.end(), [&arc](Node x, Node y) {
        return (arc.get_out_degree(x) < arc.get_out_degree(y));
      });

  qubit_vector_t available;
  for (const Qubit& q : circ.all_qubits()) {
    if (map.left.find(q) == map.left.end()) {
      available.push_back(q);
    }
  }

  for (const Node& node : unmapped_nodes) {
    if (!subgraph_remove_if_connected(arc, subarc, node)) {
      // if node can't be removed, map to first unmapped qubit
      if (available.empty())
        throw CircuitInvalidity(
            "Routing is unable to construct connected placement from partial "
            "placement using unplaced logical qubits. Please update the "
            "circuit placement to a set of connected physical qubits.");
      map.insert({available.front(), node});
      available.erase(available.begin());
    }
  }
}

qubit_mapping_t get_qmap_from_circuit(Architecture& arc, Circuit& circ) {
  qubit_vector_t all_qbs = circ.all_qubits();
  node_set_t all_nodes = arc.get_all_nodes_set();

  qubit_mapping_t qubit_map;
  for (Qubit q : all_qbs) {
    Node n(q);
    if (all_nodes.find(n) != all_nodes.end()) {
      qubit_map.insert({q, n});
    }
  }
  return qubit_map;
}

std::pair<Circuit, bool> Routing::solve(const RoutingConfig& config) {
  config_ = config;
  qubit_mapping_t qubit_map = get_qmap_from_circuit(current_arc_, circ_);
  slice_frontier_.init();
  if (slice_frontier_.slice->empty()) {
    organise_registers_and_maps();
  } else {
    // Some nodes are permanently unused due to difference in architecture nodes
    // and number of used wires in circuit To account for this, place highest
    // numbered wires (i.e. unused) into set bad nodes of architecture

    // Placement method attempts to find a good initial allocation of qubits to
    // nodes, aiming to reduce overall circuit depth. The method aims to put
    // intreacting qubits in the first few circuit timesteps on adjacent nodes
    // If no placement, qubits placed sequentially on nodes i.e. qubit 0 -> node
    // 0 etc.

    if (qubit_map.size() != 0) {
      init_map.left.insert(qubit_map.begin(), qubit_map.end());
    }
    remove_unmapped_nodes(current_arc_, init_map, circ_);
    final_map = remap(init_map);
    organise_registers_and_maps();
  }
  bool modified = circuit_modified();
  return {circ_, modified};
}

// Tidying up of qregisters and initial and final maps after SWAP adding.
void Routing::organise_registers_and_maps() {
  // Given all the new empty wires with no home, if a qubit isnt in the initial
  // map, find it an unassigned node and shove it there.
  auto all_nodes = original_arc_.get_all_nodes_vec();
  unsigned next_ind = 0;
  Node next_node = all_nodes[next_ind];

  for (const Qubit& qb : circ_.all_qubits()) {
    if (init_map.left.find(qb) == init_map.left.end()) {
      // find next free node
      while (init_map.right.count(next_node)) {
        next_node = all_nodes[++next_ind];
        if (next_ind == all_nodes.size()) {
          throw ArchitectureMismatch(circ_.n_qubits(), current_arc_.n_nodes());
        }
      }
      init_map.left.insert({qb, next_node});
      final_map.left.insert({qb, next_node});
    }
  }

  // Due to the addition of SWAP gates, a qubit path may change, and so it's
  // output boundary ordering may not match the input boundary ordering. The
  // following updates the output boundary to match the ordering of the final
  // slice frontier Make the input boundary match up to node numbering of
  // architecture.
  boundary_t new_boundary;
  qubit_mapping_t reorder_map = bimap_to_map(init_map.left);
  for (const std::pair<const Qubit, Node>& map : reorder_map) {
    Qubit target = final_map.right.at(map.second);
    new_boundary.insert(
        {map.second, circ_.get_in(map.first), circ_.get_out(target)});
    // Which makes it all nicer
  }
  // add classical bits to new_boundary
  for (auto [it, end] =
           circ_.boundary.get<TagType>().equal_range(UnitType::Bit);
       it != end; it++) {
    new_boundary.insert(*it);
  }

  circ_.boundary = new_boundary;
  circ_.update_initial_map(reorder_map);
  circ_.update_final_map(bimap_to_map(final_map.left));
}

// Remap completes the routing algorithm
// slices passed as copy as 3 pass placement needs original preserved
qubit_bimap_t Routing::remap(const qubit_bimap_t& init) {
  qmap = init;

  advance_frontier();
  // The routing algorithm:
  // 1) Slices of circuit are parallelised/packed/whatever into 'timesteps'
  // 2) Swaps are 'proposed' on edges connected to any nodes housing an
  // 'interacting' qubit (interacting -> qubit is in some two qubit interaction
  // in timestep 0) 3) A distance heuristic is used to determine whether the
  // swap proposed will bring interacting qubits closer 4) If a swap is bring
  // interacting qubits together it is compared to a held 'best swap'. The
  // comparison is achieved by applying the same distance heuristic over future
  // timesteps, until one is deemed strictly better. 5) If a succesful swap is
  // found (from 3)), the swap gate is added to the circuit, information on
  // which nodes home which qubits is updated and 1)->4) is repeated. 6) If no
  // succesful swap is found,  Dijkstra's algorithm is used to find a path in
  // the graph between two interacting qubits, which the qubits are then swapped
  // along.
  // ... The pair of interacting qubits in the first timestep with greatest path
  // distance between them is chosen. Algorithm then repeats 1)->4).
  // for(unsigned count=0;slice_frontier_.slice.size()!=0 && count<2;count++){
  while (!slice_frontier_.slice->empty()) {
    SwapResults single_swap = try_all_swaps(current_arc_.get_all_edges_vec());
    if (single_swap.success) {
      route_stats.n_try_all_swaps++;
      perform_action(single_swap.swap);
    } else {
      route_stats.n_solve_furthest++;
      if (!solve_furthest()) {
        throw RoutingFailure();
      }
    }
    advance_frontier();
  }

  qubit_bimap_t final_qmap;
  for (l_const_iterator_t it = qmap.left.begin(); it != qmap.left.end(); ++it) {
    Edge e = slice_frontier_.quantum_out_edges->get<TagKey>()
                 .find(it->first)
                 ->second;
    Vertex v = circ_.target(e);
    while (!circ_.detect_final_Op(v)) {
      e = circ_.get_next_edge(v, e);
      v = circ_.target(e);
    }
    Qubit out_q(circ_.get_id_from_out(v));
    final_qmap.insert({out_q, it->second});
  }

  return final_qmap;
}

}  // namespace tket
