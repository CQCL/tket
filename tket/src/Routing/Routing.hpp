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

#ifndef _TKET_Routing_H_
#define _TKET_Routing_H_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "Placement.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/Json.hpp"

namespace tket {

// 2 (adjacent) nodes proposed to have their concurrent qubit states swapped
typedef std::pair<Node, Node> Swap;
// node i is interacting with element (j) at i, if i==j not interacting
typedef std::map<Node, Node> Interactions;
typedef std::vector<Vertex> qubit_map_vector_t;
typedef std::pair<std::pair<bool, Node>, std::pair<bool, Node>>
    distributed_cx_info;
// TODO remove
// qubit_map_vector_t map2vec(qubit_bimap_t map, unsigned total);
struct SwapResults {  // results of try_all_swaps algorithm
  bool success;
  Swap swap;
};

/* Error Handling for Routing Circuits */
class ArchitectureMismatch : public std::logic_error {
 public:
  ArchitectureMismatch(unsigned circ_no, unsigned arch_no)
      : std::logic_error(
            std::to_string(circ_no) + " " + std::to_string(arch_no)) {
    tket_log()->error(
        "Incorrect number of nodes in the architecture. "
        "Qubits in circuit: {}, nodes in architecture: {}",
        circ_no, arch_no);
  }
};

class QMapRange : public std::logic_error {
 public:
  explicit QMapRange(const std::string &message) : std::logic_error(message) {}
};

class NodesRange : public std::logic_error {
 public:
  NodesRange(int nodes, int qubit)
      : std::logic_error(std::to_string(nodes) + " " + std::to_string(qubit)) {
    tket_log()->error(
        "Qubit indexing larger than number of available qubits."
        "Available Qubits: {}, Qubit Index: {}",
        nodes, qubit);
  }
};

class ArchitectureFull : public std::logic_error {
 public:
  ArchitectureFull()
      : std::logic_error(
            "No suitable node found in findBestNode => all nodes already "
            "used") {}
};

class NodeAlreadyActive : public std::logic_error {
 public:
  explicit NodeAlreadyActive(int node)
      : std::logic_error(std::to_string(node)) {
    tket_log()->error("Node {} already active.", node);
  }
};

class NodeInactive : public std::logic_error {
 public:
  explicit NodeInactive(int node) : std::logic_error(std::to_string(node)) {
    tket_log()->error("Node {} inactive.", node);
  }
};

class RoutingFailure : public std::logic_error {
 public:
  RoutingFailure()
      : std::logic_error(
            "Routing failed to complete. Note: Check your architecture "
            "is connected.") {}
};

class BridgeInvalid : public std::logic_error {
 public:
  explicit BridgeInvalid(const std::string &message)
      : std::logic_error(message) {}
};

class BridgePathIncorrect : public std::logic_error {
 public:
  explicit BridgePathIncorrect(int path_size)
      : std::logic_error(std::to_string(path_size)) {
    tket_log()->error("Path found has size {} which is invalid.", path_size);
  }
};

// structure of configuration parameters for routing
struct RoutingConfig {
  // circuit look ahead limit for SWAP picking
  unsigned depth_limit;
  // circuit look ahead limit for Distributed CX gate checking
  unsigned distrib_limit;
  // number of interactions considered in Distributed CX gate checking
  unsigned interactions_limit;
  // Whether to use a Distributed CX gate instead of a SWAP and a CX is
  // determined by comparing the distance between some interacting pairs of
  // qubits with and without the permutation. Changing distrib_exponent changes
  // how much later interactions are considered. distrib_exponent < 0 => less
  // effect from later interactoins, distrib_exponent > 0 => greater effect,
  // distrib_exponent = 0 => no effect
  double distrib_exponent;
  // Constructors
  RoutingConfig(
      unsigned _depth_limit, unsigned _distrib_limit,
      unsigned _interactions_limit, const double &_distrib_exponent)
      : depth_limit(_depth_limit),
        distrib_limit(_distrib_limit),
        interactions_limit(_interactions_limit),
        distrib_exponent(_distrib_exponent) {}

  RoutingConfig() : RoutingConfig(50, 75, 10, 0) {}

  bool operator==(const RoutingConfig &other) const;
};

JSON_DECL(RoutingConfig)

// stores and tracks the points of the circuit up to which has been solved
struct RoutingFrontier {
  // set of 2qb vertices which need to be solved for
  std::shared_ptr<Slice> slice;
  // Quantum Edges coming in to vertices in slice, indexed by qubit
  std::shared_ptr<unit_frontier_t> quantum_in_edges;
  // Quantum Edges leaving vertices in slice, indexed by qubit
  std::shared_ptr<unit_frontier_t> quantum_out_edges;
  // Boolean edges coming in to vertices in slice. Guarantees that all edges
  // into every vertex in slice is represented in next_cut
  std::shared_ptr<b_frontier_t> classical_in_edges;

  // reference to circuit that it acts on
  const Circuit &circ;

  explicit RoutingFrontier(const Circuit &_circ);
  // initialise at front of circuit
  void init();
  // move to next slice
  void next_slicefrontier();
};

// remove node from architecture as long as subgraph remains connected. Nodes
// not in map from architecture if possible
void remove_unmapped_nodes(
    Architecture &arc, qubit_bimap_t &map, Circuit &circ);

bool subgraph_remove_if_connected(
    Architecture &arc, const Architecture &subarc, const Node &node);

// remove nodes not in map from architecture if possible
void remove_unmapped_nodes(
    Architecture &arc, qubit_bimap_t &map, Circuit &circ);

Circuit autoroute(const Circuit &circ, const Architecture &arc);

class RoutingTester;
/* Routing class, contains solve method for transforming a circuit such that
all it's multi-qubit interactions are adjacent for some specificed architecture.
*/
class Routing {
 public:
  struct Stats {
    unsigned n_try_all_swaps;
    unsigned n_solve_furthest;
    unsigned swap_count;
    unsigned bridge_count;
    Stats()
        : n_try_all_swaps(0),
          n_solve_furthest(0),
          swap_count(0),
          bridge_count(0) {}
  };

  /* Class Constructor */
  Routing(const Circuit &_circ, const Architecture &_arc);
  /* Solve Method */
  // solve using default mapping (line_placement) and default config
  // Default RoutingConfig provides a set of parameters that use all available
  // features of Routing, but are not specialised for a certain architecture:
  // depth_limit = 50
  // distrib_limit = 75
  // interactions_limit = 10
  // distrib_exponent = 0
  // This configuration is used for any solve method that does not have config
  // specified.

  // solve with default mapping and provided config
  std::pair<Circuit, bool> solve(const RoutingConfig &_config = {});
  qubit_bimap_t remap(const qubit_bimap_t &init);
  void organise_registers_and_maps();

  // TODO:: Make relevant and useful again
  qubit_mapping_t return_final_map() const;
  qubit_mapping_t return_initial_map() const;
  /* Getters*/
  std::vector<Node> get_active_nodes() const;

  RoutingFrontier get_slicefrontier() const { return slice_frontier_; }
  Stats get_stats() const { return route_stats; }

 private:
  // Circuit being solved
  Circuit circ_;
  // RoutingFrontier tracking the position whcih has been solved up to
  RoutingFrontier slice_frontier_;
  // Configuration settings for routing
  RoutingConfig config_;

  // Architecture being solved for and the original architecture given
  Architecture current_arc_;
  Architecture original_arc_;

  // Which qubits are interacting and total distance of a board state for
  // interacting qubits
  Interactions interaction;
  // Total distance of a board state for interacting qubits
  graphs::dist_vec dist_vector;

  Stats route_stats;

  boundary_t original_boundary;

  // Various qubit mappings. Qmap is used as the algorithim proceeds, the
  // initial map is assigned from placement and the final map displays where
  // qubits end up while routed. Relative mapping is what the final mapping
  // would be if initial mapping was sequential.
  qubit_bimap_t qmap, init_map, final_map;

  /* Swap_Analysis.cpp methods */
  // Methods used in determining the best Swap for a given board state and
  // implementing it
  void increment_distance(
      graphs::dist_vec &new_dist_vector, const Swap &pair, int increment) const;
  graphs::dist_vec generate_distance_vector(const Interactions &inter) const;
  graphs::dist_vec update_distance_vector(
      const Swap &nodes, graphs::dist_vec new_dist_vector,
      const Interactions &inte) const;
  const std::pair<std::size_t, std::size_t> pair_dists(
      const Node &n1, const Node &p1, const Node &n2, const Node &p2) const;
  bool swap_decreases(const Swap &nodes, const Interactions &inte) const;
  std::vector<Swap> candidate_swaps(
      const std::vector<Architecture::Connection> &trial_edges,
      const Interactions &inte) const;
  std::vector<Swap> cowtan_et_al_heuristic(
      std::vector<Swap> &candidate_swaps, const graphs::dist_vec &base_dists,
      const Interactions &interac) const;
  SwapResults try_all_swaps(
      const std::vector<Architecture::Connection> &trial_edges);

  static void update_qmap(qubit_bimap_t &map, const Swap &swap);
  void update_central_nodes(
      const Swap &nodes, const Interactions &interac,
      distributed_cx_info &candidate_distributed_cx);

  void compare_distributed_cx_distances(
      distributed_cx_info &candidate_distributed_cx,
      const std::pair<std::vector<Node>, std::vector<Node>> &inter_node);
  distributed_cx_info check_distributed_cx(const Swap &nodes);
  void add_distributed_cx(
      const Node &control_node, const Node &target_node,
      const Node &central_node);
  void add_swap(const Swap &nodes);
  void perform_action(const Swap &nodes);

  // Dijkstras algorithm methods
  static std::vector<Swap> path_to_swaps(const std::vector<Node> &path);

  bool solve_furthest();

  /* Slice_Maniupation.cpp methods */
  // find nodes for qubits, activating if necessary
  std::vector<Node> nodes_from_qubits(const qubit_vector_t &qubs);
  // Advances slice frontier past any two_qubit operations on adjacent nodes
  bool advance_frontier();

  bool circuit_modified() const;

  friend class RoutingTester;
  // generate interaction vectors from slice_frontiers, qubits=true means return
  // qubit interactions rather than node
  Interactions generate_interaction_frontier(
      const RoutingFrontier &slice_front);
  /* Qubit_Placement.cpp methods */
  // Methods for producing a good intial qubit mapping to an architecture from
  // given circuit

  // void print_qubitlines(QubitLineList &in);

  /* Board_Analysis.cpp routing methods */
  Node find_best_inactive_node(
      const Node &target_node, const Architecture &arc) const;
  void activate_node(const Node &node);
  void reactivate_qubit(const Qubit &qb, const Qubit &target);
};

class RoutingTester {
 private:
  Routing *router;

 public:
  explicit RoutingTester(Routing *_router) : router(_router) {}

  Interactions get_interaction(const RoutingFrontier &sf);
  void set_qmap(qubit_bimap_t _qmap);
  void next_sf(RoutingFrontier &sf);
  Circuit *get_circ();
  void set_config(const RoutingConfig &_config);
  // Wrappers of private methods for testing?
  void increment_distance(
      graphs::dist_vec &new_dist_vector, const Swap &pair, int increment) const;
  graphs::dist_vec generate_distance_vector(const Interactions &inter) const;
  graphs::dist_vec update_distance_vector(
      const Swap &nodes, graphs::dist_vec new_dist_vector,
      const Interactions &inte) const;
  const std::pair<unsigned, unsigned> pair_dists(
      const Node &n1, const Node &p1, const Node &n2, const Node &p2) const;
  bool swap_decreases(const Swap &nodes, const Interactions &inte) const;
  std::vector<Swap> candidate_swaps(
      const std::vector<Architecture::Connection> &trial_edges,
      const Interactions &inte) const;
  std::vector<Swap> cowtan_et_al_heuristic(
      std::vector<Swap> &candidate_swaps, const graphs::dist_vec &base_dists,
      const Interactions &interac) const;
  void update_qmap(qubit_bimap_t &map, const Swap &swap);
  std::vector<Swap> path_to_swaps(const std::vector<Node> &path) const;
  qubit_bimap_t set_default_initial_map(
      std::optional<node_vector_t> canonical_node_order = std::nullopt);
  void initialise_slicefrontier();
  void add_distributed_cx(
      const Node &control_node, const Node &target_node,
      const Node &central_node);
  distributed_cx_info check_distributed_cx(const Swap &nodes);
  void advance_frontier();
  void set_interaction();
};

}  // namespace tket

#endif
