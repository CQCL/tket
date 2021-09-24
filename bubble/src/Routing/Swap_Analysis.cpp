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

#include "Architecture/Architectures.hpp"
#include "Circuit/CircPool.hpp"
#include "Routing/Routing.hpp"

namespace tket {

/* Routing Class Methods for picking optimal swaps */

/* Overloaded methods for generating distance vectors */
// Distance vectors comprise of information pertaining to the architectural
// distance between qubits immediately interacting

// Generates distance vector from input interaction vector
std::vector<std::size_t> Routing::generate_distance_vector(
    const Interactions &inter) const {
  const unsigned n = current_arc_.get_diameter();
  // const unsigned n = active_distance_matrix.maxCoeff();
  if (n < 1) {
    throw ArchitectureInvalidity("Architecture has diameter 0.");
  }
  std::vector<std::size_t> dv(n - 1);
  for (auto [n1, n2] : inter) {
    unsigned dist = current_arc_.get_distance(n1, n2);
    if (dist > 1) {
      ++dv[n - dist];
    }
  }
  return dv;
}

// Returns the distance between n1 and p1 and the distance between n2 and p2,
// distance ordered (greatest first)
const std::pair<std::size_t, std::size_t> Routing::pair_dists(
    const Node &n1, const Node &p1, const Node &n2, const Node &p2) const {
  std::size_t curr_dist1 = current_arc_.get_distance(n1, p1);
  std::size_t curr_dist2 = current_arc_.get_distance(n2, p2);
  return (curr_dist1 > curr_dist2) ? std::make_pair(curr_dist1, curr_dist2)
                                   : std::make_pair(curr_dist2, curr_dist1);
}

// Determines if a proposed swap brings interacting qubits closer, improving
// board state.
bool Routing::swap_decreases(
    const Swap &nodes, const Interactions &inte) const {
  Node node1 = nodes.first;
  Node pair1 = inte.at(node1);
  Node node2 = nodes.second;
  Node pair2 = inte.at(node2);

  if (pair1 == node2 || (node1 == pair1 && node2 == pair2)) {
    return false;
  }
  const std::pair<std::size_t, std::size_t> &curr_dists =
      pair_dists(node1, pair1, node2, pair2);
  const std::pair<std::size_t, std::size_t> &news_dists =
      pair_dists(node2, pair1, node1, pair2);

  return news_dists < curr_dists;
}

// Given swap and distance vector, updates distance vector to reflect increment
// change due to swaps nodes
void Routing::increment_distance(
    graph::dist_vec &new_dist_vector, const Swap &pair, int increment) const {
  const unsigned n = current_arc_.get_diameter();
  const unsigned dis_index =
      n - current_arc_.get_distance(pair.first, pair.second);
  if (dis_index < new_dist_vector.size()) {
    new_dist_vector[dis_index] += increment;
  }
}

/* Overloaded method for updating temporary distance vectors due to proposed
 * swaps */
// Updates distance vector from proposed swap using global first slice
// interaction vector solve furthest only at this point ...

// Updates distance vector from presented interaction vector
graph::dist_vec Routing::update_distance_vector(
    const Swap &nodes, std::vector<std::size_t> new_dist_vector,
    const Interactions &inte) const {
  increment_distance(new_dist_vector, {nodes.first, inte.at(nodes.first)}, -2);
  increment_distance(
      new_dist_vector, {nodes.second, inte.at(nodes.second)}, -2);
  increment_distance(new_dist_vector, {nodes.second, inte.at(nodes.first)}, 2);
  increment_distance(new_dist_vector, {nodes.first, inte.at(nodes.second)}, 2);
  return new_dist_vector;
}

// Updates qmap to reflect performed swap
void Routing::update_qmap(qubit_bimap_t &map, const Swap &swap) {
  const Qubit qb1 = map.right.at(swap.first);
  const Qubit qb2 = map.right.at(swap.second);
  map.right.erase(swap.first);
  map.right.erase(swap.second);
  map.left.insert({qb1, swap.second});
  map.left.insert({qb2, swap.first});
}

std::vector<Swap> Routing::candidate_swaps(
    const std::vector<Architecture::Connection> &trial_edges,
    const Interactions &inte) const {
  std::vector<Swap> potential_swaps;
  for (auto [node, adjacent_node] : trial_edges) {
    if (inte.at(node) != node || inte.at(adjacent_node) != adjacent_node) {
      Swap proposed = {node, adjacent_node};
      if (swap_decreases(proposed, inte)) {
        potential_swaps.push_back(proposed);
      }
    }
  }
  return potential_swaps;
}

// Move heuristic in try_all_swaps loop outside, for testing help and easy
// changing?
std::vector<Swap> Routing::cowtan_et_al_heuristic(
    std::vector<Swap> &candidate_swaps,
    const std::vector<std::size_t> &base_dists,
    const Interactions &interac) const {
  const Swap winner = candidate_swaps.back();
  candidate_swaps.pop_back();
  std::vector<std::size_t> winner_distances =
      update_distance_vector(winner, base_dists, interac);
  std::vector<Swap> smaller_set;
  smaller_set.push_back(winner);
  for (const Swap &proposed_swap : candidate_swaps) {
    const std::vector<std::size_t> proposed_distances =
        update_distance_vector(proposed_swap, base_dists, interac);
    const int comp =
        tri_lexicographical_comparison(proposed_distances, winner_distances);
    if (comp == -1) {
      smaller_set.push_back(proposed_swap);
    } else if (comp == 1) {
      smaller_set = {proposed_swap};
      winner_distances = proposed_distances;
    }
  }
  return smaller_set;
}

SwapResults Routing::try_all_swaps(const std::vector<Architecture::Connection>
                                       &trial_edges) {  // don't need to change
  std::vector<Swap> potential_swaps = candidate_swaps(trial_edges, interaction);

  if (potential_swaps.empty()) return {false, {Node(0), Node(0)}};

  RoutingFrontier high_sf = slice_frontier_;

  for (unsigned i = 0; i < config_.depth_limit && !high_sf.slice->empty() &&
                       potential_swaps.size() > 1;
       i++) {
    Interactions interac =
        (i == 0) ? interaction : generate_interaction_frontier(high_sf);
    std::vector<std::size_t> base_dists =
        (i == 0) ? dist_vector : generate_distance_vector(interac);

    potential_swaps =
        cowtan_et_al_heuristic(potential_swaps, base_dists, interac);

    high_sf.next_slicefrontier();
  }

  return {1, potential_swaps.back()};
}

std::vector<Swap> Routing::path_to_swaps(const std::vector<Node> &path) {
  const unsigned len = path.size();
  std::vector<Swap> output_swaps;
  if (len > 2) {
    unsigned halfway = len / 2;
    for (unsigned i = 0; (i < halfway) || ((halfway + 2 + i) < len); i++) {
      if (i < halfway) {
        Swap sw1 = {path[i], path[i + 1]};
        output_swaps.push_back(sw1);
      }
      if ((halfway + 2 + i) < len) {
        Swap sw2 = {path[len - i - 2], path[len - i - 1]};
        output_swaps.push_back(sw2);
      }
    }
  }
  return output_swaps;
}

// If heuristic can't settle on a suitable single swap or pair of swaps, find a
// path between the two interacting qubits at greatest distance and swap along
// it.
bool Routing::solve_furthest() {
  bool success = false;
  std::optional<Node> max_node = std::nullopt;
  unsigned max_dist = 0;
  for (auto [q1, q2] : interaction) {
    unsigned dist = current_arc_.get_distance(q1, q2);
    if (dist > max_dist) {
      max_dist = dist;
      max_node = q1;
    }
  }
  if (!max_node.has_value()) {
    throw ArchitectureInvalidity("Architecture is disconnected");
  }
  Node root = *max_node;
  if (max_dist > 1) {
    Node target = interaction.at(root);
    const std::vector<Node> path = current_arc_.get_path(root, target);
    const std::vector<Swap> swaps_to_perform = path_to_swaps(path);
    for (const Swap &swap : swaps_to_perform) {
      success = true;
      add_swap(swap);
    }
  }
  return success;
}

void Routing::update_central_nodes(
    const Swap &nodes, const Interactions &interac,
    distributed_cx_info &candidate_distributed_cx) {
  if (candidate_distributed_cx.first.first) {
    // TODO: check that there isnt a better way than get_path to do this
    std::vector<Node> path =
        current_arc_.get_path(nodes.first, interac.at(nodes.first));
    candidate_distributed_cx.first.second = path[1];
    if (interac.at(path[1]) != path[1]) {
      candidate_distributed_cx.first.first = false;
    }
  }
  if (candidate_distributed_cx.second.first) {
    // TODO: this uses solve furthest -> maybe a better way just
    // using the distance matrix alone for a speed up?
    std::vector<Node> path =
        current_arc_.get_path(nodes.second, interac.at(nodes.second));
    candidate_distributed_cx.second.second = path[1];
    if (interac.at(path[1]) != path[1]) {
      candidate_distributed_cx.second.first = false;
    }
  }
}

// Difference in distance between interacting qubits between nodes in SWAP can
// only differ by 1 Compares the difference in distance between all given
// interactions, scales them dependent on how timesteps to interaction, and
// returns whether a distributed cx is desired.
void Routing::compare_distributed_cx_distances(
    distributed_cx_info &candidate_distributed_cx,
    const std::pair<std::vector<Node>, std::vector<Node>> &inter_node) {
  std::pair<int, int> distance_check = {0, 0};
  for (unsigned i = 1; i < inter_node.first.size(); i++) {
    distance_check.first += pow(i, config_.distrib_exponent) *
                            (int(current_arc_.get_distance(
                                 inter_node.second[0], inter_node.first[i])) -
                             int(current_arc_.get_distance(
                                 inter_node.first[0], inter_node.first[i])));
  }
  for (unsigned i = 1; i < inter_node.second.size(); i++) {
    distance_check.second += pow(i, config_.distrib_exponent) *
                             (int(current_arc_.get_distance(
                                  inter_node.first[0], inter_node.second[i])) -
                              int(current_arc_.get_distance(
                                  inter_node.second[0], inter_node.second[i])));
  }
  if (distance_check.first < 0) {
    candidate_distributed_cx.first.first = false;
  }
  if (distance_check.second < 0) {
    candidate_distributed_cx.second.first = false;
  }
}

bool check_vertex_is_CX(const Circuit &circ_, const Vertex &v) {
  OpType ot = circ_.get_OpType_from_Vertex(v);
  if (ot != OpType::CX) {
    if (ot == OpType::Conditional) {
      const Conditional &b =
          static_cast<const Conditional &>(*circ_.get_Op_ptr_from_Vertex(v));
      if (b.get_op()->get_type() != OpType::CX) {
        return false;
      }
    } else {
      return false;
    }
  }
  return true;
}
// Method is supplied with a pair of nods with the intention of being swapped.
// Before this SWAP gate is added, this method considers whether a distributed
// CX gate between interacting qubits distance 2 away from eachother is a better
// option The returned bool pair instructs perform_action whether to add a
// distributed CX gate between nodes.first and its partner node and nodes.second
// and its partner node respectively
distributed_cx_info Routing::check_distributed_cx(const Swap &nodes) {
  // 1) Determine which nodes in SWAP gate could complete their CX with a
  // distributed CX gate instead
  distributed_cx_info candidate_distributed_cx = {
      {current_arc_.get_distance(nodes.first, interaction[nodes.first]) == 2,
       Node(0)},
      {current_arc_.get_distance(nodes.second, interaction[nodes.second]) == 2,
       Node(0)}};
  // 1 pt2) Is the vertex a CX gate or Conditioned CX gate?
  auto cx_check = [&](bool candidate, const Qubit &qb) {
    if (candidate)
      return check_vertex_is_CX(
          circ_,
          circ_.target(slice_frontier_.quantum_in_edges->find(qb)->second));
    return true;
  };
  if (!cx_check(
          candidate_distributed_cx.first.first, qmap.right.at(nodes.first)))
    return {{false, Node(0)}, {false, Node(0)}};
  if (!cx_check(
          candidate_distributed_cx.second.first, qmap.right.at(nodes.second)))
    return {{false, Node(0)}, {false, Node(0)}};

  if (candidate_distributed_cx.first.first ||
      candidate_distributed_cx.second.first) {
    // 2) Find number of next interactions for node in SWAP equivalent to
    // config, or reached within depth limit
    std::pair<std::vector<Node>, std::vector<Node>> inter_node = {
        {nodes.first}, {nodes.second}};
    std::pair<unsigned, unsigned> ni_limit = {0, 0};

    RoutingFrontier high_sf = slice_frontier_;

    for (unsigned i = 0; i < config_.distrib_limit && !high_sf.slice->empty() &&
                         (ni_limit.first < config_.interactions_limit ||
                          ni_limit.second < config_.interactions_limit);
         i++) {
      // Find interaction frontier for current slice and find the interacting
      // pairs of nodes for incident SWAP gate.
      Interactions interac =
          (i == 0) ? interaction : generate_interaction_frontier(high_sf);

      if (nodes.first != interac[nodes.first] &&
          ni_limit.first < config_.interactions_limit) {
        inter_node.first.push_back(interac[nodes.first]);
        ni_limit.first++;
      }
      if (nodes.second != interac[nodes.second] &&
          ni_limit.second < config_.interactions_limit) {
        inter_node.second.push_back(interac[nodes.second]);
        ni_limit.second++;
      }
      high_sf.next_slicefrontier();
    }
    if (ni_limit.first > 0 && ni_limit.second > 0) {
      // 3) Compare difference in distances between interacting qubits given the
      // permutation of qubits from added SWAP gate, or not.
      compare_distributed_cx_distances(candidate_distributed_cx, inter_node);
      if (candidate_distributed_cx.first.first ||
          candidate_distributed_cx.second.first) {
        // 4) If desirable, find the central node of the bridge.
        update_central_nodes(nodes, interaction, candidate_distributed_cx);
        return candidate_distributed_cx;
      }
    }
  }
  return {{false, Node(0)}, {false, Node(0)}};
}

// Give a node with a control qubit on it, finds its respective target node and
// node between them, and replaces the CX gate between the control and target
// with a distributed CX
void Routing::add_distributed_cx(
    const Node &cx_node_0, const Node &cx_node_1, const Node &central_node) {
  // Find interacting node for starting_node, find node between them. Also swap
  // control and target node if necessary.

  if (current_arc_.get_distance(cx_node_0, cx_node_1) != 2) {
    throw BridgeInvalid("Bridge Nodes are not distance 2 apart.");
  }
  if (current_arc_.get_distance(cx_node_0, central_node) != 1 ||
      current_arc_.get_distance(cx_node_1, central_node) != 1) {
    throw BridgeInvalid(
        "Central BRIDGE node not adjacent to Control and Target "
        "nodes.");
  }

  route_stats.bridge_count++;
  Edge edge_0 =
      slice_frontier_.quantum_in_edges->find(qmap.right.at(cx_node_0))->second;
  Edge edge_1 =
      slice_frontier_.quantum_in_edges->find(qmap.right.at(cx_node_1))->second;

  // Assign control and target nodes from cx_node_0 and cx_node_1
  // Depends on the port ordering of the cx_node_0 and cx_node_1 corresponding
  // edges attached to the CX vertex
  Node control_node, target_node;
  if (circ_.get_ports(edge_1).second < circ_.get_ports(edge_0).second) {
    control_node = cx_node_1;
    target_node = cx_node_0;
  } else {
    control_node = cx_node_0;
    target_node = cx_node_1;
  }

  // Find qubits associated to each node
  const Qubit control_qb = qmap.right.at(control_node);
  const Qubit central_qb = qmap.right.at(central_node);
  const Qubit target_qb = qmap.right.at(target_node);

  // Initialize variables appropriate for substituting Conditionals with CX
  // gates to Conditionals with BRIDGE gates.
  Op_ptr new_bridge_ptr;
  EdgeVec b_in_edges = {};
  OpType gate_op;
  std::vector<std::tuple<Vertex, port_t, port_t>> classical_edge_info = {};

  Vertex to_be_replaced = slice_frontier_.circ.target(
      slice_frontier_.quantum_in_edges->find(control_qb)->second);
  // If OpType is a Conditional{CX}, replace with Conditional{BRIDGE} instead
  if (circ_.get_OpType_from_Vertex(to_be_replaced) == OpType::Conditional) {
    Op_ptr pt = circ_.get_Op_ptr_from_Vertex(to_be_replaced);
    const Conditional &b = static_cast<const Conditional &>(
        *circ_.get_Op_ptr_from_Vertex(to_be_replaced));
    gate_op = b.get_op()->get_type();
    new_bridge_ptr = std::make_shared<Conditional>(
        get_op_ptr(OpType::BRIDGE, std::vector<Expr>(), 3), b.get_width(),
        b.get_value());
    // Also collect any classical in edges
    b_in_edges = circ_.get_in_edges_of_type(to_be_replaced, EdgeType::Boolean);
    for (Edge e : b_in_edges) {
      classical_edge_info.push_back(
          {circ_.source(e), circ_.get_source_port(e),
           circ_.get_target_port(e)});
    }
  } else {  // else make normal bridge
    new_bridge_ptr = get_op_ptr(OpType::BRIDGE);
    gate_op = circ_.get_OpType_from_Vertex(to_be_replaced);
  }

  if (gate_op != OpType::CX) {
    throw BridgeInvalid(
        "OpType::BRIDGE being substituted for a vertex that isn't "
        "OpType::CX. Please rebase two-qubit primitive to CX gate.");
  }
  // Collect all required Quantum edge information

  Edge control_in_edge =
      slice_frontier_.quantum_in_edges->find(control_qb)->second;
  Edge control_out_edge =
      slice_frontier_.quantum_out_edges->find(control_qb)->second;
  Edge central_edge =
      slice_frontier_.quantum_in_edges->find(central_qb)->second;
  Edge target_in_edge =
      slice_frontier_.quantum_in_edges->find(target_qb)->second;
  Edge target_out_edge =
      slice_frontier_.quantum_out_edges->find(target_qb)->second;

  VertPort control_pred = {
      circ_.source(control_in_edge), circ_.get_source_port(control_in_edge)};
  VertPort central_pred = {
      circ_.source(central_edge), circ_.get_source_port(central_edge)};
  VertPort target_pred = {
      circ_.source(target_in_edge), circ_.get_source_port(target_in_edge)};

  VertPort control_succ = {
      circ_.target(control_out_edge), circ_.get_target_port(control_out_edge)};
  VertPort central_succ = {
      circ_.target(central_edge), circ_.get_target_port(central_edge)};
  VertPort target_succ = {
      circ_.target(target_out_edge), circ_.get_target_port(target_out_edge)};

  //  remove old vertex, add new vertex
  circ_.remove_vertex(
      to_be_replaced, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  Vertex bridge_vert = circ_.add_vertex(new_bridge_ptr);
  // add Boolean edges
  for (std::tuple<Vertex, port_t, port_t> vpp : classical_edge_info) {
    circ_.add_edge(
        {std::get<0>(vpp), std::get<1>(vpp)}, {bridge_vert, std::get<2>(vpp)},
        EdgeType::Boolean);
  }

  unsigned num_classicals = classical_edge_info.size();
  // add control qubit in edge
  Edge control_in = circ_.add_edge(
      control_pred, {bridge_vert, num_classicals}, EdgeType::Quantum);
  // add control qubit out edge
  Edge control_out = circ_.add_edge(
      {bridge_vert, num_classicals}, control_succ, EdgeType::Quantum);
  // add central qubit in edge
  Edge central_in = circ_.add_edge(
      central_pred, {bridge_vert, num_classicals + 1}, EdgeType::Quantum);
  // add central qubit out edge
  Edge central_out = circ_.add_edge(
      {bridge_vert, num_classicals + 1}, central_succ, EdgeType::Quantum);
  // add target qubit in edge
  Edge target_in = circ_.add_edge(
      target_pred, {bridge_vert, num_classicals + 2}, EdgeType::Quantum);
  // add target qubit out edge
  Edge target_out = circ_.add_edge(
      {bridge_vert, num_classicals + 2}, target_succ, EdgeType::Quantum);

  // Remove central_edge which is now going through BRIDGE vertex
  circ_.remove_edge(central_edge);

  unit_frontier_t::iterator control_qb_in_it =
      slice_frontier_.quantum_in_edges->find(control_qb);
  unit_frontier_t::iterator central_qb_in_it =
      slice_frontier_.quantum_in_edges->find(central_qb);
  unit_frontier_t::iterator target_qb_in_it =
      slice_frontier_.quantum_in_edges->find(target_qb);

  slice_frontier_.quantum_in_edges->replace(
      control_qb_in_it, {control_qb, control_in});
  slice_frontier_.quantum_in_edges->replace(
      central_qb_in_it, {central_qb, central_in});
  slice_frontier_.quantum_in_edges->replace(
      target_qb_in_it, {target_qb, target_in});

  // Update slice frontier out edges
  unit_frontier_t::iterator control_qb_out_it =
      slice_frontier_.quantum_out_edges->find(control_qb);
  unit_frontier_t::iterator central_qb_out_it =
      slice_frontier_.quantum_out_edges->find(central_qb);
  unit_frontier_t::iterator target_qb_out_it =
      slice_frontier_.quantum_out_edges->find(target_qb);

  slice_frontier_.quantum_out_edges->replace(
      control_qb_out_it, {control_qb, control_out});
  slice_frontier_.quantum_out_edges->replace(
      central_qb_out_it, {central_qb, central_out});
  slice_frontier_.quantum_out_edges->replace(
      target_qb_out_it, {target_qb, target_out});

  // Remove CX vertex from Slice (i.e. VertexVec) in slice_frontier-
  slice_frontier_.slice->erase(
      std::remove(
          slice_frontier_.slice->begin(), slice_frontier_.slice->end(),
          to_be_replaced),
      slice_frontier_.slice->end());
  slice_frontier_.slice->push_back(bridge_vert);
}

// Suitable swap found, amend all global constructs
void Routing::add_swap(const Swap &nodes) {
  route_stats.swap_count++;
  const Qubit qb1 = qmap.right.at(nodes.first);
  const Qubit qb2 = qmap.right.at(nodes.second);

  update_qmap(qmap, nodes);

  // ---   --X--\ /--
  //     =   |   X
  // ---   --X--/ \--
  // So we insert a SWAP gate and perform the wire swap by changing the output
  // ports

  // find edges using qubits
  EdgeVec preds = {
      slice_frontier_.quantum_in_edges->find(qb1)->second,
      slice_frontier_.quantum_in_edges->find(qb2)->second};

  Vertex swap_v = circ_.add_vertex(OpType::SWAP);
  circ_.rewire(swap_v, preds, {EdgeType::Quantum, EdgeType::Quantum});
  EdgeVec swap_outs = circ_.get_all_out_edges(swap_v);

  circ_.dag[swap_outs[0]].ports.first = 1;
  circ_.dag[swap_outs[1]].ports.first = 0;
  unit_frontier_t::iterator qb1_in_it =
      slice_frontier_.quantum_in_edges->find(qb1);
  slice_frontier_.quantum_in_edges->replace(qb1_in_it, {qb1, swap_outs[0]});
  unit_frontier_t::iterator qb2_in_it =
      slice_frontier_.quantum_in_edges->find(qb2);
  slice_frontier_.quantum_in_edges->replace(qb2_in_it, {qb2, swap_outs[1]});
  unit_frontier_t::iterator qb1_out_it =
      slice_frontier_.quantum_out_edges->find(qb1);
  unit_frontier_t::iterator qb2_out_it =
      slice_frontier_.quantum_out_edges->find(qb2);
  if (preds[0] == qb1_out_it->second) {
    slice_frontier_.quantum_out_edges->replace(qb1_out_it, {qb1, swap_outs[0]});
  } else if (preds[1] == qb2_out_it->second) {
    slice_frontier_.quantum_out_edges->replace(qb2_out_it, {qb2, swap_outs[1]});
  }
}

void Routing::perform_action(const Swap &nodes) {
  distributed_cx_info cdcx = check_distributed_cx(nodes);
  if (cdcx.first
          .first) {  // in current heuristic, both nodes in SWAP being distance
                     // two from target, and closer to next interaction if not
                     // permuted is exceptionally rare (never so far...)
    Node temp = nodes.first;
    add_distributed_cx(temp, interaction[nodes.first], cdcx.first.second);
  } else if (cdcx.second.first) {
    Node temp = nodes.second;
    add_distributed_cx(temp, interaction[nodes.second], cdcx.second.second);
  } else {
    add_swap(nodes);
  }
}
}  // namespace tket
