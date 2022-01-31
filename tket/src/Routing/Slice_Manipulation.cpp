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

#include <algorithm>
#include <numeric>

#include "Routing.hpp"

namespace tket {

RoutingFrontier::RoutingFrontier(const Circuit& _circ) : circ(_circ) { init(); }
void RoutingFrontier::init() {
  VertexVec input_slice;
  quantum_in_edges = std::make_shared<unit_frontier_t>();
  classical_in_edges = std::make_shared<b_frontier_t>();

  for (const Qubit& qb : circ.all_qubits()) {
    Vertex input = circ.get_in(qb);
    input_slice.push_back(input);
    Edge candidate = circ.get_nth_out_edge(input, 0);
    quantum_in_edges->insert({qb, circ.skip_irrelevant_edges(candidate)});
  }
  for (const Bit& bit : circ.all_bits()) {
    Vertex input = circ.get_in(bit);
    EdgeVec candidates = circ.get_nth_b_out_bundle(input, 0);
    classical_in_edges->insert({bit, candidates});
  }

  CutFrontier next_cut = circ.next_cut(quantum_in_edges, classical_in_edges);
  slice = next_cut.slice;
  quantum_out_edges = next_cut.u_frontier;
}

void RoutingFrontier::next_slicefrontier() {
  quantum_in_edges = std::make_shared<unit_frontier_t>();
  classical_in_edges = std::make_shared<b_frontier_t>();
  for (const std::pair<UnitID, Edge>& pair : quantum_out_edges->get<TagKey>()) {
    Edge new_e = circ.skip_irrelevant_edges(pair.second);
    quantum_in_edges->insert({pair.first, new_e});
    Vertex targ = circ.target(new_e);
    EdgeVec targ_classical_ins =
        circ.get_in_edges_of_type(targ, EdgeType::Boolean);
    classical_in_edges->insert(
        {Bit("frontier_bit", pair.first.index()), targ_classical_ins});
  }

  CutFrontier next_cut = circ.next_cut(quantum_in_edges, classical_in_edges);
  slice = next_cut.slice;
  quantum_out_edges = next_cut.u_frontier;
}

std::vector<Node> Routing::nodes_from_qubits(const qubit_vector_t& qubs) {
  std::vector<Node> nodes;
  unsigned start = 0;
  if (qmap.empty()) {
    Node node0 = *(original_arc_.max_degree_nodes().begin());
    activate_node(node0);
    qmap.left.insert({qubs[0], node0});
    init_map.left.insert({qubs[0], node0});
    nodes.push_back(node0);
    start++;
  }

  for (unsigned i = start; i < qubs.size(); i++) {
    l_const_iterator_t node_find = qmap.left.find(qubs[i]);
    if (node_find == qmap.left.end()) {
      if (i < qubs.size() - 1 &&
          qmap.left.find(qubs[i + 1]) !=
              qmap.left.end()) {  // TODO: Could this if condition cause some
                                  // nasty non determinism?
        reactivate_qubit(qubs[i], qubs[i + 1]);
        nodes.push_back(qmap.left.at(qubs[i]));
      } else {
        if (i != 0) {
          reactivate_qubit(qubs[i], qubs[0]);
          nodes.push_back(qmap.left.at(qubs[i]));
        } else {
          reactivate_qubit(qubs[i], qmap.begin()->left);
          nodes.push_back(qmap.left.at(qubs[i]));
        }
      }
    } else {
      nodes.push_back(node_find->second);
    }
  }
  return nodes;
}

/*
Advances slice frontier past any two_qubit operations on adjacent nodes
*/
bool Routing::advance_frontier() {
  bool found_adjacent_op = true;
  while (found_adjacent_op && !slice_frontier_.slice->empty()) {
    found_adjacent_op = false;
    for (const Vertex& vert : *slice_frontier_.slice) {
      qubit_vector_t qubs;
      for (const Edge& q_out :
           circ_.get_out_edges_of_type(vert, EdgeType::Quantum)) {
        for (const std::pair<UnitID, Edge>& pair :
             slice_frontier_.quantum_out_edges->get<TagKey>()) {
          if (pair.second == q_out) {
            qubs.push_back(Qubit(pair.first));
            break;
          }
        }
      }
      // Find OpType. If OpType is a Conditional, unpack to find vertex inside.
      // If it's nested, this will fail.
      OpType vert_type = circ_.get_OpType_from_Vertex(vert);
      if (vert_type == OpType::Conditional) {
        const Conditional& b = static_cast<const Conditional&>(
            *circ_.get_Op_ptr_from_Vertex(vert));
        vert_type = b.get_op()->get_type();
      }

      // the vertex must be two qubits or a bridge, which we can skip past

      if (qubs.size() != 2 && vert_type != OpType::BRIDGE &&
          vert_type != OpType::Barrier) {
        throw(CircuitInvalidity(
            "Vertex has " + std::to_string(qubs.size()) +
            " qubits, expected 2."));
      }
      // BRIDGE gates are guaranteed to be across 3 adjacent nodes,
      // already mapped so can just be read directly from the qmap
      // otherwise, qubits may need to be activated first
      std::vector<Node> nods = nodes_from_qubits(qubs);

      bool all_qbs_adjacent = true;
      for (unsigned i = 0; i < nods.size() - 1; i++) {
        all_qbs_adjacent &=
            (current_arc_.get_distance(nods[i], nods[i + 1]) == 1);
      }
      if (all_qbs_adjacent ||
          vert_type == OpType::Barrier) {  // if by eachother
        found_adjacent_op = true;          // i.e. at least one 2qb gate has
                                           // been able to run
        // for all qubits skip subsequent single qubit vertices to move
        // in edges to be prior to next multiqubit vertex
        for (const Qubit& qub : qubs) {
          Edge new_e = circ_.skip_irrelevant_edges(
              slice_frontier_.quantum_out_edges->find(qub)->second);
          slice_frontier_.quantum_in_edges->replace(
              slice_frontier_.quantum_in_edges->find(qub), {qub, new_e});
          Vertex targ = circ_.target(new_e);
          EdgeVec targ_classical_ins =
              circ_.get_in_edges_of_type(targ, EdgeType::Boolean);
          Bit b("frontier_bit", qub.index());
          if (slice_frontier_.classical_in_edges->find(b) ==
              slice_frontier_.classical_in_edges->end()) {
            slice_frontier_.classical_in_edges->insert({b, targ_classical_ins});
          } else {
            slice_frontier_.classical_in_edges->replace(
                slice_frontier_.classical_in_edges->find(b),
                {b, targ_classical_ins});
          }
        }
      }
    }
    if (found_adjacent_op) {
      CutFrontier next_cut = circ_.next_cut(
          slice_frontier_.quantum_in_edges, slice_frontier_.classical_in_edges);
      slice_frontier_.slice = next_cut.slice;
      slice_frontier_.quantum_out_edges = next_cut.u_frontier;
      slice_frontier_.classical_in_edges = std::make_shared<b_frontier_t>();
      for (const std::pair<UnitID, Edge>& pair :
           slice_frontier_.quantum_in_edges->get<TagKey>()) {
        Vertex targ = circ_.target(pair.second);
        EdgeVec targ_classical_ins =
            circ_.get_in_edges_of_type(targ, EdgeType::Boolean);
        Bit b("frontier_bit", pair.first.index());
        slice_frontier_.classical_in_edges->insert({b, targ_classical_ins});
      }
    }
  }

  interaction = generate_interaction_frontier(slice_frontier_);  // reset
  dist_vector = generate_distance_vector(interaction);
  return found_adjacent_op;
}

Interactions Routing::generate_interaction_frontier(
    const RoutingFrontier& slice_front) {
  Interactions inter;
  for (const UnitID& uid : current_arc_.nodes()) {
    Node n(uid);
    inter.insert({n, n});
  }
  for (const Vertex& vert : *slice_front.slice) {
    qubit_vector_t qubs;
    for (const Edge& q_out :
         circ_.get_out_edges_of_type(vert, EdgeType::Quantum)) {
      for (const std::pair<UnitID, Edge>& pair :
           slice_front.quantum_out_edges->get<TagKey>()) {
        if (pair.second == q_out) {
          qubs.push_back(Qubit(pair.first));
          break;
        }
      }
    }
    // if generate_interaction_frontier called with slice_frontier_ no ops with
    // more than two qubits will be present if generate_interaction_frontier
    // called with frontier made in try_all_swaps or check_distributed_cx,
    // Barrier Op possible. If barrier op in slice, don't add qubits in barrier
    // interaction to Interactions
    if (qubs.size() != 2) {
      if (circ_.get_OpType_from_Vertex(vert) == OpType::Barrier) continue;
      throw CircuitInvalidity(
          "Vertex has " + std::to_string(qubs.size()) + " qubits, expected 2.");
    }

    l_const_iterator_t node0_find = qmap.left.find(qubs[0]);
    l_const_iterator_t node1_find = qmap.left.find(qubs[1]);
    if (node0_find != qmap.left.end() && node1_find != qmap.left.end()) {
      Node one = node0_find->second;
      Node two = node1_find->second;
      inter[one] = two;
      inter[two] = one;
    }
  }
  return inter;
}

}  // namespace tket
