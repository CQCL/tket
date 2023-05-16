// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Mapping/MappingManager.hpp"

#include "tket/Architecture/BestTsaWithArch.hpp"

namespace tket {

MappingManager::MappingManager(const ArchitecturePtr& _architecture)
    : architecture_(_architecture) {}

bool MappingManager::route_circuit(
    Circuit& circuit,
    const std::vector<RoutingMethodPtr>& routing_methods) const {
  // we make bimaps if not present
  unit_bimaps_t maps;
  for (const Qubit& qubit : circuit.all_qubits()) {
    maps.initial.left.insert({qubit, qubit});
    maps.final.left.insert({qubit, qubit});
  }
  return this->route_circuit_with_maps(
      circuit, routing_methods, std::make_shared<unit_bimaps_t>(maps));
}

bool MappingManager::route_circuit_with_maps(
    Circuit& circuit, const std::vector<RoutingMethodPtr>& routing_methods,
    std::shared_ptr<unit_bimaps_t> maps) const {
  if (circuit.n_qubits() > this->architecture_->n_nodes()) {
    std::string error_string =
        "Circuit has" + std::to_string(circuit.n_qubits()) +
        " logical qubits. Architecture has " +
        std::to_string(this->architecture_->n_nodes()) +
        " physical qubits. Circuit to be routed can not have more "
        "qubits than the Architecture.";
    throw MappingManagerError(error_string);
  }

  /**
   * Before mapping we find the set of Circuit UnitID that are both
   * not also Architecture UnitiD, but have no Architecture dependent gates
   * and so could be mapped without being assigned.
   * We then assign each of these UnitID to nominally "bad" Architecture UnitID.
   * In this case a "bad" UnitID is one with low out-degree and large
   * total distance to other UnitID.
   */
  // get unplaced UnitID
  std::set<UnitID> unplaced, placed;
  std::set<Node> subarch_nodes = this->architecture_->nodes();
  for (const Qubit& qubit : circuit.all_qubits()) {
    Node n(qubit);
    if (!this->architecture_->node_exists(n)) {
      unplaced.insert(n);
    } else {
      subarch_nodes.erase(n);
    }
  }
  // get qubit with no physical qubit constraint
  Vertex v;
  Edge e;
  std::set<UnitID> to_assign;
  for (const UnitID& unitid : unplaced) {
    v = circuit.get_in(Qubit(unitid));
    e = circuit.get_nth_out_edge(v, 0);
    while (true) {
      if (circuit.n_out_edges_of_type(v, EdgeType::Quantum) > 1 &&
          circuit.get_OpType_from_Vertex(v) != OpType::Barrier)
        break;
      v = circuit.target(e);
      if (circuit.detect_final_Op(v)) break;
      e = circuit.get_next_edge(v, e);
    }
    if (circuit.detect_final_Op(v)) {
      to_assign.insert(unitid);
    }
  }
  // assign these edge cases to "bad" Architecture Node
  bool modify_at_start = false;
  std::set<Node> reassignable_nodes;
  if (!to_assign.empty()) {
    modify_at_start = true;
    // Use a copy as it finds iteratively bad nodes by removing
    std::vector<Node> subarch_nodes_v;
    std::copy(
        subarch_nodes.begin(), subarch_nodes.end(),
        std::back_inserter(subarch_nodes_v));
    // as we produce a subarch, these "bad" nodes may in practice be fine
    // However we allow these single qubit assignments to be overwritten by
    // dynamic assignment
    Architecture subarc = this->architecture_->create_subarch(subarch_nodes_v);
    std::set<Node> nodes = subarc.remove_worst_nodes(to_assign.size());
    TKET_ASSERT(nodes.size() == to_assign.size());
    // note that as we do this pre MappingFrontier we just need to update the
    // Circuit and bimaps
    std::map<Qubit, Node> relabelling_map;
    auto it = nodes.begin();
    auto jt = to_assign.begin();
    while (it != nodes.end()) {
      relabelling_map.insert({Qubit(*jt), Node(*it)});
      reassignable_nodes.insert(Node(*it));
      ++it;
      ++jt;
    }
    circuit.rename_units(relabelling_map);
    // relabel maps
    for (const std::pair<const Qubit, Node>& q_n : relabelling_map) {
      auto q_initial_it = maps->initial.right.find(q_n.first);
      TKET_ASSERT(q_initial_it != maps->initial.right.end());
      UnitID original_first = q_initial_it->second;
      auto q_final_it = maps->final.left.find(original_first);
      TKET_ASSERT(q_final_it != maps->final.left.end());
      maps->initial.right.erase(q_initial_it);
      maps->final.left.erase(q_final_it);
      maps->initial.left.insert({original_first, q_n.second});
      maps->final.left.insert({original_first, q_n.second});
    }
  }
  // mapping_frontier tracks boundary between routed & un-routed in circuit
  // when initialised, boundary is over output edges of input vertices
  MappingFrontier_ptr mapping_frontier;
  if (!maps->initial.empty() && !maps->final.empty()) {
    mapping_frontier = std::make_shared<MappingFrontier>(circuit, maps);
  } else {
    mapping_frontier = std::make_shared<MappingFrontier>(circuit);
  }
  mapping_frontier->reassignable_nodes_ = reassignable_nodes;
  // updates routed/un-routed boundary
  mapping_frontier->advance_frontier_boundary(this->architecture_);

  /**
   * Criteria for Routing being finished.
   * Each linear edge has reached end of Circuit.
   */
  auto check_finish = [&mapping_frontier]() {
    for (const std::pair<UnitID, VertPort>& pair :
         mapping_frontier->linear_boundary->get<TagKey>()) {
      Edge e = mapping_frontier->circuit_.get_nth_out_edge(
          pair.second.first, pair.second.second);
      Vertex v = mapping_frontier->circuit_.target(e);
      OpType ot = mapping_frontier->circuit_.get_OpType_from_Vertex(v);
      if (!is_final_q_type(ot) && ot != OpType::ClOutput) {
        return false;
      }
    }
    return true;
  };

  bool circuit_modified = !check_finish();
  while (!check_finish()) {
    // The order methods are passed in std::vector<RoutingMethod> is
    // the order they are run
    // If a method performs better but only on specific subcircuits,
    // rank it earlier in the passed vector
    bool valid_methods = false;
    for (const auto& rm : routing_methods) {
      // true => can use held routing method
      std::pair<bool, unit_map_t> bool_map =
          rm->routing_method(mapping_frontier, this->architecture_);
      if (bool_map.first) {
        valid_methods = true;
        if (bool_map.second.size() > 0) {
          std::map<Node, Node> node_map;
          for (const auto& x : bool_map.second) {
            node_map.insert({Node(x.first), Node(x.second)});
          }
          for (const std::pair<Node, Node>& swap :
               BestTsaWithArch::get_swaps(*this->architecture_, node_map)) {
            mapping_frontier->add_swap(swap.first, swap.second);
          }
        }
        break;
      }
    }
    if (!valid_methods) {
      throw MappingManagerError(
          "No RoutingMethod suitable to map given subcircuit.");
    }
    // find next routed/unrouted boundary given updates
    mapping_frontier->advance_frontier_boundary(this->architecture_);
  }
  return (circuit_modified || modify_at_start);
}
}  // namespace tket