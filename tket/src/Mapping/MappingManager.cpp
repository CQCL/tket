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

#include "Mapping/MappingManager.hpp"

#include "Architecture/BestTsaWithArch.hpp"

namespace tket {

MappingManager::MappingManager(const ArchitecturePtr& _architecture)
    : architecture_(_architecture) {}

bool MappingManager::route_circuit(
    Circuit& circuit, const std::vector<RoutingMethodPtr>& routing_methods,
    bool label_isolated_qubits) const {
  return this->route_circuit_with_maps(
      circuit, routing_methods, std::make_shared<unit_bimaps_t>(),
      label_isolated_qubits);
}

bool MappingManager::route_circuit_with_maps(
    Circuit& circuit, const std::vector<RoutingMethodPtr>& routing_methods,
    std::shared_ptr<unit_bimaps_t> maps, bool label_isolated_qubits) const {
  if (circuit.n_qubits() > this->architecture_->n_nodes()) {
    std::string error_string =
        "Circuit has" + std::to_string(circuit.n_qubits()) +
        " logical qubits. Architecture has " +
        std::to_string(this->architecture_->n_nodes()) +
        " physical qubits. Circuit to be routed can not have more "
        "qubits than the Architecture.";
    throw MappingManagerError(error_string);
  }

  // mapping_frontier tracks boundary between routed & un-routed in circuit
  // when initialised, boundary is over output edges of input vertices
  MappingFrontier_ptr mapping_frontier;
  if (!maps->initial.empty() && !maps->final.empty()) {
    mapping_frontier = std::make_shared<MappingFrontier>(circuit, maps);
  } else {
    mapping_frontier = std::make_shared<MappingFrontier>(circuit);
  }
  // updates routed/un-routed boundary

  mapping_frontier->advance_frontier_boundary(this->architecture_);

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
  // there may still be some unlabelled qubits
  if (label_isolated_qubits) {
    circuit_modified = true;

    unit_map_t placement;
    qubit_vector_t to_place;
    std::vector<Node> placed;

    // Find which/if any qubits need placing
    for (const Qubit& q : mapping_frontier->circuit_.all_qubits()) {
      Node n(q);
      if (!this->architecture_->node_exists(n)) {
        // Ancilla qubits can be assigned during routing
        // If some qubits are unplaced then its possible the returned circuit
        // has more qubits than the architecture has nodes, which is bad instead
        // at least assign any unlabelled qubits to any ancilla nodes to prevent
        // this
        if (mapping_frontier->ancilla_nodes_.size() > 0) {
          circuit_modified = true;
          Node ancilla = *mapping_frontier->ancilla_nodes_.begin();
          mapping_frontier->merge_ancilla(q, ancilla);
          mapping_frontier->ancilla_nodes_.erase(
              mapping_frontier->ancilla_nodes_.begin());
          placed.push_back(n);
        } else {
          to_place.push_back(n);
        }
      } else {
        placed.push_back(n);
        // if already placed, make sure qubit retains placement
        placement.insert({n, n});
      }
    }
    // avoid doing std::set_difference unless qubits need to be placed
    unsigned n_placed = to_place.size();
    if (n_placed > 0) {
      std::vector<Node> difference,
          architecture_nodes = this->architecture_->get_all_nodes_vec();
      std::set_difference(
          architecture_nodes.begin(), architecture_nodes.end(), placed.begin(),
          placed.end(), std::inserter(difference, difference.begin()));
      // should always be enough remaining qubits to assign unplaced qubits to
      TKET_ASSERT(difference.size() >= n_placed);
      for (unsigned i = 0; i < n_placed; i++) {
        // naively assign each qubit to some free node
        placement.insert({to_place[i], difference[i]});
        mapping_frontier->update_bimaps(
            mapping_frontier->get_qubit_from_circuit_uid(to_place[i]),
            difference[i]);
      }

      mapping_frontier->update_linear_boundary_uids(placement);
      mapping_frontier->circuit_.rename_units(placement);
    }
  }

  return circuit_modified;
}
}  // namespace tket