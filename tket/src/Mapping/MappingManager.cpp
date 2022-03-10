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
    Circuit& circuit,
    const std::vector<RoutingMethodPtr>& routing_methods) const {
  return this->route_circuit_with_maps(
      circuit, routing_methods, std::make_shared<unit_bimaps_t>());
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
  return circuit_modified;
}
}  // namespace tket