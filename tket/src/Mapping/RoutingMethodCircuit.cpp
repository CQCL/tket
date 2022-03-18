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

#include "RoutingMethodCircuit.hpp"

namespace tket {

RoutingMethodCircuit::RoutingMethodCircuit(
    const std::function<std::tuple<bool, Circuit, unit_map_t, unit_map_t>(
        const Circuit&, const ArchitecturePtr&)>
        _route_subcircuit,
    unsigned _max_size, unsigned _max_depth)
    : route_subcircuit_(_route_subcircuit),
      max_size_(_max_size),
      max_depth_(_max_depth){};

std::pair<bool, unit_map_t> RoutingMethodCircuit::routing_method(
    MappingFrontier_ptr& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  // Produce subcircuit and circuit
  Subcircuit frontier_subcircuit = mapping_frontier->get_frontier_subcircuit(
      this->max_depth_, this->max_size_);
  Circuit frontier_circuit =
      mapping_frontier->circuit_.subcircuit(frontier_subcircuit);
  frontier_circuit.rename_units(
      mapping_frontier->get_default_to_linear_boundary_unit_map());

  // get routed subcircuit
  std::tuple<bool, Circuit, unit_map_t, unit_map_t> routed_subcircuit =
      this->route_subcircuit_(frontier_circuit, architecture);

  if (!std::get<0>(routed_subcircuit)) {
    return {false, {}};
  }

  // update unit id at boundary in case of relabelling
  // The route_subcircuit_ method populates its initial map
  // with unit ids from the circuit. e.g. Initial map from frontier ==
  // q[0]:unplaced[0], circuit.all_qubits() == unplaced[0]. Then the produced
  // initial map == unplaced[0]:node[0] We have to update the initial map to
  // q[0]:node[0]
  mapping_frontier->update_linear_boundary_uids(std::get<2>(routed_subcircuit));
  for (const auto& pair : std::get<2>(routed_subcircuit)) {
    mapping_frontier->update_bimaps(
        mapping_frontier->get_qubit_from_circuit_uid(pair.first), pair.second);
  }

  unit_map_t swap_permutation;
  for (const auto& pair : std::get<2>(routed_subcircuit)) {
    if (pair.first != pair.second &&
        architecture->node_exists(Node(pair.first))) {
      swap_permutation.insert(pair);
    }
  }
  // permute edges held by unitid at out boundary due to swaps
  mapping_frontier->permute_subcircuit_q_out_hole(
      std::get<3>(routed_subcircuit), frontier_subcircuit);

  // substitute old boundary with new cirucit
  std::get<1>(routed_subcircuit).flatten_registers();
  mapping_frontier->circuit_.substitute(
      std::get<1>(routed_subcircuit), frontier_subcircuit);
  // return initial unit_map_t incase swap network required
  return {true, swap_permutation};
}

}  // namespace tket