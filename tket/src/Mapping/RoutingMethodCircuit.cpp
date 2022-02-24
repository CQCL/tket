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
    const std::function<std::tuple<Circuit, unit_map_t, unit_map_t>(
        const Circuit&, const ArchitecturePtr&)>
        _route_subcircuit,
    const std::function<bool(const Circuit&, const ArchitecturePtr&)>
        _check_subcircuit,
    unsigned _max_size, unsigned _max_depth)
    : route_subcircuit_(_route_subcircuit),
      check_subcircuit_(_check_subcircuit),
      max_size_(_max_size),
      max_depth_(_max_depth){};

bool RoutingMethodCircuit::check_method(
    const std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  // Get circuit, pass to held check method
  Subcircuit frontier_subcircuit = mapping_frontier->get_frontier_subcircuit(
      this->max_depth_, this->max_size_);
  Circuit frontier_circuit =
      mapping_frontier->circuit_.subcircuit(frontier_subcircuit);
  frontier_circuit.rename_units(
      mapping_frontier->get_default_to_linear_boundary_unit_map());

  return this->check_subcircuit_(frontier_circuit, architecture);
}

unit_map_t RoutingMethodCircuit::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  // Produce subcircuit and circuit
  Subcircuit frontier_subcircuit = mapping_frontier->get_frontier_subcircuit(
      this->max_depth_, this->max_size_);
  Circuit frontier_circuit =
      mapping_frontier->circuit_.subcircuit(frontier_subcircuit);
  frontier_circuit.rename_units(
      mapping_frontier->get_default_to_linear_boundary_unit_map());

  // get routed subcircuit
  std::tuple<Circuit, unit_map_t, unit_map_t> routed_subcircuit =
      this->route_subcircuit_(frontier_circuit, architecture);
  unit_map_t new_labelling = std::get<1>(routed_subcircuit);

  // update unit id at boundary in case of relabelling
  mapping_frontier->update_linear_boundary_uids(new_labelling);

  unit_map_t swap_permutation;
  for (const auto& pair : new_labelling) {
    if (pair.first != pair.second &&
        architecture->node_exists(Node(pair.first))) {
      swap_permutation.insert(pair);
    }
  }
  // permute edges held by unitid at out boundary due to swaps
  mapping_frontier->permute_subcircuit_q_out_hole(
      std::get<2>(routed_subcircuit), frontier_subcircuit);

  // substitute old boundary with new cirucit
  std::get<0>(routed_subcircuit).flatten_registers();
  mapping_frontier->circuit_.substitute(
      std::get<0>(routed_subcircuit), frontier_subcircuit);
  // return initial unit_map_t incase swap network required
  return swap_permutation;
}

}  // namespace tket