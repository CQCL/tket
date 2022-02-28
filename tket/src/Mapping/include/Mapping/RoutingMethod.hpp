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

#pragma once

#include "Mapping/MappingFrontier.hpp"
#include "Utils/Json.hpp"

namespace tket {

class RoutingMethod {
 public:
  RoutingMethod(){};
  virtual ~RoutingMethod() {}

  /**
   * routing_method modifies circuit held in mapping_frontier with gates for the
   * purpose of moving circuit closer to one physically permitted by given
   * architecture. Returns a pair with a bool returning whether any modification
   * was made and a new initial mapping of qubits in case permutation via swap
   * network is then required, or new ancilla qubits are added. This is
   * completed by converting boundary subcircuit in mapping frontier to a
   * Circuit object which is then passed to route_subcircuit_ as defined in the
   * constructor.
   *
   * Overloaded parameter mapping_frontier contains boundary of routed/unrouted
   * circuit for modifying.
   * Overloaded parameter architecture provides physical constraints
   *
   * @return Whether circuit is modified and Logical to Physical mapping at
   * boundary due to modification.
   *
   */
  virtual std::pair<bool, unit_map_t> routing_method(
      MappingFrontier_ptr& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const {
    return {false, {}};
  }

  virtual nlohmann::json serialize() const {
    nlohmann::json j;
    j["name"] = "RoutingMethod";
    return j;
  }
};

typedef std::shared_ptr<const RoutingMethod> RoutingMethodPtr;

}  // namespace tket