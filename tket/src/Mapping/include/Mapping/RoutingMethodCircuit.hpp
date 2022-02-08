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

#include "Mapping/RoutingMethod.hpp"

namespace tket {

class RoutingMethodCircuit : public RoutingMethod {
 public:
  virtual ~RoutingMethodCircuit() {}
  /**
   * RoutingMethodCircuit objects hold methods for partially routing subcircuits
   * in the incremental routing of full circuits.
   *
   * @param _route_subcircuit Function ptr for partial routing method
   * @param _check_subcircuit Function ptr for confirming if method sufficient
   * @param _max_size Max number of gates in partial routing circuit
   * @param _max_depth Max depth of partial routing circuit
   */
  RoutingMethodCircuit(
      const std::function<std::tuple<Circuit, unit_map_t, unit_map_t>(
          const Circuit&, const ArchitecturePtr&)>
          _route_subcircuit,
      const std::function<bool(const Circuit&, const ArchitecturePtr&)>
          _check_subcircuit,
      unsigned _max_size, unsigned _max_depth);

  /**
   * @param mapping_frontier Contains boundary of gates to be checked for method
   * @param architecture Architecture method would work with if permitted
   * @return true if method can route subcircuit, false if not
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const;

  /**
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const;

 private:
  const std::function<std::tuple<Circuit, unit_map_t, unit_map_t>(
      const Circuit&, const ArchitecturePtr&)>
      route_subcircuit_;
  const std::function<bool(const Circuit&, const ArchitecturePtr&)>
      check_subcircuit_;
  unsigned max_size_, max_depth_;
};

JSON_DECL(RoutingMethod);

}  // namespace tket
