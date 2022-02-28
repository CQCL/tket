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

#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"

namespace tket {

class AASLabellingMethod : public RoutingMethod {
 public:
  /**
   * Checking and Routing methods redefined for dynamically assigning qubits to
   * some Architecture.
   */
  AASLabellingMethod(){};

  /**
   * will place all the qubits of the given circuit that are not placed at the
   * moment. All nodes assigend to placed qubits will not be changed
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return bool if the method has been executed and logical to Physical
   * mapping at boundary due to modification.
   *
   */
  std::pair<bool, unit_map_t> routing_method(
      std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const override;

  nlohmann::json serialize() const override;

  static AASLabellingMethod deserialize(const nlohmann::json& j);
};
}  // namespace tket