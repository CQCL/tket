// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <optional>

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
class RNG;
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;

/** Out of all the unassigned pattern v, which one should we choose next
 * to assign? Should ONLY be called with fully reduced domains, etc.
 */
class VariableOrdering {
 public:
  /** Our choice of variable PV to assign (if possible). */
  struct Result {
    // If null, we couldn't find an unassigned variable.
    // So EITHER we've got a full solution, OR the problem is insoluble
    // (there is an empty domain).
    std::optional<VertexWSM> variable_opt;
    bool empty_domain;
  };

  /** Choose a variable PV to assign, if possible.
   * @param accessor An object giving us access to domains data of the current
   * node.
   * @param rng A random number generator, in case we need to choose between
   * multiple equally good choices.
   * @return Our choice of PV to assign next
   */
  Result get_variable(DomainsAccessor& accessor, RNG& rng);

 private:
  std::vector<VertexWSM> m_pv_list;
  std::vector<VertexWSM> m_work_vector;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
