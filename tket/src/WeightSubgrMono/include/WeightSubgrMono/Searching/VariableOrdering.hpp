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
#include <optional>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
class RNG;
namespace WeightedSubgraphMonomorphism {

/** Out of all the unassigned pattern v, which one should we choose next
 * to assign? Should ONLY be called with fully reduced domains, etc.
 */
class VariableOrdering {
 public:
  struct Result {
    // If null, we couldn't find an unassigned variable.
    // So EITHER we've got a full solution, OR the problem is insoluble
    // (there is an empty domain).
    std::optional<VertexWSM> variable_opt;
    bool empty_domain;
  };

  /** We prefer, first, vertices in the candidate set.
   * Only if no unassigned candidates exist do we allow other vertices.
   * Next, we prefer those with smallest domains.
   */
  Result get_variable(
      const PossibleAssignments& possible_assignments,
      const std::set<VertexWSM>& candidate_vertices, RNG& rng);

 private:
  std::vector<VertexWSM> m_pv_list;

  /** Returns false if an empty domain is detected. */
  bool check_candidate(
      VertexWSM pv, std::size_t domain_size, std::size_t& min_domain_size);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
