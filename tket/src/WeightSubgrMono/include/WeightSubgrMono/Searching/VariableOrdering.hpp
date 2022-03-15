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
#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct SearchNode;
struct SharedData;

/** Out of all the unassigned pattern v, which one should we choose next
 * to assign? Should ONLY be called with fully reduced domains, etc.
 */
class VariableOrdering {
 public:

  /** 
   * @param node The node containing data on the current domains.
   * @param assignments All current assignments, not just those in the current
   * node.
   * @param shared_data Contains all extra data which could be useful for making
   * the decision.
   * @return The pattern vertex PV which we should assign next (choosing the
   * associated TV is a separate task).
   */
  VertexWSM choose_next_variable(
      const SearchNode& node, const Assignments& assignments,
      SharedData& shared_data);

 private:

  std::vector<VertexWSM> m_pattern_vertices_with_smallest_domain;

  // Fills m_pattern_vertices_with_smallest_domain,
  // such that we are unassigned AND adjacent
  // to an assigned p-vertex if possible.
  void fill_pattern_vertices_with_smallest_domain(
      const SearchNode& node, const Assignments& assignments,
      SharedData& shared_data);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
