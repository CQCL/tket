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

#include "DistancesRawReducer.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** "DistancesRawReducer" only knows about single assignments.
 * However we want to check them in a specific order for best performance,
 * so this class keeps track of which reductions have already been done,
 * and what to do next. (Might seem a bit complicated, but actually
 * negligible in comparison to the reduction calculation).
 */
class DistancesReducer {
 public:
  DistancesReducer(
      const NeighboursData& pattern_ndata,
      NearNeighboursData& pattern_near_ndata,
      const NeighboursData& target_ndata,
      NearNeighboursData& target_near_ndata);

  /** Checks if pv->tv appears to be valid at first glance. */
  bool check(
      const std::pair<VertexWSM, VertexWSM>& assignment,
      unsigned distance) const;

  /** When we start reducing a node, the assignments are stored in a vector
   * which the node controls. This class only needs to know the associated
   * sizes, and gets the assignments data from the node.
   */
  void reset(unsigned distance_value);

  struct Result {
    std::optional<std::pair<VertexWSM, VertexWSM>> impossible_assignment;
    bool new_assignments_created;
    bool nogood_found;
  };

  /** Tries to reduce everything, but breaks off early
   * if new assignments arise in the node due to reductions. */
  Result operator()(NodeWSM& node);

 private:
  DistancesRawReducer m_raw_reducer;

  // Element[i] gives data for the assignments still to be checked
  // for distance i+1.
  // It is simply the next index in the assignments list
  // which should be checked.
  std::vector<std::size_t> m_data;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
