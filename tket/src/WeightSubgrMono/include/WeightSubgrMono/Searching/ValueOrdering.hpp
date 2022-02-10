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
#include <set>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct SharedData;

/** Given an unassigned PV, choose a possible target vertex TV from Dom(PV)
 * so that the search can proceed with trying pv->tv.
 * Different heuristics may give better performance.
 * It's interesting that there are multiple competing heuristics, even in the
 * unweighted case, so the idea is to try different ones
 * (or even run several in parallel) to improve performance.
 */
class ValueOrdering {
 public:
  /** By default, just choose a TV with maximum degree.
   * Override with better heuristics for better performance.
   * @param possible_values The domain of the pattern vertex PV.
   * @param shared_data Data about graphs, etc. etc. to assist with the
   * decision.
   * @return The chosen TV from Dom(PV).
   */
  virtual VertexWSM get_target_value(
      const std::set<VertexWSM>& possible_values,
      const SharedData& shared_data);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
