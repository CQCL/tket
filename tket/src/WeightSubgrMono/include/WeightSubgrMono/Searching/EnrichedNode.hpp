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
#include "NodeWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Contains extra data necessary for searching,
 * which nodes themselves don't need to know about.
 */
struct EnrichedNode {
  NodeWSM node;

  /** When we initialise another node from this one by assigning pv->tv,
   * set this to false if it becomes invalid after erasing the assignment.
   * Thus, when backtracking, we can pass over this node quickly.
   */
  bool superficially_valid;

  std::set<VertexWSM> pvs_adjacent_to_newly_assigned_vertices;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
