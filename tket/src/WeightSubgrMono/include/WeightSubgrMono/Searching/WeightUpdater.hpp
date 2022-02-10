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
#include <map>
#include <optional>
#include <set>
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class SearchNodeWrapper;
struct FixedData;

/** Some new pv->tv assignments have been made in a search node.
 * Updates the total weight (scalar product) of this node,
 * whenever new edges become assigned.
 */
class WeightUpdater {
 public:
  /** Returns false if too heavy, and stops updating.
   * @param fixed_data Contains data (including neighbours data) necessary to
   * compute edges and weight.
   * @param assignments All current assignments, including previous nodes.
   * @param node A wrapper around the current node; only the total weight will
   * be changed.
   * @param number_of_assignments_previously_processed_in_this_node We process
   * the assignments in order; all new assignments since the previous call for
   * this node are processed.
   * @param max_weight The maximum weight (scalar product) which we allow;
   * ifthis is exceeded, the calculation breaks off early.
   * @return False if the current weight exceeds the maximum weight (so that
   * we're at a nogood; treat this the same as a graph-theoretic nogood, i.e. no
   * complete monomorphism is possible from extending the current state).
   */
  bool operator()(
      const FixedData& fixed_data, const Assignments& assignments,
      SearchNodeWrapper& node,
      std::size_t number_of_assignments_previously_processed_in_this_node,
      WeightWSM max_weight) const;

 private:
  // Necessary to avoid double counting edges.
  // TODO: but, is there a trick to avoid this?
  mutable std::set<std::pair<VertexWSM, VertexWSM>> m_p_edges_processed;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
