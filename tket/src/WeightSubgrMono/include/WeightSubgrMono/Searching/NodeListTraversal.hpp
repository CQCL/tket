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
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct NodesRawData;
class NodesRawDataWrapper;

/** This is concerned with moving up and down a NodesRawData object,
 * i.e. used for traversing the complete search tree;
 * "vertical" motion, as explained in NodesRawData.hpp.
 */
class NodeListTraversal {
 public:
  /** The wrapped NodesRawData object will be directly altered.
   * @param raw_data_wrapper A wrapper around the raw data, which will be
   * directly manipulated by this NodeListTraversal object. The raw data is
   * shared with other manipulation objects, such as DomainsAccessor.
   */
  explicit NodeListTraversal(NodesRawDataWrapper& raw_data_wrapper);

  /** Simply return all TVs which occur in some domain, somewhere. */
  std::set<VertexWSM> get_used_target_vertices() const;

  /** Moves up (i.e., decreases current_node_level) UNTIL it reaches
   * an apparently valid node; returns false if it cannot
   * (runs out of valid nodes).
   * @return False if it cannot move up any further to reach a valid node (which
   * means that the search is finished).
   */
  bool move_up();

  /** Make the given assignment PV->TV, which must be valid,
   * and move down to a new node, taking care of new assignments.
   * @param p_vertex A pattern vertex PV, for the assignment PV->TV (which must
   * be valid).
   * @param t_vertex A target vertex tV, for the assignment PV->TV (which must
   * be valid).
   */
  void move_down(VertexWSM p_vertex, VertexWSM t_vertex);

  enum class ImpossibleAssignmentAction {
    IGNORE_CURRENT_NODE,
    PROCESS_CURRENT_NODE
  };

  /** We've newly discovered that PV->TV is always impossible.
   * Erase it from all data, so it never arises again.
   * Returns false if the current node is checked
   * and found to be impossible.
   * If the caller knows that the current node is already a nogood,
   * it would be a waste of time checking the current node,
   * as we are about to move up.
   * @param impossible_assignment An assignment PV->TV which should be erased
   * from ALL data.
   * @param action Should we include the current node, or ignore it? (If we're
   * about to move up, it's pointless to include it).
   * @return False if the current node becomes a nogood.
   */
  bool erase_impossible_assignment(
      std::pair<VertexWSM, VertexWSM> impossible_assignment,
      ImpossibleAssignmentAction action);

 private:
  NodesRawData& m_raw_data;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
