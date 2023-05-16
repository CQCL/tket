// Copyright 2019-2023 Cambridge Quantum Computing
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
 *
 * Don't think of the nodes as "histories".
 * Instead, think of them as states still to be explored
 * (which happen to be explored in a particular order, like traversing
 * a tree depth first; only needs linear size. Other orderings are possible,
 * but more complicated to implement, needing some kind of priority queue,
 * and could need exponential size).
 * Thus, when we move down from a node by assigning PV->TV,
 * we REMOVE TV from Domain(PV) before creating a new node.
 * This is simpler than trying to remember which assignment was made
 * and removing TV next time we return to the node.
 * Thus the node left behind is NOT a previous state, but is rather
 * a new CANDIDATE state to be further reduced.
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
  boost::dynamic_bitset<> get_used_target_vertices() const;

  /** Keep popping back the internal nodes data UNTIL it reaches
   * an apparently valid node (i.e., a "good" node, i.e. a node which is
   * NOT nogood); returns false if it cannot.
   * However, does not reduce the node (so it may still end up being invalid
   * after we perform reductions).
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

  /** We've newly discovered that PV->TV is always impossible.
   * Erase it from all data, so it never arises again.
   * Does not check if the current node becomes a nogood.
   * @param impossible_assignment An assignment PV->TV which will be erased
   * from ALL data.
   */
  void erase_impossible_assignment(
      std::pair<VertexWSM, VertexWSM> impossible_assignment);

 private:
  NodesRawData& m_raw_data;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
