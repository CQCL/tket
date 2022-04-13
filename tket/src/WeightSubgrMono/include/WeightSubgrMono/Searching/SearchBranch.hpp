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
#include <memory>

#include "../Reducing/DerivedGraphsReducer.hpp"
#include "../Reducing/HallSetReducer.hpp"
#include "AssignmentChecker.hpp"
#include "EnrichedNode.hpp"
#include "WeightUpdater.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DistancesReducer;
class NeighboursData;
class WeightChecker;

/** Represents a depth-first traversal.
 * However, it's really about reducing nodes etc.,
 * not the search logic. This contains the list of nodes so far,
 * which may be viewed as a stack, but the caller is responsible for
 * moving up and down the list.
 */
class SearchBranch {
 public:
  /** Constructs node[0], i.e. the initial node.*/
  SearchBranch(
      PossibleAssignments initial_pattern_v_to_possible_target_v,
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
      DistancesReducer& distances_reducer);

  /** Extra parameters needed to configure how we reduce a single node
   * with all the extra components stored in the SearchComponents object.
   */
  struct ReductionParameters {
    WeightWSM max_weight;
    unsigned max_distance_reduction_value;
  };

  /** This is at the branch level, not node level, BECAUSE this branch has
   * m_outstanding_impossible_assignments, which nodes don't know about.
   * If it returns false, the node is in an invalid state; but that's OK
   * because the caller should immediately discard it and move up.
   */
  bool reduce_current_node(const ReductionParameters& parameters);

  /** Keep moving up one level and reducing, until EITHER it's reduced
   * (returning true), OR we can't move up any more (and return false).
   */
  bool backtrack(const ReductionParameters& parameters);

  /** ASSUMING that the current node has been reduced,
   * make the given assignment [which must be valid],
   * and move down to a new node.
   * @param p_vertex A pattern vertex, pv.
   * @param t_vertex A target vertex, tv. We are making the new assignment
   * pv->tv.
   */
  void move_down(VertexWSM p_vertex, VertexWSM t_vertex);

  /** Inform the branch that pv->tv is (and always was)
   * impossible (regardless of the state of the search,
   * i.e. this would be safe to share with other searches even at the start).
   * The assignment will be removed from all domains, so it can never occur.
   * However, IGNORE the current node.
   */
  void register_impossible_assignment(
      const std::pair<VertexWSM, VertexWSM>& assignment);

  const NodeWSM& get_current_node() const;

  /** This is the set of vertices adjacent to any existing assigned vertex,
   * to be used as candidates when choosing a new variable to assign.
   */
  std::set<VertexWSM>& get_current_node_candidate_variables();

  std::set<VertexWSM> get_used_target_vertices() const;

  void activate_weight_checker(WeightWSM total_p_edge_weights);

 private:
  const NeighboursData& m_pattern_ndata;
  const NeighboursData& m_target_ndata;
  DistancesReducer& m_distances_reducer;
  DerivedGraphsReducer m_derived_graphs_reducer;
  AssignmentChecker m_assignment_checker;

  /** To save reallocation, we DON'T resize m_enriched_nodes;
   * instead, we just remember the index and overwrite/resize as appropriate.
   */
  std::size_t m_enriched_nodes_index;

  const WeightUpdater m_weight_updater;
  HallSetReducer m_hall_set_reducer;

  std::vector<EnrichedNode> m_enriched_nodes;

  std::unique_ptr<WeightChecker> m_weight_checker_ptr;

  /** As long as this is nonempty, will keep trying to erase these assignments
   * completely from the nodes, until they're all gone.
   */
  std::vector<std::pair<VertexWSM, VertexWSM>>
      m_outstanding_impossible_assignments;

  std::vector<VertexWSM> m_outstanding_impossible_target_vertices;

  NodeWSM& get_current_node_nonconst();

  /** Reduces the current node aggressively, returning false if it becomes
   * impossible. */
  bool attempt_to_clear_outstanding_impossible_data();
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
