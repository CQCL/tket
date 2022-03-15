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

#include "../WeightPruning/WeightNogoodDetectorManager.hpp"
#include "EnrichedNode.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
struct SharedData;

/** Represents a depth-first traversal.
 * However, it's really about reducing nodes etc.,
 * not the search logic. This contains the list of nodes so far,
 * which may be viewed as a stack, but the caller is responsible for
 * moving up and down the list.
 */
class SearchBranch {
 public:

  SearchBranch();

  /** Fills in the top node, but might not FULLY reduce.
   * @param shared_data Contains data shared across ALL search branches (we may
   * be running several "in parallel", i.e. interleaved), which may include
   * previously searched parts of the tree, etc. etc.
   * @return What happens when we reduce the top node; does it give a complete
   * solution, or an inconsistency (meaning that the search is finished?)
   */
  ReductionResult initialise(SharedData& shared_data);

  /** For the search logic, it's necessary to know if we are at
   * the very first search or not.
   */
  bool move_down_has_been_called() const;

  /** Move up one level, but do NOT reduce the node.
   * Returns false if finished, i.e. the search is OVER!
   */
  bool backtrack();

  /** ASSUMING that the current node has been fully reduced,
   * make the given assignment [which must be valid; hence this cannot fail]
   * and move down one level to a new node.
   * The new node, however, will not yet be reduced.
   * @param p_vertex A pattern vertex, pv.
   * @param t_vertex A target vertex, tv. We are making the new assignment
   * pv->tv.
   */
  void move_down(VertexWSM p_vertex, VertexWSM t_vertex);

  /** Reduce the current node as much as possible,
   * until it is valid (or invalid), i.e. not in a "partial" state;
   * if successfully reduced, we are confident that
   * there are no edge/vertex clashes, the current weight
   * is up to date,...
   * Does NOT take any further action with full/partial solutions, though.
   * @param shared_data The data shared across ALL search branches.
   * @param max_weight The maximum permitted weight (scalar product) for any
   * solution.
   * @return What happened with the reduction.
   */
  ReductionResult reduce_current_node(
      SharedData& shared_data, WeightWSM max_weight);

  /** Simply returns the read-only accessor for the current search node. */
  const SearchNodeWrapper& get_current_node_wrapper() const;

  typedef std::vector<EnrichedNode> EnrichedNodes;

  /** Returns all the data for the search so far, and the current level
   * by reference (i.e., index of current enriched node;
   * NOT the number of nodes, which is level+1).
   * The CALLER can decide what to do with partial/complete solutions, etc.
   * @param level The current level, returned by reference.
   * @return The current list of all search nodes. NOTE that there may be extra
   * unused data on the end, which should be ignored; it would be wasteful to
   * keep resizing and deallocating/reallocating memory for search nodes. This
   * is why the current level is also returned.
   */
  const EnrichedNodes& get_data(std::size_t& level) const;

  /** Returns the current assignments pv->tv in the search.
   * @return All current pv->tv assignments made in this search branch down to
   * and including the current node.
   */
  const Assignments& get_assignments() const;

  /** Should be rare; returns a writeable version of the assignments.
   * @return All current pv->tv assignments made in this search branch,
   * non-const.
   */
  Assignments& get_assignments_mutable();

  /** If we've decided that PV->TV is an impossible assignment,
   * try to erase TV from every Dom(PV). However, for simplicity,
   * do not introduce extra assignments, i.e. if Dom(PV) = { TV, u }
   * at some level, leave Dom(PV) unchanged, as otherwise we'd need
   * more complicated code to force PV=u; it will naturally take care
   * of itself in time. Returns true if all trace of PV->TV is removed
   * from all domains.
   * @param pv A pattern vertex.
   * @param tv A target vertex.
   * @return True if tv has been successfully erased from every Dom(pv), at every level.
   */
  bool erase_assignment(VertexWSM pv, VertexWSM tv);


 private:
  std::size_t m_level;
  EnrichedNodes m_enriched_nodes;

  // KEY: PV, VALUE: TV
  // Assignments across the WHOLE branch (not just the current node).
  // Note that we will maintain the AllDiff constraint,
  // always deleting TV from domains immediately;
  // so we NEVER have to check all assignments, only
  // the domains of the unassigned variables.
  Assignments m_assignments;

  bool m_move_down_has_been_called = false;

  // A workset.
  std::set<VertexWSM> m_values_assigned_in_this_node;

  WeightNogoodDetectorManager m_weight_nogood_detector_manager;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
