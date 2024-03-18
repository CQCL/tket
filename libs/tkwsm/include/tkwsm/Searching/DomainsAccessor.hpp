// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <boost/dynamic_bitset.hpp>
#include <optional>
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct NodesRawData;
class NodesRawDataWrapper;

/** This is passed to Reducer, Checker, Updater objects,
 * to obtain information about the domains in the current node only,
 * and alter them.
 */
class DomainsAccessor {
 public:
  /** The wrapped NodesRawData object will be directly altered. */
  explicit DomainsAccessor(NodesRawDataWrapper& raw_data_wrapper);

  unsigned get_number_of_pattern_vertices() const;

  /** Every unassigned p-vertex (i.e., with size Domain(pv) > 1) in the
   * current node is included in here.
   * However, this may also include some assigned vertices.
   * @return A vector which definitely includes all pattern vertices which are
   * unassigned, in the current node. However, it may include other
   * pattern vertices.
   */
  const std::vector<VertexWSM>& get_unassigned_pattern_vertices_superset()
      const;

  bool current_node_is_valid() const;

  /** This may be a reference to the SAME object returned by
   * get_unassigned_pattern_vertices_superset(), or it may be different.
   * The caller is free to overwrite it at the END (after processing).
   * A vector::swap is most efficient, of course.
   * @return A vector which can be directly overwritten at the END, in order to
   * fill out all the unassigned vertices in the current node.
   */
  std::vector<VertexWSM>&
  get_unassigned_pattern_vertices_superset_to_overwrite();

  /** Return Domain(pv) in the current node, i.e. the set of all target
   * vertices which pv could be mapped to, as we extend the current mapping.
   * @param pv A vertex in the pattern graph.
   * @return The domain of pv in the current search node.
   */
  const boost::dynamic_bitset<>& get_domain(VertexWSM pv) const;

  /** The number of TV currently in Domain(PV). */
  std::size_t get_domain_size(VertexWSM pv) const;

  /** Returns true if Domain(pv) is different in the current search node
   * and the previous node. */
  bool domain_created_in_current_node(VertexWSM pv) const;

  /** The newly made assignments, just in the current node.
   * Note that this list might be cleared after they have been processed;
   * thus it need not be the complete list of ALL assignments which occurred
   * in this node.
   */
  const std::vector<std::pair<VertexWSM, VertexWSM>>& get_new_assignments()
      const;

  /** Once a domain is fully reduced and all assignments processed, clear the
   * data; there is no use for it any more.
   */
  void clear_new_assignments();

  /** Returns the scalar product  sum_e w(e).w(f(e))  over all p-edges e
   * in the which have actually been assigned, i.e. for which
   * both end vertices are assigned in the current node (although of course
   * they may have been assigned first in previous nodes).
   */
  WeightWSM get_scalar_product() const;

  /** Simply overwrite the scalar product in the current node with
   * the new value; the caller is assumed to know how to calculate it
   * correctly (which requires, of course, carefully considering the new
   * assignments in the current node, and adding the values appropriately).
   * @param scalar_product The new scalar product value; simply overwrites the
   * existing one.
   * @return This object, for chaining.
   */
  DomainsAccessor& set_scalar_product(WeightWSM scalar_product);

  /** Returns the sum of all pattern edge weights for those edges
   *  which have been assigned so far (i.e., both end vertices
   * have been assigned).
   * @return The sum of weights of all assigned pattern edges.
   */
  WeightWSM get_total_p_edge_weights() const;

  /** Simply overwrite the total sum of pattern edge weights in the current
   * node with the new value. The caller is responsible for doing this
   * calculation correctly.
   * @param value The new value; simply overwrites the existing one.
   * @return This object, for chaining.
   */
  DomainsAccessor& set_total_p_edge_weights(WeightWSM value);

  /** Assuming that the caller has already processed the given number
   * of new assignments in this current node (without clearing the
   * new assignments list - which stores them in order of creation),
   * go through ALL remaining assignments and process them, by applying
   * alldiff propagation. (I.e., if PV->y, then y must be erased from EVERY
   * other domain in this node).
   * The caller must keep track of how many new_assignments are processed.
   * @param n_assignments_already_processed The number of assignments in the
   * current new assignments list which were previously processed by this
   * function in the current node. The caller must keep track of this
   * information. After returning, EITHER all will have been processed, OR a
   * nogood is found.
   * @return False if some domain becomes empty (so, we are at a nogood: an
   * invalid node, meaning that we must backtrack in the search. We return early
   * as soon as this occurs, since there's no point in continuing with an
   * invalid node).
   */
  bool alldiff_reduce_current_node(std::size_t n_assignments_already_processed);

  struct IntersectionResult {
    ReductionResult reduction_result;

    // Only bother filling this in if it's NOT a nogood or new assignment.
    std::size_t new_domain_size;

    bool changed;
  };

  /** Simply intersect Domain(PV) with the specified set of TV, represented
   * by a bitset. For speed, uses the "swap" function (similar to vector::swap
   * and set::swap - time O(1), no copies) - which thus alters domain_mask.
   * We should have an intersect function WITHOUT doing this -
   * the operator -= function should work.
   * We check for and update new assignments if necessary.
   * @param pattern_v The pattern vertex.
   * @param domain_mask A set of target vertices to intersect with the domain;
   * will be altered.
   * @return The result of changing the domain.
   */
  IntersectionResult intersect_domain_with_swap(
      VertexWSM pattern_v, boost::dynamic_bitset<>& domain_mask);

 private:
  NodesRawData& m_raw_data;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
