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
//#include <optional>
#include <set>
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** This is concerned entirely with domain reductions.
 * More complex logic involving searching, backtracking, etc.
 * and associated extra bookkeeping data
 * is the responsibility of other classes.
 */
class NodeWSM {
 public:
  NodeWSM();

  const PossibleAssignments& get_possible_assignments() const;

  /** This also checks for new assignments. */
  void set_possible_assignments(PossibleAssignments possible_assignments);

  // std::optional<VertexWSM> get_target_vertex(VertexWSM pattern_vertex) const;

  /** Lists (pv,tv) for pv->tv in order of assignment, but ONLY in this node. */
  const std::vector<std::pair<VertexWSM, VertexWSM>>& get_new_assignments()
      const;

  void set_scalar_product(WeightWSM scalar_product);
  WeightWSM get_scalar_product() const;

  void set_total_pattern_edge_weights(WeightWSM new_weight);
  WeightWSM get_total_pattern_edge_weights() const;

  bool alldiff_reduce(std::size_t n_assignments_processed);

  /** Doesn't check for validity (new_domain being a subset of the existing
   * one). But, takes care of new assignments.
   */
  void overwrite_domain(VertexWSM pattern_v, std::set<VertexWSM> new_domain);

  /** Overwrite with values in a vector instead of a set. Note that the vector
   * need not be sorted, and can contain duplicates.
   */
  void overwrite_domain(
      VertexWSM pattern_v, const std::vector<VertexWSM>& new_domain);

  void clear_new_assignments();

  /** Without any checks, simply set Dom(pv) = {tv}
   * and update the new assignments list.
   */
  void force_assignment(const std::pair<VertexWSM, VertexWSM>& assignment);

  struct ErasureResult {
    bool valid;
    bool assignment_was_possible;
  };

  /** Erase tv from Dom(pv). */
  ErasureResult erase_assignment(
      const std::pair<VertexWSM, VertexWSM>& assignment);

  enum class SoftErasureResult {
    TV_WAS_NOT_PRESENT,
    TV_WAS_ERASED,
    TV_REMAINS
  };

  /** Erase tv from Dom(pv), but ONLY if it would not make a new assignment.
   * Returns TRUE if tv is no longer present (i.e., either it was erased,
   * OR it was not present initially).
   */
  SoftErasureResult attempt_to_erase_assignment(
      const std::pair<VertexWSM, VertexWSM>& assignment);

  std::string str() const;

 private:
  // friend class SearchBranch;

  /** The total weight (scalar product) over all edges currently assigned. */
  WeightWSM m_scalar_product;

  /** The sum of the pattern edge weights which have currently been assigned.
   * (To decide between incomplete solution, we may care about which one has
   * most assigned total weight).
   */
  WeightWSM m_total_p_edge_weights;

  /** The KEY is a pattern vertex x.
   * The VALUE is Dom(x), i.e. all possible target vertices which x
   * may be mapped to.
   * We store ALL Dom(x), including those of size 1.
   */
  PossibleAssignments m_pattern_v_to_possible_target_v;

  std::vector<std::pair<VertexWSM, VertexWSM>> m_new_assignments;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
