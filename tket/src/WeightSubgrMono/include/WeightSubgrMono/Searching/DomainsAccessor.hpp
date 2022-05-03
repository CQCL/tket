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

/** This is passed to Reducer, Checker, Updater objects,
 * to obtain information about the domains in the current node only,
 * and alter them.
 */
class DomainsAccessor {
 public:
  /** The wrapped NodesRawData object will be directly altered. */
  explicit DomainsAccessor(NodesRawDataWrapper& raw_data_wrapper);

  const std::vector<VertexWSM>& get_pattern_vertices() const;

  /** Every unassigned p-vertex (i.e., with size Domain(pv) > 1) in the
   * current node is included in here.
   * However, this may also include some assigned vertices.
   */
  const std::set<VertexWSM>& get_unassigned_pattern_vertices_superset() const;

  const std::set<VertexWSM>& get_domain(VertexWSM pv) const;

  /** Is Dom(pv) different in the previous node? */
  bool domain_created_in_current_node(VertexWSM pv) const;

  /** We prefer these vertices, when assigning.
   * But, this may be empty, in which case we'll have to look at all
   * unassigned vertices. However, also provide direct access.
   */
  std::set<VertexWSM>& get_candidate_vertices_for_assignment_nonconst();
  const std::set<VertexWSM>& get_candidate_vertices_for_assignment() const;

  /** The newly made assignments, just in the current node. */
  const std::vector<std::pair<VertexWSM, VertexWSM>>& get_new_assignments()
      const;

  /** Once a domain is fully reduced and all assignments processed, clear the
   * data. */
  void clear_new_assignments();

  WeightWSM get_scalar_product() const;
  DomainsAccessor& set_scalar_product(WeightWSM);

  WeightWSM get_total_p_edge_weights() const;
  DomainsAccessor& set_total_p_edge_weights(WeightWSM);

  std::string str(bool full = false) const;

  /** The caller should keep track of how many
   * new_assignments are processed. */
  bool alldiff_reduce_current_node(std::size_t n_assignments_already_processed);

  /** The new domain should be a subset of the old domain,
   * but this is not fully checked - only trivial O(1) checks.
   */
  ReductionResult overwrite_domain(
      VertexWSM pv, const std::set<VertexWSM>& new_domain);

  /** "new domain" has already been filled, and we are happy to alter it
   * (it is wrok data, i.e. "dummy" reusable data).
   * Thus, do set::swap instead of copying.
   */
  ReductionResult overwrite_domain_with_set_swap(
      VertexWSM pv, std::set<VertexWSM>& new_domain);

  /** Vectors are faster than sets for some purposes; we assume (but don't
   * check) that the vector elements are distinct, but unsorted.
   */
  ReductionResult overwrite_domain(
      VertexWSM pattern_v, const std::vector<VertexWSM>& new_domain);

  struct IntersectionResult {
    ReductionResult reduction_result;

    // Only bother filling this in if it's NOT a nogood or new assignment.
    std::size_t new_domain_size;

    bool changed;
  };
  /** We erase all of the given TV from Dom(PV). */
  IntersectionResult intersect_domain_with_complement_set(
      VertexWSM pattern_v,
      const std::set<VertexWSM>& forbidden_target_vertices);

  std::set<VertexWSM>&
  get_current_node_unassigned_pattern_vertices_superset_to_overwrite();

 private:
  NodesRawData& m_raw_data;

  std::set<VertexWSM>& get_domain_nonconst(VertexWSM pv);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
