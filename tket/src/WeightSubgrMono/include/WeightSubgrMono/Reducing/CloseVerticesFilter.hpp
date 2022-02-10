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

#include "../GraphTheoretic/CloseNeighboursCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
class SearchNodeWrapper;

/** For a new pv->tv assignment, check the vertices near pv
// and require their targets to be close to tv,
// i.e. reduce the domains of all pattern vertices
// within a small distance of the vertex.
 */
class CloseVerticesFilter {
 public:
  CloseVerticesFilter();

  /** Go through and check that no duplicate target vertices occurred.
   * Also, remove them from the domains.
   * NOTE: we only need to handle the current node
   * with its current assignments and domains,
   * because of symmetry. If pv1->tv1 conflicted with pv2->tv2,
   * and pv1->tv1 was in a previous node,
   * then tv2 should already have been removed from Dom(pv2)
   * in a previous node.
   * Returns false if a nogood is detected.
   * @param fixed_data The data (graph theoretic, etc. calculated initially and
   * valid for all time).
   * @param assignments The up-to-date assignments over the COMPLETE search
   * branch (not just the ones made in this search node).
   * @param number_of_assignments_previously_processed_in_this_node How many new
   * assignments, just in THIS node, were previously checked by this class? It
   * will possibly also add some new assignments, but check all of them by the
   * time it returns.
   * @param node_wrapper An object to access the search node data.
   * @return False if it detects an impossible state (so the search must
   * backtrack).
   */
  bool reduce(
      const FixedData& fixed_data, Assignments& assignments,
      std::size_t number_of_assignments_previously_processed_in_this_node,
      SearchNodeWrapper& node_wrapper);

 private:
  unsigned m_num_levels;
  CloseNeighboursCalculator m_close_neighbours_calculator;

  // KEY: the vertex in a graph. VALUE: data about its nearby neighbours.
  std::map<VertexWSM, CloseNeighboursCalculator::CloseVertices> m_pattern_data;
  std::map<VertexWSM, CloseNeighboursCalculator::CloseVertices> m_target_data;

  std::vector<VertexWSM> m_reduced_domain;

  // We WANT to intersect the current domain of a variable
  // with a UNION of other containers, namely the other t-vertices
  // which are sufficiently close to a nearby vertex.
  // So, pass in the t-vertices one-by-one to check
  // if they're in the intersection.
  // This seems quite inefficient; we should improve.
  bool t_vertex_is_suitable(
      const CloseNeighboursCalculator::CloseVertices& close_target_vertices,
      VertexWSM original_target_vertex, unsigned current_p_level);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
