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
#include <optional>
#include <set>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;
class NeighboursData;

/** Some new pv->tv assignments have been made in a search node.
 * Updates the total weight (scalar product) of this node,
 * whenever new edges become assigned.
 * This does NOT have to be attached to a specific node,
 * BUT the caller must keep track of how many assignments in THIS node
 * were already processed in a previous call (similar to reducing).
 */
class WeightCalculator {
 public:
  struct Result {
    WeightWSM scalar_product;
    // It's convenient to return the additional edge weights,
    // rather than the total.
    WeightWSM total_extra_p_edge_weights;
  };

  /** Calculates the new scalar product and total p-edges weight
   * from all new p-edges which have become assigned
   * (i.e., both end vertices are now assigned).
   * @param pattern_ndata Data about the pattern graph.
   * @param target_ndata Data about the target graph.
   * @param accessor Object to retrieve information about the current search
   * node (inclusing new assignments).
   * @param number_of_processed_assignments The number of assignments stored at
   * the start of the new assignments vector, of the current node, which this
   * class has already processed (and thus, need not process again).
   * @param max_scalar_product The maximum scalar product which is allowed; if
   * we exceed this, we break off immediately: it's a nogood.
   * @return Null if the current weight exceeds the maximum weight,
   * OR we detect a graph-theoretic nogood (i.e., ignoring the weights).
   */
  std::optional<Result> operator()(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
      const DomainsAccessor& accessor,
      std::size_t number_of_processed_assignments,
      WeightWSM max_scalar_product) const;

 private:
  // Necessary to avoid double counting edges,
  // since newly_assigned_p_vertices is not sorted and thus not
  // efficiently searchable.
  // A p-edge is only added if exactly one end vertex was already seen.
  // If both were seen, it would have been added
  // when the first vertex was seen.
  mutable std::set<VertexWSM> m_p_vertices_seen;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
