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
#include <set>

#include "../GraphTheoretic/DerivedGraphsCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
class SearchNodeWrapper;

class DerivedGraphsContainer;

/** We take a new PV->TV assignment and use derived graphs to reduce the
 * domains of other vertices (e.g., neighbours of PV in a derived graph
 * must map to neighbours of TV in the corresponding derived target graph).
 */
class DerivedGraphsReducer {
 public:

  /** Go through the new assignments PV->TV and reduce all associated domains
   * using derived graphs.
   * @param fixed_data Contains, e.g. neighbours data for original pattern and target graphs.
   * @param assignments All assignments pv->tv made so far.
   * @param number_of_assignments_previously_processed_in_this_node Used to avoid processing older assignments multiple times.
   * @param node_wrapper The object containing the node to be updated.
   * @param derived_graphs The object containing the necessary data for the derived graphs.
   * @return True if the search is still ongoing, false if a nogood is detected (so the search must backtrack).
   */
  bool reduce_domains(
      const FixedData& fixed_data, Assignments& assignments,
      std::size_t number_of_assignments_previously_processed_in_this_node,
      SearchNodeWrapper& node_wrapper,
      DerivedGraphsContainer& derived_graphs);

 private:

  std::vector<VertexWSM> m_reduced_domain;

  bool reduce_with_derived_weighted_neighbours(
          const DerivedGraphsCalculator::NeighboursAndCounts& pattern_neighbours,
          const DerivedGraphsCalculator::NeighboursAndCounts& target_neighbours,
          Assignments& assignments,
          SearchNodeWrapper& node_wrapper);

  // Within the derived graph D3, look at distance-2 vertices,
  // calculated lazily.

  std::set<VertexWSM> m_work_set;

  std::map<VertexWSM, std::vector<VertexWSM>> m_pattern_neighbours_at_distance_two_in_d2;
  std::map<VertexWSM, std::vector<VertexWSM>> m_target_neighbours_at_distance_two_in_d2;
  std::map<VertexWSM, std::vector<VertexWSM>> m_pattern_neighbours_at_distance_two_in_d3;
  std::map<VertexWSM, std::vector<VertexWSM>> m_target_neighbours_at_distance_two_in_d3;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
