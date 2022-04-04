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

#include "../GraphTheoretic/DerivedGraphStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
class SearchNodeWrapper;
struct DerivedGraphs;

/** We take a new PV->TV assignment and use derived graphs to reduce the
 * domains of other vertices (e.g., neighbours of PV in a derived graph
 * must map to neighbours of TV in the corresponding derived target graph).
 */
class DerivedGraphsReducer {
 public:
  /** Go through the new assignments PV->TV and reduce all associated domains
   * using derived graphs.
   * @param fixed_data Contains, e.g. neighbours data for original pattern and
   * target graphs.
   * @param assignments All assignments pv->tv made so far.
   * @param number_of_assignments_previously_processed_in_this_node Used to
   * avoid processing older assignments multiple times.
   * @param node_wrapper The object containing the node to be updated.
   * @param derived_pattern_graphs The object containing the necessary data for the
   * derived p-graphs.
   * @param derived_target_graphs The object containing the necessary data for the
   * derived t-graphs.
   * @return True if the search is still ongoing, false if a nogood is detected
   * (so the search must backtrack).
   */
  bool reduce_domains(
      const FixedData& fixed_data, Assignments& assignments,
      std::size_t number_of_assignments_previously_processed_in_this_node,
      SearchNodeWrapper& node_wrapper, DerivedGraphs& derived_pattern_graphs,
      DerivedGraphs& derived_target_graphs);

 private:
  // A work vector.
  std::vector<VertexWSM> m_reduced_domain;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
