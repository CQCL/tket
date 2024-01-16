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
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** A full solution to a weighted subgraph monomorphism problem (WSM).
 */
struct SolutionWSM {
  /** The (pv->tv) assignment pairs. If nonempty,
   * then this should be a complete valid solution.
   * Will be sorted with increasing PV.
   */
  std::vector<std::pair<VertexWSM, VertexWSM>> assignments;

  /** The total weight, i.e. sum w(e).w(f(e)) over all pattern edges e,
   * where f(e) is the corresponding target edge
   * (which exists by the definition of f being a subgraph
   * monomorphism).
   */
  WeightWSM scalar_product = 0;

  /** sum w(e) over all pattern edges e currently assigned.
   * Of course, all valid solutions should have the same value.
   */
  WeightWSM total_p_edges_weight = 0;

  /** For testing purposes.
   * Check thoroughly that the solution is valid (if complete),
   * return a helpful error message if not. Returns an empty string if no
   * errors. So, if complete and this solution returns no errors,
   * we can be sure that it's valid.
   * @param pattern_edges_and_weights Data for the pattern graph
   * @param target_edges_and_weights Data for the target graph
   * @return A string, empty if there are no errors.
   */
  std::string get_errors(
      const GraphEdgeWeights& pattern_edges_and_weights,
      const GraphEdgeWeights& target_edges_and_weights) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
