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
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** A solution (maybe only partial) to a weighted
 * subgraph monomorphism problem (WSM).
 */
struct SolutionWSM {
  /** It may be incomplete or even empty
   * (in which case, the assignments might be invalid,
   * and the total p edge weights and scalar product
   * are not so meaningful - more of a guide,
   * only as a last resort for the caller).
   */
  bool complete;

  /** The total weight, i.e. sum w(e).w(f(e)) over all pattern edges e
   * which have currently been assigned, where f(e) is the corresponding
   * target edge (which exists by the definition of f being a subgraph
   * monomorphism).
   */
  WeightWSM total_scalar_product_weight;

  /** sum w(e) over all pattern edges e currently assigned. */
  WeightWSM total_p_edges_weight;

  /** All the (pv->tv) assignment pairs, in order of being made.
   * NOT checked for validity, i.e. some tv may be duplicated;
   * the CALLER should decide what to do with incomplete/invalid solutions.
   */
  std::vector<std::pair<VertexWSM, VertexWSM>> assignments;

  SolutionWSM();

  /** For testing.
   * Checks thoroughly that the solution is valid (if complete),
   * returns a helpful error message if not. Returns an empty string if no
   * errors. So, if complete and this solution returns no errors, we can be sure
   * it's valid.
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
