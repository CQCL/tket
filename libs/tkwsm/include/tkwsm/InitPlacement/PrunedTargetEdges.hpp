// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
class NeighboursData;
namespace InitialPlacement {

/** Experiment to get good default values! */
struct TargetEdgePruningParameters {
  /** As a fraction x/1024, how many additional target edges
   * should we allow at most?
   * So, e.g., if our input solution uses 14 target edges,
   * and we have the value 512,
   * then our final target graph will have up to
   *   14*(1 + (512/1024)) = 21  edges.
   */
  unsigned max_additional_number_of_target_edges_factor_per_kilo = 1024;

  /** As a fraction x/1024, what is the minimum number of target edges
   * which we want to leave still unused at the end, as a fraction
   * of the input unused edges?
   * So, e.g., if our input solution uses 18 target edges, but there are
   * 12 target vertices, giving a maximum of 12*11/2 = 66
   * possible target edges, then 66-18 = 48 potential target edges are
   * currently unused by the input solution.
   * So if the varaible value were 256, say, then since 256/1024 = 1/4,
   * then we would ensure that at least 48*(1/4) = 12 target edges were
   * not part of the new target graph, i.e. there could be at most
   * 66-12 = 54 target edges in the new graph.
   * Thus, since our input solution uses 18 target edges, we could add
   * up to  54-18=36 new target edges (but we might add fewer).
   */
  unsigned min_implicit_unused_number_of_target_edges_factor_per_kilo = 500;
};

/** Once we've found some solution (not necessarily optimal) to a WSM
 * problem with COMPLETE target graph, remove some unused target edges
 * to get a new WSM problem which is still soluble (because we've GIVEN
 * a solution - although maybe not optimal).
 * It's important to note that the explicit_target_ndata only gives the
 * explicit target edges, NOT all the edges (since the target graph is
 * COMPLETE). All additional possible target edges not explicitly mentioned
 * are taken to have the given implicit weight.
 *
 * Again, this needs more experiments and far more test problems
 * to tweak the algorithm, as currently (in combination with MCCT)
 * it doesn't give useful WSM problems.
 */
GraphEdgeWeights get_new_target_graph_data(
    const NeighboursData& pattern_ndata,
    const NeighboursData& explicit_target_ndata,
    WeightWSM implicit_target_weight,
    // assigned_target_vertices[pv] = tv, where pv->tv is used in the input
    // solution.
    const std::vector<unsigned>& assigned_target_vertices,
    const TargetEdgePruningParameters& parameters = {});

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
