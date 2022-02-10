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
#include <set>

#include "WeightSubgrMono/GraphTheoretic/GeneralStructs.hpp"

namespace tket {

/** Responsible for adding extra weighted edges to the target graph.
 * Note that the weight of a target edge is regarded as directly proportional
 * to the probability of a 2-qubit gate on that edge giving an error.
 * Therefore, the total expected number of errors is proportional to
 * the sum of the edge weights, multiply counted,
 * each target edge being counted the number of times
 * it is used to apply a gate.
 */
struct TargetGraphData {
  /** A target graph with edges, which could be used in a WSM problem. */
  WeightedSubgraphMonomorphism::GraphEdgeWeights final_data;

  /** The sorted list of nonisolated target vertices (isolated vertices
   * are just discarded), since it's calculated anyway,
   * so we might as well store it here.
   */
  std::vector<WeightedSubgraphMonomorphism::VertexWSM> sorted_vertices;

  /** Parameters used to add extra edges and weights to the original target
   * graph.
   */
  struct Parameters {
    /** This is effectively +infinity for weights;
     * each edge weight will never be allowed to go beyond this.
     */
    std::optional<WeightedSubgraphMonomorphism::WeightWSM> max_edge_weight;

    /** If N>1, then a SWAP gate on 2 vertices is regarded as
     * composed of N single primitive 2-qubit gates,
     * each of which incurs the edge weight cost.
     */
    unsigned swap_gate_count = 3;

    /** Only add new edges where the original graph distance between the
     * vertices is <= this value.
     */
    unsigned max_path_length_for_new_edges = 5;

    /** Also used to constrain max_edge_weight,
     * by looking at the largest weight that already exists.
     * Don't allow any new weight to be more than this multiple
     * of the largest weight.
     */
    WeightedSubgraphMonomorphism::WeightWSM
        max_edge_weight_largest_weight_ratio = 100;

    /** Don't allow any new weight to be more than this multiple
     * of the smallest existing nonzero weight.
     */
    WeightedSubgraphMonomorphism::WeightWSM
        max_edge_weight_smallest_weight_ratio = 100000;

    /** If an added edge weight would go beyond the maximum, do we
     * cap it at the maximum (so that the edge still exists),
     * or not add the edge at all?
     */
    bool remove_high_edge_weights = true;

    /** Suppose that we have a very low-fidelity edge [v1,v2],
     * and we discover a longer path [v1,v3,...,v2] which is actually better
     * (would have a lower edge weight if we added the edge),
     * because it moves along high fidelity gates.
     * Should we REPLACE the original edge weight with the lower weight?
     */
    bool replace_low_fidelity_primitive_gates_with_longer_paths = true;
  };

  /** Pass in the initial target edges and weights, which SHOULD
   * correspond to actual edges and fidelity data;
   * and pass in parameters for computing extra edges and weights.
   */
  TargetGraphData(
      WeightedSubgraphMonomorphism::GraphEdgeWeights data,
      Parameters parameters);
};

}  // namespace tket
