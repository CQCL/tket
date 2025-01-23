// Copyright Quantinuum
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
#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
class NeighboursData;
struct VertexRelabelling;
namespace InitialPlacement {

WeightWSM get_scalar_product_with_complete_target(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    WeightWSM implicit_target_weight,
    // assignments[PV] = TV
    const std::vector<unsigned>& assignments);

GraphEdgeWeights get_relabelled_graph_data(
    const GraphEdgeWeights& graph_data, const VertexRelabelling& relabelling);

/** Check if some of the weight values are so large that integer overflow
 * might occur, when calculating the scalar product for a pattern->target
 * embedding with complete target graph. If so, throw.
 * This should not occur unless you have stupidly large weights,
 * but if this does occur then you simply have to reduce the weights.
 *
 * Otherwise, returns a crude upper bound. Of course, no way to know
 * the ACTUAL largest possible value without solving the WSM problem.
 *
 * @param pattern_ndata Data about the pattern graph.
 * @param explicit_target_ndata Data about the EXPLICIT target graph.
 * @param implicit_target_weight For each (v1,v2) of distinct target vertices,
 * if the edge is not mentioned in the EXPLICIT target graph, then it is taken
 * to exist, with this edge weight.
 * @param extra_safety_factor When converting unsigned to signed int types, we
 * want the unsigned value to fit within the signed value. So require that
 * multiplying by this extra factor still does not overflow.
 * @return A guaranteed upper bound for the scalar product (but maybe far too
 * big). Throws if it cannot get such a bound.
 */
WeightWSM get_scalar_product_upper_bound_for_complete_target_graph(
    const NeighboursData& pattern_ndata,
    const NeighboursData& explicit_target_ndata,
    WeightWSM implicit_target_weight, WeightWSM extra_safety_factor = 3);

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
