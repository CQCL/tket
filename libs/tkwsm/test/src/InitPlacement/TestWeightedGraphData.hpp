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
#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
class RNG;
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

// Returns data for a random graph GUARANTEED to use exactly the vertices
// {0,1,2,...,num_vertices-1}, with a mixture of weights.
// Note that the final number of edges, and min/max weights,
// are likely but not guaranteed to hold, because we don't bother
// to fiddle with the random selection to guarantee them.
GraphEdgeWeights get_graph_data(
    RNG& rng, unsigned num_vertices, unsigned approx_edges,
    WeightWSM approx_min_weight, WeightWSM approx_max_weight);

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
