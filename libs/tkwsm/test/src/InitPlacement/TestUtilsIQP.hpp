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
#include "PlacementCostModelInterface.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

// The first is the actual cost of applying the gates.
// The second is the placement (logical qubit -> physical qubit pairs).
typedef std::pair<WeightWSM, std::vector<std::pair<VertexWSM, VertexWSM>>>
    CostedIQPSolution;

// Check that each element of "solutions" does indeed give a valid placement
// for the physical qubits contained within the PlacementCostModelInterface
// object, with the given cost, for the specified gates to apply.
void test_known_solutions(
    const std::vector<CostedIQPSolution>& solutions,
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates,
    PlacementCostModelInterface& placement_cost_model);

// Takes the given gates and PlacementCostModelInterface object,
// and runs IQP with WSM to find a placement.
// Note that the pattern_graph needs to be constructed by the caller
// using the given gates somehow - usually done with a "PatternGraphData"
// object.
// It will only place those PV mentioned in "gates".
void run_end_to_end_IQP_and_check_solution(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates,
    const GraphEdgeWeights& pattern_graph,
    PlacementCostModelInterface& placement_cost_model, WeightWSM cost,
    unsigned timeout_ms = 200, bool verbose = false);

// Given a list of vertices v0, v1, v2, ...,
// successively work out the paths from v[i] to v[i+1]
// and put them in a string, for easy copy/paste and manual inspection.
std::string get_many_paths_test_str(
    const PlacementCostModelInterface& graph,
    const std::vector<VertexWSM>& vertices);

std::string do_token_swaps_and_check_placements(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates,
    PlacementCostModelInterface& graph);

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
