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

#include "TestUtilsIQP.hpp"

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/InitPlacement/EndToEndIQP.hpp>
#include <tkwsm/InitPlacement/InputStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

void test_known_solutions(
    const std::vector<CostedIQPSolution>& solutions,
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates,
    PlacementCostModelInterface& placement_cost_model) {
  for (const auto& cost_and_solution : solutions) {
    placement_cost_model.initialise_with_qubit_placement(
        cost_and_solution.second);
    CHECK(placement_cost_model.get_cost(gates) == cost_and_solution.first);
  }
}

void run_end_to_end_IQP_and_check_solution(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates,
    const GraphEdgeWeights& pattern_graph,
    PlacementCostModelInterface& placement_cost_model, WeightWSM cost,
    unsigned timeout_ms, bool verbose) {
  if (verbose) {
    std::cerr << "\nRunning IQP with " << gates.size() << " gates, "
              << pattern_graph.size() << " pattern edges, expected cost "
              << cost << "; timeout " << timeout_ms;
  }
  // Just use the defaults - but should experiment!
  const IQPParameters iqp_parameters;

  const IQPResult iqp_result(
      pattern_graph, placement_cost_model.get_graph_data(), timeout_ms,
      iqp_parameters);

  // It should be sorted by PV.
  for (unsigned ii = 1; ii < iqp_result.initial_qubit_placement.size(); ++ii) {
    REQUIRE(
        iqp_result.initial_qubit_placement[ii - 1].first <
        iqp_result.initial_qubit_placement[ii].first);
  }
  placement_cost_model.initialise_with_qubit_placement(
      iqp_result.initial_qubit_placement);
  CHECK(placement_cost_model.get_cost(gates) == cost);
}

std::string get_many_paths_test_str(
    const PlacementCostModelInterface& graph,
    const std::vector<VertexWSM>& vertices) {
  std::stringstream ss;
  for (unsigned ii = 1; ii < vertices.size(); ++ii) {
    ss << graph.get_path_str(vertices[ii - 1], vertices[ii]);
  }
  return ss.str();
}

std::string do_token_swaps_and_check_placements(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates,
    PlacementCostModelInterface& graph) {
  std::stringstream ss;
  // We won't directly change this, although the reference will remain valid;
  // the token swaps will indirectly change the values.
  const auto& token_to_vertex_map = graph.get_current_placement();
  for (const auto& swap : gates) {
    const VertexWSM v1 = token_to_vertex_map.at(swap.first);
    const VertexWSM v2 = token_to_vertex_map.at(swap.second);
    const auto cost =
        graph.do_token_swapping_and_apply_gate(swap.first, swap.second);
    ss << "\nTOKEN Swap (" << swap.first << "," << swap.second
       << ") between vertices " << v1 << " " << v2 << "; cost " << cost
       << graph.get_path_str(v1, v2) << "\nNOW, placement: { ";
    for (const auto& entry : token_to_vertex_map) {
      ss << entry.first << "->" << entry.second << " ";
    }
    ss << "}\n";
  }
  return ss.str();
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
