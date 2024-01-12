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
#include <optional>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

/** This contains all the input config parameters for initial qubit placement;
 * HOPEFULLY set to reasonable defaults, but more testing and experimentation
 * needs to be done.
 */
struct IQPParameters {
  // Very crude, improve! Should be set according to number of PV, TV,
  // edges, etc. etc.
  unsigned max_wsm_iterations = 10000000;
};

/** Using all the config IQPParameters, actually calculate a solution. */
struct IQPResult {
  /** There are various strategies for computing the pattern graph weights;
   * see the struct PatternGraphData.
   * But they're irrelevant to this class. Just directly pass in
   * the pattern graph weights.
   */
  IQPResult(
      const GraphEdgeWeights& pattern_graph_weights,
      const GraphEdgeWeights& target_architecture_with_error_weights,
      unsigned timeout_ms, const IQPParameters& parameters);

  /** This will be sorted (i.e., in increasing order of PV).
   * In a sense, this is the only relevant output data.
   * However, note that isolated pattern vertices
   * (i.e., logical qubits with no 2-qubit gates acting on them)
   * will not occur here.
   */
  std::vector<std::pair<VertexWSM, VertexWSM>> initial_qubit_placement;

  // For testing, can be useful to have the timings and other information,
  // broken down into the different phases.

  unsigned long long mcct_time_ms = 0;
  unsigned mcct_iterations = 0;
  WeightWSM mcct_scalar_product = 0;

  unsigned long long wsm_init_time_ms = 0;
  unsigned long long wsm_solve_time_ms = 0;
  unsigned wsm_number_of_pruned_tv = 0;
  unsigned wsm_number_of_pruned_t_edges = 0;
  unsigned wsm_iterations = 0;
  std::optional<WeightWSM> wsm_scalar_product_opt;

  unsigned long long total_time_ms = 0;
};

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
