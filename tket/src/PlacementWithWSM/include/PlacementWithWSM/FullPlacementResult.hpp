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
#include "PlacementAndStatistics.hpp"

namespace tket {

/** An internal class; the current strategy is a multistep process;
 * FIRST try to solve with adding few extra target graph edges;
 * then only if that fails, add more target graph edges
 * with higher weights to make a solution possible (even if bad).
 */
struct FullPlacementResult {
  /** Includes the placement, and extra information about it. */
  PlacementAndStatistics result;

  enum class Pass { INITIAL, COMPLETE_TARGET_GRAPH };

  /** The solver is multipass. We begin trying to solve the problem
   * with the given data, and if it fails, we successively relax
   * the constraints (weight constraints, or adding extra target edges)
   * until a solution IS found. This gives the pass at which it occurred.
   */
  Pass pass = Pass::INITIAL;

  /** Only the best result from all passes is retained,
   * but of course more passes might have occurred.
   */
  unsigned number_of_passes = 0;

  /** For testing purposes and reproducibility, we want to be able to
   * terminate a solve at an exact point, if it takes too long.
   * We cannot use timeouts for this.
   * Instead, use the number of iterations.
   */
  std::size_t iterations_for_pass = 0;

  /** How long, in milliseconds, did the initialisation steps take, in total?
   * (Useful for testing. Ideally initialisation steps should be very fast,
   * because most calculations can be made dynamic and lazy, i.e. you only
   * calculate if PV->TV may be possible the first time you try to make
   * the assignment; this way, many PV->TV assignments might never need
   * to be hecked, because they'll be pruned by the search strategy
   * before they occur).
   */
  long long total_init_time_ms = 0;

  /** How long, in milliseconds, did the main search steps take, in total? */
  long long total_search_time_ms = 0;

  /** Parameters to configure the algorithms. */
  struct Parameters {
    /** The timeout in milliseconds. */
    unsigned timeout_ms = 10000;

    /** For repeatability in tests, if desired to terminate early,
     * specify the pass and number of iterations.
     */
    std::optional<std::pair<Pass, std::size_t>> pass_data_opt;

    /** If not null, a constraint on how many search iterations to perform.
     * Of course the meaning is algorithm-dependent so this is really only
     * useful for testing purposes.
     */
    std::optional<std::size_t> max_iterations_opt;
  };

  /** An empty constructor. */
  FullPlacementResult();

  /** The main constructor; does the full solving and fills the data.
   * @param pattern_graph The original pattern graph, including edge weights.
   * @param original_target_graph The original target graph with edge weights,
   * but WITHOUT any added edges.
   * @param enlarged_target_graph The original_target_graph, but with some extra
   * edges and weights added.
   * @param gates List of the gates (more precisely, the logical qubits, i.e.
   * pattern vertices, involved in each single gate), in time order, which was
   * used to construct the pattern graph and weights.
   * @param parameters Extra options to control the algorithm.
   */
  FullPlacementResult(
      const WeightedSubgraphMonomorphism::GraphEdgeWeights& pattern_graph,
      const WeightedSubgraphMonomorphism::GraphEdgeWeights&
          original_target_graph,
      // "enlarged_target_graph" should be original_target_graph
      // with extra edges and weights added.
      const WeightedSubgraphMonomorphism::GraphEdgeWeights&
          enlarged_target_graph,
      const std::vector<std::set<WeightedSubgraphMonomorphism::VertexWSM>>&
          gates,
      const Parameters& parameters = {});

  /** For testing, returns a string with the data.
   * @param print_times If true, the timings in milliseconds are included.
   * @return A string representation of the data.
   */
  std::string str(bool print_times = false) const;
};

}  // namespace tket
