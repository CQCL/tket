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

#include "WeightSubgrMono/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
struct SolutionWSM;
}

/** Given a complete or partial solution to a WSM problem,
 * and the original data about the gates, computes some statistics
 * and ensures that a VALID assignment is produced (i.e., no two p-vertices
 * being assigned to the same t-vertex).
 */
struct PlacementAndStatistics {
  // Mainly results for two qubit gates (other gates are mostly ignored).

  /** Number of 2-qubit gates where qubits were assigned to
   * adjacent vertices in the original architecture.
   */
  unsigned n_gates_with_original_edges = 0;

  /** Number of 2-qubit gates where qubits were assigned to
   * nearby vertices in the original architecture, i.e. they are adjacent
   * in the target graph used for the WSM problem (which is the original
   * architecture, plus extra edges).
   */
  unsigned n_gates_with_some_token_swapping = 0;

  /** If n_gates_with_some_token_swapping is nonzero,
   * the sum of the total target edge weights just for those gates.
   */
  WeightedSubgraphMonomorphism::WeightWSM total_weights_with_token_swapping = 0;

  /** The number of 1-qubit gates; these are just ignored. */
  unsigned single_qubit_gates = 0;

  /** The number of n-qubit gates with n>2; these are ignored
   * when selecting the best solution, BUT involved in the initial
   * problem because they are artificially decomposed into
   * a sequence of simultaneous 2-qubit gates.
   */
  unsigned n_many_qubit_gates = 0;

  /** The number of n-qubit gates with n>2, with at least 1 qubit
   * being unassigned.
   */
  unsigned n_many_qubit_gates_unassigned = 0;

  /** Number of 2-qubit gates where the qubits were assigned,
   * but to far away vertices (i.e., they are not adjacent
   * in the target graph, even though this has added extra edges to
   * the original architecture).
   * Therefore, so much token swapping is involved that our simple
   * cost model is probably not very meaningful for them.
   */
  unsigned n_poor_gates = 0;

  /** Number of 2-qubit gates where at least one qubit
   * was assigned.
   */
  unsigned n_unassigned_gates = 0;

  /** PV->TV assignments which are at least valid (even if poor).
   * Note that, UNLIKE SolutionWSM, these are actually checked
   * and ensured to be valid. Invalid assignments in SolutionWSM
   * are simply discarded.
   * (A better algorithm would try to find the best possible
   * consistent assignments, maybe using bipartite matching,
   * but this just uses the first assignment seen).
   */
  std::map<
      WeightedSubgraphMonomorphism::VertexWSM,
      WeightedSubgraphMonomorphism::VertexWSM>
      valid_assignments;

  PlacementAndStatistics();

  /** Compute and fill all the data.
   * @param pattern_graph The original input for the WSM problem.
   * @param original_target_graph The original architecture, with no extra edges
   * added (so, NOT the input for the WSM problem).
   * @param enlarged_target_graph The original_target_graph with extra edges and
   * weights added, used as in input for a WSM problem. Thus adjacent vertices
   * in this graph are regarded as "close", if not adjacent in the original
   * graph.
   * @param gates The original time-sliced gate interactions, in order.
   * @param solution The best solution to a WSM problem (maybe only partial).
   */
  PlacementAndStatistics(
      const WeightedSubgraphMonomorphism::GraphEdgeWeights& pattern_graph,
      const WeightedSubgraphMonomorphism::GraphEdgeWeights&
          original_target_graph,
      const WeightedSubgraphMonomorphism::GraphEdgeWeights&
          enlarged_target_graph,
      const std::vector<std::set<WeightedSubgraphMonomorphism::VertexWSM>>&
          gates,
      const WeightedSubgraphMonomorphism::SolutionWSM& solution);

  /** This is a bit subjective, where partial solutions are concerned;
   * if neither this nor the other has found a complete solution,
   * which should we choose to return?
   */
  bool prefer_other_solution(const PlacementAndStatistics& other) const;

  /** A string representation, for debugging. */
  std::string str() const;
};

}  // namespace tket
