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
#include <set>

#include "WeightSubgrMono/GraphTheoretic/GeneralStructs.hpp"

namespace tket {

/** Used to convert the raw gate interaction data into the pattern graph with
 * weights, for a WSM problem.
 */
struct PatternGraphTimeSlices {
  /** Element i is for time i.
   * It simply lists all edges which are supposed to occur in the
   * interaction graph, at time i.
   * Note that edges are NOT required to be independent at each time
   * (multiqubit interactions, with >= 3 qubits at once,
   * are modelled simply as a list of 2-qubit interactions stuck together).
   * Also, the same edge can occur at multiple times
   * (as we may have multiple gates).
   * All edges (v1,v2) will have v1<v2.
   */
  std::vector<std::vector<std::pair<
      WeightedSubgraphMonomorphism::VertexWSM,
      WeightedSubgraphMonomorphism::VertexWSM>>>
      time_sliced_data;

  /** Pass in the interactions one-by-one, in time order.
   * @param gates_in_order For each gate, specify the set of all logical qubits
   * involved in that gate.
   */
  explicit PatternGraphTimeSlices(
      const std::vector<std::set<WeightedSubgraphMonomorphism::VertexWSM>>&
          gates_in_order);

  /** Takes into account "time decay"; a gate occurring at a later time
   * contributes smaller weights, because it's less likely that
   * the qubits will actually be on or close to their original assigned
   * physical qubits, as many extra swaps may have been added.
   * The final weight of an edge is the SUM of all the weights
   * at each time where that edge occurs in a gate.
   */
  struct WeightParameters {
    /** A single edge at time 0 has this much weight. */
    WeightedSubgraphMonomorphism::WeightWSM time_zero_edge_weight = 1000;

    /** A single edge at the final time has this much weight. */
    WeightedSubgraphMonomorphism::WeightWSM final_time_edge_weight = 200;

    /** Does simple (integer-valued) linear interpolation (PLUS checks
     * for overflow!), given the desired number of weights
     * (i.e., number of time slices).
     * Note that exponential decay would be nice, BUT is considerably
     * more complicated without using doubles and std::exp, std::log
     * (which would make results non-portable; we only want integer
     * operations).
     * @param size the size of the desired vector.
     * @return A vector of weights at each time, of given size, interpolating
     * between the first and last weights.
     */
    std::vector<WeightedSubgraphMonomorphism::WeightWSM>
    get_single_edge_weights(unsigned size) const;
  };

  /** Once the initial calculations have been done in the constructor,
   * complete the calculation to get the final pattern graph with weights.
   * @param parameters Parameters to control the calculation of the weights at
   * each time.
   * @return The raw pattern graph with weights, for a WSM problem.
   */
  WeightedSubgraphMonomorphism::GraphEdgeWeights get_weights(
      const WeightParameters& parameters) const;

  /** Use custom weights at each time, rather than generating them
   * by linear interpolation.
   * @param single_edge_weights_at_all_times A vector of weights at each time,
   * of the correct size.
   * @return The raw pattern graph with weights, for a WSM problem.
   */
  WeightedSubgraphMonomorphism::GraphEdgeWeights get_weights(
      const std::vector<WeightedSubgraphMonomorphism::WeightWSM>&
          single_edge_weights_at_all_times) const;
};

}  // namespace tket
