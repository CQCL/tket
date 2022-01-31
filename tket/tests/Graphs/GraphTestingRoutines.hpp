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

#include <cstddef>
#include <utility>
#include <vector>

namespace tket {
namespace graphs {

class AdjacencyData;
struct GraphColouringResult;

namespace tests {

/** Simple routines to test graphs. */
struct GraphTestingRoutines {
  /**
   * Check that the colouring is valid, even if suboptimal
   * (without knowing the graph it came from).
   */
  static void require_valid_colouring_without_graph(
      const GraphColouringResult& colouring_result,
      bool require_no_colour_gaps = true);

  /**
   * As well as checking the colouring object alone,
   * also check that it is valid for the initial specified number of edges.
   * We do this because we might have generated a colouring only for
   * some of the edges.
   */
  static void require_valid_colouring(
      const GraphColouringResult& colouring_result,
      const std::vector<std::pair<std::size_t, std::size_t>>& edges,
      std::size_t number_of_edges_to_check);

  /**
   * Check that the colouring is valid for the graph (even if suboptimal).
   */
  static void require_valid_suboptimal_colouring(
      const GraphColouringResult& colouring_result,
      const AdjacencyData& graph_data);

  /**
   * Recalculate a colouring for this graph, check it is valid,
   * AND check that it doesn't use more colours than the given colouring
   * (although it might use fewer colours; the given colouring
   * is allowed to be suboptimal).
   */
  static void calculate_and_check_optimal_colouring(
      const GraphColouringResult& known_possibly_suboptimal_colouring,
      const AdjacencyData& graph_data);
};

}  // namespace tests
}  // namespace graphs
}  // namespace tket
