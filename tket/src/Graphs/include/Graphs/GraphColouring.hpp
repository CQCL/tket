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
#include <string>
#include <vector>

namespace tket {
namespace graphs {

class AdjacencyData;

/**
 * The calculated colouring for a graph.
 */
struct GraphColouringResult {
  /**
   * The number of colours used; if nonzero, it should be correct.
   */
  std::size_t number_of_colours = 0;

  /**
   * For testing/debugging, a printable string version of the result.
   */
  std::string to_string() const;

  /**
   * The calculated colours. Element i is the colour of vertex i.
   * The colours will be 0,1,2,...,m.
   */
  std::vector<std::size_t> colours;

  GraphColouringResult();

  /**
   * Automatically sets the number_of_colours variable.
   * @param colours The calculated colouring, colours[i] being the colour of
   * vertex i.
   */
  explicit GraphColouringResult(const std::vector<std::size_t>& colours);
};

/**
 * It's expected that more routines will be added over time!
 */
struct GraphColouringRoutines {
  /**
   * The main end-to-end colouring function.
   * @param adjacency_data The graph to be coloured.
   */
  static GraphColouringResult get_colouring(
      const AdjacencyData& adjacency_data);
};

}  // namespace graphs
}  // namespace tket
