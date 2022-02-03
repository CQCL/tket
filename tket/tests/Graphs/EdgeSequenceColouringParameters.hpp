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
#include <limits>

namespace tket {
namespace graphs {
namespace tests {

struct EdgeSequence;
struct RandomGraphParameters;

/**
 * Parameters to specify how to run tests. We want to generate graphs
 * by adding edges one-by-one.
 * We then colour the graphs in sequence, and check that the colouring
 * is valid for EVERY previous graph, and also
 * that the number of colours used is an increasing sequence.
 * (which does not GUARANTEE that the colouring is optimal, but does
 * give some extra confidence that it might be).
 */
struct EdgeSequenceColouringParameters {
  /**
   * Colour only every k graphs, so value 1 means colour EVERY graph.
   * Checking a colouring is much cheaper than generating it,
   * so we test each calculated colouring against EVERY previous graph.
   */
  std::size_t step_skip_size = 1;

  /** Stop the test when you've done this many colourings. */
  std::size_t max_number_of_colourings =
      std::numeric_limits<std::size_t>::max();

  /** Stop the test when the graph has this many edges. */
  std::size_t max_number_of_edges = std::numeric_limits<std::size_t>::max();

  /**
   * This counts how many times the graph colouring function was called.
   * I think it's good practice to check the value at the end of a test.
   * The actual value is unimportant, just whether it changes between commits.
   * The reasons are:
   * (i) Check that we are actually testing something; it's all too easy
   *     to have bugs in test data which mean we're not testing anything,
   *     or as much as we thought.
   * (ii) If something changes unexpectedly, either due to the algorithms
   *      or test data, then we have a chance of spotting it.
   * (iii) Although it could give a false sense of security, it's always nice
   *       to see that a test is doing 1000 colourings, rather than just 5, say.
   */
  std::size_t total_number_of_colourings = 0;

  /**
   * Take the parameters for generating random graphs,
   * actually generate them, then test the colourings.
   * @param parameters An object with functions to generate random graphs
   *  of a certain type, edge-by-edge.
   * @param edge_sequence The object which actually calls the "parameters"
   * object to construct the edges, and store them.
   * @return The number of calculated colourings in this call. (Useful as an
   * extra check).
   */
  std::size_t test_colourings(
      RandomGraphParameters& parameters, EdgeSequence& edge_sequence);
};

}  // namespace tests
}  // namespace graphs
}  // namespace tket
