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
#include <map>
#include <vector>

#include "Utils/RNG.hpp"

namespace tket {
namespace graphs {
namespace tests {

/**
 * For testing purposes only, not of much independent interest
 * (and definitely an inefficient implementation).
 * Start with an NxN grid of squares; merge regions one-by-one across
 * a square edge, to make planar regions; the four colour theorem
 * then guarantees that a four colouring is possible.
 * (Still one of the longest mathematical proofs ever found,
 * requiring a computer!)
 */
class RandomPlanarGraphs {
 public:
  /** Make a fresh new grid. */
  void reset(std::size_t width = 20);

  /**
   * Remove a single dividing edge between two 1x1 squares,
   * possibly making fewer regions (i.e. two different regions become merged).
   *
   * @param rng The random number generator to provide the randomness.
   * @return The number of regions after the operation
   * (which thus forms a DECREASING sequence).
   */
  std::size_t merge_squares(RNG& rng);

  /**
   * This is really just the edges of the planar graph, but the vertices
   * are regions, not squares.
   * @return Element i is for region i, and lists all other regions which touch
   * it.
   */
  std::vector<std::vector<std::size_t>> get_region_data() const;

 private:
  /** Element i is the region to which 1x1 square i belongs. */
  std::vector<std::size_t> m_region_ids;

  /** The width, N, of the NxN square grid. */
  std::size_t m_width = 0;

  std::size_t m_number_of_regions = 0;

  /** Call each edge of the 1x1 square a "gate".
   * This tells us the possible "edges" of the planar graph
   * (the regions being the vertices).
   */
  struct Gate {
    std::size_t vertex1;
    std::size_t vertex2;
  };

  std::vector<Gate> m_remaining_gates;

  mutable std::vector<std::vector<std::size_t>> m_region_data;
  mutable std::map<std::size_t, std::size_t> m_old_id_to_new_id;

  /**
   * Used when computing m_region_data.
   * Register the fact that the two squares touch each other.
   */
  void register_touching_regions(
      std::size_t square1, std::size_t square2) const;
};

}  // namespace tests
}  // namespace graphs
}  // namespace tket
