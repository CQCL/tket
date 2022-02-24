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

#include <map>

#include "Placement.hpp"
#include "TokenSwapping/SwapFunctions.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/RNG.hpp"

namespace tket {

/**
 * @brief Given a placement map generates `n` nearby placement maps.
 *
 * Based on an architecture and a placement map, generates random
 * placements that can be achieved with `m` swaps.
 *
 * Optionally uses token swapping optimisations to try to ensure
 * that the generated placements cannot be obtained in less than `m`
 * swaps, but this cannot be guaranteed.
 */
class NeighbourPlacements {
 public:
  using SwapVec = std::vector<Swap>;
  using NodeSwap = std::pair<Node, Node>;
  using NodeSwapVec = std::vector<NodeSwap>;
  struct Result {
    qubit_mapping_t map;
    NodeSwapVec swaps;
  };
  using ResultVec = std::vector<Result>;

  /**
   * @brief Construct a new Swap Placement object.
   *
   * @param arc The architecture defining the allowed swaps.
   * @param init_map The initial Qubit => Node map.
   */
  NeighbourPlacements(const Architecture& arc, const qubit_mapping_t& init_map);

  /**
   * @brief Generate `n` distinct placement maps using `dist` swaps for each map
   *
   * The sequences of swaps are generated randomly. Note that it cannot be
   * guaranteed that the generated placement cannot be obtained in less than
   * `dist` swaps. When optimise=true (default), we attempt to simplify
   * chains of swaps to make it more likely that `dist` swaps are indeed
   * necessary for the generated placement maps.
   *
   * If optimise=true, it is also possible that placements `dist` swaps away
   * do not exist. `max_tries` controls the number of attempts to generate
   * placements.
   *
   * If it is impossible (or very hard) to generate `n` distinct placement maps
   * of distance `dist` swaps away, then this method will raise a warning
   * and return fewer results and/or results with less than `dist` swaps.
   *
   * @param dist The number of swaps allowed on the architecture.
   * @param n The number of placement maps to generate (default n=1).
   * @param optimise Simplify the generated swap sequences (default true).
   * @param seed Seed for random number generator (default seed=5489).
   * @param max_tries Number of tries before aborting placement map generation
   *                  (default max_tries=10).
   * @return ResultVec A vector of the generated maps and swaps
   */
  ResultVec get(
      unsigned dist, unsigned n = 1, bool optimise = true, unsigned seed = 5489,
      unsigned max_tries = 10);

 private:
  // generate a single Result
  Result gen_result(
      unsigned dist, bool optimise = true, unsigned max_tries = 10);
  // generate a single swap
  Swap gen_swap();
  // apply swap list to init_map and return new placement map
  Result convert_to_res(const SwapVec& swaps);

  Architecture arc_;
  qubit_mapping_t init_map_;
  boost::bimap<unsigned, Node> u_to_node_;
  RNG rng_;
};

}  // namespace tket