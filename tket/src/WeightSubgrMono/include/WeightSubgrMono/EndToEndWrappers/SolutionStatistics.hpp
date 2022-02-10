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

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Basic statistics about a solution and the current state of the search
 * after a single "solve" call, without details of the actual solution.
 */
struct SolutionStatistics {
  /** If true, the search is over;
   * either we've found the OPTIMAL solution,
   * or we've proved that there is NO solution.
   *
   * Note, that if the "terminate_with_first_full_solution"
   * solver option was chosen and a complete solution was found,
   * so that it terminated early, this will still be set to "false",
   * so that the caller can continue searching for more solutions
   * if desired.
   */
  bool finished = false;

  /** The total search time in milliseconds. */
  long long search_time_ms = 0;

  /** The total initialisation time in milliseconds (should be small). */
  long long initialisation_time_ms = 0;

  /** A simple lower bound for the total weight (scalar product) any complete
   * valid solution would have - although note that we might not know
   * if a solution actually exists!
   */
  WeightWSM trivial_weight_lower_bound = 0;

  /** The number of search iterations taken. */
  std::size_t iterations = 0;

  /** If not null, a simple upper bound for the total weight
   * that any valid solution can have.
   * (Of course, it cannot detect if solutions actually exist
   * or not).
   */
  std::optional<WeightWSM> trivial_weight_initial_upper_bound;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
