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

struct FixedData;

/** Various simple checks on the weights, and other inconsistencies,
 * including that the weight addition and multiplication
 * does not overflow. (If they do, then an IntegerOverflow exception
 * is thrown).
 */
struct CheckedWeightBounds {
  // If these lower, upper bounds
  // are null then no solution is possible.
  // Except for trivial cases (e.g. empty pattern graph or full
  // target graph), we cannot know if a solution exists
  // without doing the full search.
  // However we CAN do simple estimates on the possible weight
  // of any solution.
  //
  // These would be useful in a solving strategy which begins with
  // a very strict imposed upper bound weight constraint on any solution,
  // in the hope that it speeds up the search:
  // if the imposed constraint is possible, it should find the optimal
  // solution faster because higher-weight solutions will be pruned early.
  // If instead the imposed constraint is impossible, then (if weight pruning
  // is good) it would hopefully be able to prove quickly that no solution
  // satisfying the weight constraint exists.

  /** If non-null, a guaranteed lower bound on the best possible
   * solution weight, IF a solution exists (but it may not).
   */
  std::optional<WeightWSM> lower_bound;

  /** If non-null, a guaranteed upper bound on the weight of
   * any solution, IF any exists (but it may not).
   */
  std::optional<WeightWSM> upper_bound;

  /** If true, then no solution exists. */
  bool other_inconsistency_occurred;

  /*
  // The most accurate, but also expensive, check.
  CheckedWeightBounds(const FixedData& fixed_data,
      const GraphEdgeWeights& pattern_edges);
  */

  /**We do progressively more expensive but more accurate checks,
   * trying to get (upper bound * the safety factor) not to overflow.
   */
  CheckedWeightBounds(
      const FixedData& fixed_data, const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges, WeightWSM extra_safety_factor);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
