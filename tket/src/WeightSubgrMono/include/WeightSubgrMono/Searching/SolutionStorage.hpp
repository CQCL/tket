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
#include <optional>
#include <set>
#include <string>

#include "SearchBranch.hpp"
#include "SolutionWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Store the best solution found so far (which might be partial). */
class SolutionStorage {
 public:
  SolutionStorage();

  /** Used for "squeezing": any solution, full or partial, MUST have
   * scalar product <= this, or it will be rejected.
   *
   * It can be decreased at any time, but if ever increased,
   * incorrect results (i.e., incorrectly cutting out some parts
   * of the search tree) may result.
   * @param weight The maximum allowed weight of any new solution.
   * @return This object, for chaining.
   */
  SolutionStorage& set_pruning_weight(WeightWSM weight);

  /** For testing, may be removed in future. If set to level>0,
   * may print debug information to std::cerr.
   * @param log_level A level; any number>0 may result in writing to std::cerr
   * @return This object, for chaining.
   */
  SolutionStorage& set_log_level(unsigned log_level);

  /** ASSUMES that the current assignments within the SearchBranch
   * form a complete solution.
   * @param branch A SearchBranch with a current complete list of assignments.
   * @return True if it was accepted (if it is better than the current solution,
   * so replaces it).
   */
  bool add_full_solution(const SearchBranch& branch);

  /** Should only be called with solutions KNOWN to be incomplete
   * (including, maybe, nogoods - don't worry in detail what to do
   * if the assignment is impossible, let the caller decide later).
   * Will never overwrite an existing full solution.
   * @param branch A SearchBranch with a current complete list of assignments.
   * @return True if it was accepted (if it is better than the current solution,
   * so replaces it).
   */
  bool add_partial_solution(const SearchBranch& branch);

  /** @return The best solution so far, stored within this class. */
  const SolutionWSM& best_solution() const;

  /** If not null, any solution, full or partial, MUST have scalar product
   * <= this value, or it will be rejected.
   * @return If not null, the maximum scalar product (weight) which a solution
   * can have, to be accepted.
   */
  std::optional<WeightWSM> get_acceptable_scalar_product() const;

 private:
  unsigned m_log_level;
  WeightWSM m_pruning_weight;
  SolutionWSM m_solution;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
