// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Input parameters to configure the solving. All set to sensible defaults. */
struct MainSolverParameters {
  /** The timeout in milliseconds, but only for the current call to "solve",
   * NOT the cumulative time.
   */
  long long timeout_ms;

  /** Mainly for internal testing; the same as timeout_ms,
   * but counts the actual raw number of ADDITIONAL search iterations
   * during the "solve" function (NOT cumulative number of iterations),
   * rather than the time. Thus, will stop as soon as the extra iterations
   * hit this limit, giving exactly repeatable and portable results.
   */
  std::size_t iterations_timeout;

  /** If set to TRUE, it will break off early when a complete solution
   * is found, no matter how poor (but still count as UNFINISHED,
   * in case the caller wants to continue looking for more solutions).
   */
  bool terminate_with_first_full_solution;

  /** If set to 0 (the default), ignore;
   * then at most one solution (the best one) will be stored
   * (and guaranteed to be optimal if the solve runs to completion).
   *
   * Otherwise, we will store ALL full solutions found, regardless
   * of how good or bad they are,
   * and terminate as soon as this many different full solutions
   * are stored. However, we make no guarantees about the order of solutions.
   */
  std::size_t for_multiple_full_solutions_the_max_number_to_obtain;

  /** If non-null, it means that we artificially impose
   * an upper bound on the total weight for a full solution,
   * i.e. an extra constraint on the solution.
   * This could speed up solving, BUT may prevent us
   * from finding any solution.
   */
  std::optional<WeightWSM> weight_upper_bound_constraint;

  /** For domain initialisation, how far out should we filter?
   * This is a balancing act: higher values will lead to smaller
   * initial domains and hence faster searching, but will also
   * need more computation; the SUM of initialisation
   * and search times is relevant.
   */
  unsigned max_distance_for_domain_initialisation_distance_filter;

  /** This is used during search; as for distance filtering,
   * it's a balancing act how large to make it.
   */
  unsigned max_distance_for_distance_reduction_during_search;

  /** Just set the timeout in milliseconds; the most common parameter. */
  explicit MainSolverParameters(long long timeout_ms = 1000);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
