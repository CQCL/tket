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

/** Simple input parameters to configure the solving. */
struct MainSolverParameters {
  /** The timeout only for the current call to "solve",
   * not the cumulative time.
   */
  long long timeout_ms;

  /** If set to TRUE, it will break off early when a complete solution
   * is found, no matter how poor (but still count as UNFINISHED,
   * in case the caller wants to continue looking for more solutions).
   */
  bool terminate_with_first_full_solution;

  /** For internal testing; set an upper bound on
   * the number of search iterations.
   */
  std::size_t max_iterations;

  /** If non-null, it means that we artificially impose
   * an upper bound on the total weight for a full solution,
   * i.e. an extra constraint on the solution.
   * This could speed up solving, BUT may prevent us
   * from finding any solution.
   */
  std::optional<WeightWSM> weight_upper_bound_constraint;

  MainSolverParameters();

  /** Just set the timeout in milliseconds; the most common parameter. */
  explicit MainSolverParameters(long long timeout_ms);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
