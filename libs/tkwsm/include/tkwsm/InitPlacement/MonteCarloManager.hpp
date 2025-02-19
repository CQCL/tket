// Copyright Quantinuum
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
#include <tkrng/RNG.hpp>

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

struct MonteCarloManagerParameters {
  // How many iterations do we demand as a minimum before we can reset?
  unsigned min_iterations_for_change = 20;

  // There are 2 types of progress possible
  // (i.e., a changed solution is accepted):
  //  "weak progress" means that we have just decreased the current cost,
  //    EVEN THOUGH the cost might still be higher than the
  //    best one found so far.
  // "record breaker" means that we really have found a new cost lower than
  //    the previous all-time lowest cost.

  // The "per kilo fractions" are x/1024 for an integer x
  // (i.e., "per 1024" rather than "per 100", for percentages).
  // This is because dividing by 1024 might be a bit faster than
  // dividing by 100 (and definitely won't be slower!)
  // (Apparently int division can be surprisingly slower than
  // int multiplication or addition).

  // How many extra iterations without any progress do we allow,
  // as a fraction of the existing number, before we demand a reset?
  unsigned per_kilo_fraction_of_allowed_extra_iterations_without_weak_progress =
      500;

  // How many extra iterations without any new reocrd breaker do we allow,
  // as a fraction of the existing number, before we demand a reset?
  unsigned
      per_kilo_fraction_of_allowed_extra_iterations_without_record_breakers =
          1000;

  unsigned max_runs_without_record_breaking = 10;
  unsigned max_runs_without_progress = 10;
};

/** For use in MonteCarloCompleteTargetSolution.
 * Even for a simplified jumping algorithm, simpler than simulated annealing
 * because we only allow improving moves, it's still unclear how to choose
 * the best parameters.
 * The difficulty of simulated annealing and other similar random jumping
 * algorithms is NOT in the algorithm itself; it is in the fully automatic
 * selection of parameter values and termination criteria,
 * which still doesn't seem to have a fully satisfactory solution.
 */
class MonteCarloManager {
 public:
  /** It should be clear by now that all these values are rather arbitrary!
   * We need to do more experiments to find the best values.
   */

  explicit MonteCarloManager(
      const MonteCarloManagerParameters& parameters =
          MonteCarloManagerParameters());

  enum class Action {
    CONTINUE_WITH_CURRENT_SOLUTION,
    RESET_TO_NEW_SOLUTION,
    TERMINATE
  };

  /** We've just improved upon our current solution
   * (which may or may not be a new record breaker).
   * What should we do now?
   */
  Action register_progress(WeightWSM new_cost, unsigned iteration);

  /** At this iteration number, we have not made any progress
   * (weak progress OR record breaking).
   */
  Action register_failure(unsigned iteration);

 private:
  const MonteCarloManagerParameters m_parameters;

  unsigned m_runs_without_record_breaking;
  unsigned m_runs_without_progress;

  // What is the lowest weight found so far?
  WeightWSM m_best_cost;

  unsigned m_next_iteration_to_reset_if_no_progress;
  unsigned m_next_iteration_to_reset_if_no_record_breaker;

  void update_after_reset(unsigned iteration);

  void update_after_weak_progress(unsigned iteration);
};

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
