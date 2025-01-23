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
#include <map>
#include <set>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
class SearchBranch;

/** The basic idea is simple: during a search, if:
 * (1) we can QUICKLY calculate
 * a guaranteed lower bound on the additional scalar product which would
 * arise if the partial solution were to be completed, WHATEVER new
 * assignments are made; and if
 * (2) that would make the final scalar product larger than in the current
 * best solution;
 * ...THEN we can prune the search, and hence save a lot of time.
 *
 * The PROBLEM is: if (2) doesn't hold, then all of this calculation
 * is pointless. (1) is quite fast, but not free.
 *
 * Thus, we wrap within ANOTHER estimator (this class); we estimate
 * (very quickly) the lower bound roughly, and only if it seems likely
 * to detect a nogood do we actually do the calculation.
 *
 * Uses simple heuristics and feedback control to decide
 * when to activate the weight nogood detector,
 * and when to skip (trying to guess if the extra calculation
 * is worth it).
 *
 * It doesn't try to TIME or estimate AMOUNT of calculation
 * by the detector, or use fancy graph theory. It doesn't really
 * know anything about detection itself.
 * It just looks at results: a calculation was done, and a nogood
 * detected or not.
 * If not, we know by how much the lower bound fell short
 * of the required amount, so we can have "OK" or "bad" failures.
 *
 * Detecting a nogood is a SUCCESS; it prunes, and we ASSUME
 * that the time saved by pruning is more than the time spent calculating
 * to detect the nogood originally.
 *
 * Not detecting a nogood is a FAILURE. We spent time calculating,
 * looking for a nogood but not detecting one (whether or not
 * it really is a nogood; we just don't know at this stage).
 *
 * Really, almost ALL calculations (looking for node reduction,
 * nogood detection, etc. etc.) should be wrapped inside a manager object.
 * It would be great to have a general AI-type manager with a few simple
 * adjustable (or self-adjusting!) parameters to control such things
 * automatically; but that's probably very hard to do well.
 */
class WeightNogoodDetectorManager {
 public:
  explicit WeightNogoodDetectorManager(WeightWSM total_p_edge_weights);

  /** Decide, based upon the current data and previous results,
   * whether to try detecting or not.
   * If this returns true, the caller MUST follow the advice and activate
   * the detector, and immediately report back the result (see the "register"
   * functions).
   * @param current_weight The total weight (scalar product)
   * @param max_weight The total allowed weight of any solution
   * @param current_sum_of_p_edge_weights The sum of p-weights currently
   * assigned
   * @param current_number_of_assigned_vertices How many p-vertices are
   * assigned.
   * @param current_number_of_unassigned_vertices How many p-vertices are
   * unassigned.
   * @return True if the weight nogood detector should be used; if so, the
   * caller MUST report the result immediately using a "register" function.
   */
  bool should_activate_detector(
      WeightWSM current_weight, WeightWSM max_weight,
      WeightWSM current_sum_of_p_edge_weights,
      std::size_t current_number_of_assigned_vertices,
      std::size_t current_number_of_unassigned_vertices);

  // The caller is responsible for calling the weight nogood detector,
  // then MUST immediately report back what happened.

  /** The detector did detect a nogood (i.e., it found a guaranteed lower
   * bound which was so big that we immediately could backtrack).
   */
  void register_success();

  /** The detector calculated a lower bound L, but it was too small,
   * so it did not detect a nogood. (I.e., whether or not it really IS a
   * nogood, we couldn't be sure; so we had to proceed with the search).
   * @param current_weight The total weight so far
   * @param max_weight The allowed maximum total weight
   * @param extra_weight_lower_bound A weight L, such that the total weight
   * would definitely increase by at least L for any valid assignment.
   */
  void register_lower_bound_failure(
      WeightWSM current_weight, WeightWSM max_weight,
      WeightWSM extra_weight_lower_bound);

 private:
  const WeightWSM m_total_p_edge_weights;

  // The below values of the parameters are little more than
  // (slightly intelligent) GUESSES, with a very small amount of
  // experimental evidence. It would be great to have some really good theory
  // to pick sensible values; but probably no such theory exists today.

  // "pk" numbers are "per K", where K=1024,
  // rather than "per 100" (ordinary percentages).
  // This is preferable to a power of ten because
  // (1) it can give faster int division
  // (arbitrary int division apparently can be surprisingly slow),
  // (2) allows lossless calculations in most cases.

  // These adjustable controls specify the behaviour;
  // but experimentation is needed.
  struct ControlParameters {
    unsigned min_weight_pk_to_activate = 80;
    unsigned max_weight_pk_to_activate = 1024;
    // m_min_weight_pk_to_activate will decrease upon success,
    // meaning that our lower bounds are good, so try triggering them earlier.
    unsigned success_pk_to_activate_growth_factor_pk = 320;
    unsigned ok_failure_pk_to_activate_growth_factor_pk = 1400;
    unsigned bad_failure_pk_to_activate_growth_factor_pk = 1600;

    // If the number of assigned vertices drops below a certain amount,
    // reset the history as we've backtracked a lot.
    unsigned drop_below_n_assigned_vertices_reset_pk = 256;

    // Don't allow resets again until the number of assigned vertices
    // has risen above a certain amount.
    unsigned rise_above_n_assigned_vertices_allow_reset_pk = 800;

    unsigned min_final_weight_estimate_pk = 800;
    unsigned max_final_weight_estimate_pk = 1280;

    // If we detected a nogood, it's working; so allow more detector activation,
    // i.e. decrease final_weight_estimate_pk_to_activate.
    unsigned final_weight_estimate_pk_success_growth_pk = 720;
    // If we failed to detect a nogood, but only just,
    // reduce detector activation a bit.
    unsigned final_weight_estimate_pk_ok_failure_growth_pk = 1200;
    unsigned final_weight_estimate_pk_bad_failure_growth_pk = 1600;

    unsigned bad_failure_skip = 2;
    unsigned ok_failure_skip = 2;
  };

  struct State {
    // Only activate if the weight so far is above
    // a minimum fraction of the max weight.
    unsigned min_weight_pk_to_activate = 512;

    // Only activate if a crude estimate of the final weight
    // is above a certain amount of the max weight.
    unsigned final_weight_estimate_pk_to_activate = 1024;

    // "Resetting" the state is for when the conditions of the current
    // search node have changed so much that we should regard it as
    // a new, fresh environment.
    bool can_reset = false;
    unsigned remaining_skips = 0;
  };

  ControlParameters m_parameters;
  State m_state;

  void clamp_state_values();
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
