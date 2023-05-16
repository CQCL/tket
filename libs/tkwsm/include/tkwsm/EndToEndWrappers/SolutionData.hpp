// Copyright 2019-2023 Cambridge Quantum Computing
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
#include "SolutionWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** These are mainly useful for testing. It is important that
 * they are all cheap to calculate.
 */
struct ExtraStatistics {
  std::size_t number_of_pattern_vertices = 0;
  std::size_t number_of_target_vertices = 0;

  /** Count how many PV->TV assignments were "maybe possible" at the start,
   * i.e. not excluded by the domain initialisation (although we may
   * subsequently discover that some are, and always were, impossible).
   */
  std::size_t initial_number_of_possible_assignments = 0;

  /** When the weight nogood detector is initialised with the set of all used
   * tarbet vertices, record their number here.
   */
  std::optional<std::size_t> n_tv_initially_passed_to_weight_nogood_detector;

  /** How many target vertices were still under consideration by
   * the weight nogood detector?
   */
  std::optional<std::size_t> n_tv_still_valid_in_weight_nogood_detector;

  /** How many PV->TV assignments were actually carried out
   * during the search? This is of course counting a subset of those
   * counted in initial_number_of_possible_assignments.
   * The smaller this number is, the better the pruning which took place.
   */
  std::size_t total_number_of_assignments_tried = 0;

  /** How many PV->TV assignments were excluded during the search,
   * which domain initialisation had not originally excluded?
   * These are of course a subset of those
   * counted in initial_number_of_possible_assignments.
   * However, there may some some PV->TV which are counted in NEITHER
   */
  std::size_t total_number_of_impossible_assignments = 0;

  /** Occasionally, a target vertex TV is found to be impossible;
   * NOTHING can actually map to it, even though some initial domains
   * included it. This is very rare, but record them here.
   */
  std::vector<VertexWSM> impossible_target_vertices;
};

struct SolutionData {
  /** If true, the search is over;
   * EITHER we've found a (joint) OPTIMAL solution,
   * OR we've proved that there is NO solution.
   * If false, then our solution (if any) is merely the best
   * found so far, not necessarily optimal.
   * But, if the "terminate_with_first_full_solution"
   * solver option was chosen and a solution was found,
   * so that it terminated early, this will still be set to "false",
   * so that the caller can continue searching for more solutions
   * if desired.
   */
  bool finished = false;

  /** The total cumulative search time in milliseconds. */
  long long search_time_ms = 0;

  /** The initialisation time in milliseconds. */
  long long initialisation_time_ms = 0;

  // If upper/lower bounds for the scalar product are equal,
  // then it's effectively the standard pure subgraph monomorphism problem;
  // there's no point in evaluating scalar products for different solutions
  // because they're all equal,
  // and no point using a weight nogood detector.
  //
  // But if we allow zero weights, then this is NOT quite equivalent to
  // "the pattern weights are constant, and the target weights are constant".
  // e.g. consider p-weights = [0,1,1] and t-weights = [1,1,1].
  // Even though the p-weights are NOT constant,
  // any valid assignment will give scalar product 2 exactly,
  // and (in a simple case like this) it will calculate L=U=2.

  /** A simple lower bound for the total weight (scalar product) any complete
   * valid solution would have. But we might not know if a solution exists!
   */
  WeightWSM trivial_weight_lower_bound = 0;

  /** A simple upper bound for the total weight
   * that any full valid solution can have.
   * (But, we might not know if a solution exists!)
   */
  WeightWSM trivial_weight_initial_upper_bound = 0;

  /** The total number of search iterations taken. */
  std::size_t iterations = 0;

  /** Does the target graph contain every possible edge?
   * Obviously, if the target graph is complete
   * (at least, for all the target vertices mentioned in the edge weights)
   * then no graph theoretic vertex filtering is possible;
   * EVERY tv lies in every PV domain.
   */
  bool target_is_complete = false;

  /** Simply store the sum of all p edge weights. */
  WeightWSM total_p_edge_weights;

  std::vector<SolutionWSM> solutions;

  ExtraStatistics extra_statistics;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
