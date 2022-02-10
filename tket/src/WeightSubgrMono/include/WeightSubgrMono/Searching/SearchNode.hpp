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

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Raw data about current assignments and future possibilities.
 * A search branch consists of a list (logically, a stack) of search nodes
 * which we can move up and down (logically, pushing and popping).
 * The nodes know about domains for variables and current total weights,
 * but none of the extra bookkeeping needed to reduce/change domains.
 */
struct SearchNode {
  /** The total weight (scalar product) over all edges currently assigned. */
  WeightWSM current_scalar_product;

  /** The sum of the pattern edge weights which have currently been assigned.
   * (To decide between incomplete solution, we may care about which one has
   * most assigned total weight).
   */
  WeightWSM total_p_edge_weights;

  /** All (pv, tv) pairs which have been locked in.
   * Thus, the first element represents the choice that was made
   * to construct this search node initially.
   * All subsequent ones are inevitable (and thus, should NOT be listed
   * in the nogoods).
   */
  std::vector<std::pair<VertexWSM, VertexWSM>> chosen_assignments;

  /** The KEY is a pattern vertex x.
   * The VALUE is Dom(x), i.e. all possible target vertices which x
   * may be mapped to. (Always, Dom(x) has size >= 2).
   * TODO: store std::sets elsewhere, and REUSE, replacing this
   * with a std::map<VertexWSM, size_t> and lazy copy-on-write, etc. etc.
   * Operations on std::map<VertexWSM, size_t> should be a lot faster
   * because size_t objects are very cheap to copy.
   */
  PossibleAssignments pattern_v_to_possible_target_v;

  /** Interesting for testing; return log10 of
   * the total search space size (with naive exhaustive search:
   * i.e., the product of the domain sizes).
   * Of course, there is only a very slight connection between this number
   * and the hardness of a particular WSM problem.
   * This knows nothing about "rigidity", i.e. once one assignment is made,
   * how many more assignments inevitably follow with minimal calculation
   * on average? Seems like there should be interesting theory behind this,
   * as yet mostly undiscovered.
   */
  double get_log10_search_space_size() const;

  /** For testing, a string representation. */
  std::string str() const;

  /** The new assignment pv->tv, ASSUMED to be valid, has been made,
   * to a previous search node.
   * Then, copy over and initialise the data in THIS node.
   * However, does NOT remove tv from the domains of all variables
   * (this is left to a reducer object elsewhere).
   * @param pattern_v A pattern vertex pv.
   * @param target_v A pattern vertex tv.
   * @param previous_node The state of a previous node, from which we construct
   * this node by making the assignment pv->tv.
   */
  void initialise_from_assignment(
      VertexWSM pattern_v, VertexWSM target_v, const SearchNode& previous_node);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
