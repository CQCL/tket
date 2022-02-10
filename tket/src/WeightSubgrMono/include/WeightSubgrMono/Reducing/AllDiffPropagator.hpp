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

class SearchNodeWrapper;

/** Some pv->tv assignments have been made. This means that
 * no other pattern vertices are allowed to map to any of these tv.
 * Thus, we must both check for inconsistencies (duplicated tv)
 * and reduce domains (erase tv from every other domain).
 * This is not completely simple, because erasing TV from
 * some domain Dom(v) may reduce Dom(v) = {y} to size 1,
 * meaning that a new assignment v->y must then be made,
 * and recursively propagated.
 */
class AllDiffPropagator {
 public:
  /** Go through all chosen_assignments of the node which
   * have not yet been processed (according to
   * number_of_assignments_previously_processed_in_this_node);
   * update that variable;
   * remove the TV from all domains.
   *
   * Return false if ever a domain becomes empty
   * (so, we've encountered a nogood).
   *
   * If ever a domain reduces to size 1,
   * delete the PV key and move the PV->TV assignment
   * across to "assignments".
   * @param assignments All assignments made; will be updated
   * @param node_wrapper Contains the search node to be updated.
   * @param number_of_assignments_previously_processed_in_this_node The
   * assignments are stored in a vector, which has new assignments added to the
   * end; this keeps track of how many have been processed so far.
   * @return False if an inconsistency occurs (so, our search has hit an
   * impossible state - a "dead end" - and must backtrack).
   */
  bool reduce(
      Assignments& assignments, SearchNodeWrapper& node_wrapper,
      std::size_t& number_of_assignments_previously_processed_in_this_node)
      const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
