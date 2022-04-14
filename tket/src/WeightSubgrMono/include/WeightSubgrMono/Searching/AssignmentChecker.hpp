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

class DerivedGraphsReducer;
class DistancesReducer;

/** Each time an assignment PV->TV is made, do various graph theoretic checks
 * to try to find a contradiction, i.e. prove that it is actually impossible.
 * Note that we DON'T cache negative results:
 * the caller is responsible for erasing impossible assignments
 * as soon as they occur; thus, they will be erased from ALL data,
 * and never arise again. (And thus, not be checked twice).
 */
class AssignmentChecker {
 public:
  /** Note that many checks also come from reducer objects.
   * Usually, a check (i.e. returning yes/no for a specific PV->TV,
   * but without altering domains) is cheaper than a reduce
   * (i.e., reducing Domain(PV) by erasing TV which are incompatible
   * with ALL the assignments made so far, including other PV).
   * Thus we do the checks here, but store the reducer objects elsewhere
   * to be used later.
   */
  AssignmentChecker(
      DerivedGraphsReducer& derived_graphs_reducer,
      DistancesReducer& distances_reducer);

  /** Positive results are stored, so it won't waste time checking twice
   * if it's already passed.
   */
  bool operator()(
      const std::pair<VertexWSM, VertexWSM>& assignment,
      unsigned distances_reducer_max_dist);

 private:
  DerivedGraphsReducer& m_derived_graphs_reducer;
  DistancesReducer& m_distances_reducer;

  /** Every pv->tv assignment in here has previously passed all
   * the simple checks for validity, so we need not repeat.
   * (Also, we don't need to store the impossible ones,
   * since the caller should erase them from all data as soon as they fail
   * the checks).
   */
  PossibleAssignments m_checked_assignments;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
