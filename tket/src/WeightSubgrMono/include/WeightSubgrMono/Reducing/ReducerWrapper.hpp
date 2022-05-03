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

class DomainsAccessor;

class ReducerInterface {
 public:
  virtual void clear();

  /** Check if pv->tv may be valid, considered in isolation from all other
   * assignments. This should be cheaper than a reduction. By default, just
   * returns true always.
   */
  virtual bool check(std::pair<VertexWSM, VertexWSM> assignment);

  /** Given that PV->TV is a new assignment, reduces the domains
   * of all affected vertices. Breaks off early
   * if new assignments arise in the node due to reductions.
   * Of course, a nogood here does NOT mean that PV->TV is invalid always;
   * just that it is invalid IN COMBINATION with the complete collection
   * of domains.
   */
  virtual ReductionResult reduce(
      std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
      std::set<VertexWSM>& work_set);

  virtual ~ReducerInterface() = default;

 protected:
  /*
  Most natural reducers have a symmetry. Suppose that

  pv1 -> tv1   gives a map (function)
          M(reducer,pv1,tv1) : {all PV} -> PowerSet({all TV}),

  meaning (by definition) that: for every pv2,

  Domain(pv2) is a subset of S = M(reducer,pv1,tv1)[pv2]

  (i.e., as soon as we know that the assignment  pv1 -> tv1
  has been made, we can intersect the current Domain(pv2) with S).

  [We can of course let S be the set of all TV for some pv2,
  in which case this reducer has no effect upon Dom(pv2),
  i.e. pv1, pv2 do not affect each other, as far as this reducer knows].

  THEN, it is often the case that if subsequently we reduce further
  (in any way - perhaps from other unrelated reducers),
  and obtain a new assignment  pv2 -> tv2,
  (so that tv2 is in M(reducer,pv1,tv1)[pv2]),
  then our reducer is automatically consistent with pv1;
  we do not need to check if  pv1 -> tv1  is valid (for this reducer).

  Another way to phrase this: the reducer satisfies

  tv2 in M(reducer,pv1,tv1)[pv2]   <==>
  tv1 in M(reducer,pv2,tv2)[pv1].

  E.g., if "reducer" is the distance reducer with parameter d, then:

  M(., pv1, tv1) maps every  pv2  with dist(pv1, pv2) = d
  to the set  { tv2 : dist(tv1, tv2) <= d }.

  Thus by definition:

  tv2 in M(.,pv1,tv1)[pv2]   <==>   dist(pv1,pv2)=d,  dist(tv1,tv2)<=d.

  It's obvious from this that we can swap pv1,pv2 and tv1,tv2.
  */

  /** For use within "reduce", for reducers possessing this symmetry
   * (otherwise this function is meaningless).
   * Check if the reduction can be skipped for the other pv to save time.
   * @param other_domain Domain(pv2). We will check if Domain(pv2) = {tv2}.
   * @param accessor An object to retrieve data about domains.
   * @param this_vertex Vertex pv1, where we are checking the assignment
   * pv1->tv1.
   * @param other_vertex Vertex pv2, which may already have had Domain(pv2)
   * reduced to a singleton.
   * @return True if there is no need to reduce Domain(pv2).
   */
  bool other_vertex_reduction_can_be_skipped_by_symmetry(
      const std::set<VertexWSM>& other_domain, const DomainsAccessor& accessor,
      VertexWSM this_vertex, VertexWSM other_vertex);
};

/** A class to wrap a raw Reducer object, and keep track of which assignments
 * have already been processed.
 */
class ReducerWrapper {
 public:
  explicit ReducerWrapper(ReducerInterface& reducer_interface);

  /** Call at the start, when we are about to begin reducing a node. */
  void clear();

  /** Checks if the given PV->TV assignment appears to be valid,
   * separately from others (i.e., in isolation).
   * Does NOT keep track of whether it was checked before.
   */
  bool check(std::pair<VertexWSM, VertexWSM> assignment);

  /** Keeps track of previously processed assignments and doesn't repeat them.
   * @param accessor An object to get information about the current domains, and
   * alter them if necessary.
   * @param work_set An object for reuse, to avoid expensive memory
   * reallocations.
   * @return The result of reducing the current node with ALL new assignments
   * not yet processed. Can be resumed later if breaking off early due to new
   * assignments occurring.
   */
  ReductionResult reduce(
      DomainsAccessor& accessor, std::set<VertexWSM>& work_set);

 private:
  ReducerInterface& m_reducer;
  std::size_t m_number_of_processed_assignments;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
