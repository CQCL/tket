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
#include <boost/dynamic_bitset.hpp>
#include <optional>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;

/** A general interface for an object which knows, using various
 * graph-theoretic tools, that if a new assignment PV->TV is made,
 * then certain other domains can be reduced,
 * thus possibly pruning the search space and speeding up the search.
 * However, these objects do NOT remember which assignments were previously
 * processed.
 *
 * The reason for having a "check" and "reduce" function is that, in most
 * natural cases, the SAME graph theory which enables a particular reduction
 * ALSO can be applied to give a useful "check" function.
 */
class ReducerInterface {
 public:
  /** Check if pv->tv may be valid, considered in isolation from all other
   * assignments. This should be cheaper than a reduction. By default, just
   * returns true always.
   * Of course it is a waste of time to call this multiple times with the same
   * assignment; the CALLER is responsible for arranging things so as to avoid
   * duplicated calls.
   * Note that this doesn't actually alter any domains (it doesn't even know
   * about the current domains). The caller should erase PV->TV completely
   * from ALL data whenever this returns false.
   * @param assignment An assignment PV->TV under consideration.
   * @return False if the assignment is and was ALWAYS invalid (NOT just because
   * of the other domains at this exact moment in the search). Thus, TV can be
   * erased from EVERY Domain(PV) in EVERY node of EVERY search, not just the
   * current node).
   */
  virtual bool check(std::pair<VertexWSM, VertexWSM> assignment) = 0;

  /** Given that PV->TV is a new assignment, reduces the domains
   * of all affected vertices. Breaks off early
   * if new assignments arise in the node due to reductions.
   * Of course, UNLIKE when the "check" function returns false, a nogood here
   * does NOT mean that PV->TV is invalid always;
   * just that it is invalid IN COMBINATION with the complete collection
   * of domains in the CURRENT node.
   * @param assignment An assignment PV->TV under consideration.
   * @param accessor An object providing read/write access to the complete
   * collection of domains in the CURRENT node.
   * @param work_set A reusable object which may or may not be used by this
   * class, to avoid unnecessary memory reallocation.
   * @return After reducing the domains, information about what happened; did we
   * hit a nogood? Did we create a new assignment?
   */
  virtual ReductionResult reduce(
      std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
      boost::dynamic_bitset<>& work_set) = 0;

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
  static bool other_vertex_reduction_can_be_skipped_by_symmetry(
      const boost::dynamic_bitset<>& other_domain,
      const DomainsAccessor& accessor, VertexWSM this_vertex,
      VertexWSM other_vertex);
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
   * For convenience, just wrap around the "check" function of the
   * ReducerInterface object.
   * @param assignment An assignment PV->TV under consideration.
   * @return False if the assignment is and was ALWAYS invalid, in ALL search
   * nodes.
   */
  bool check(std::pair<VertexWSM, VertexWSM> assignment);

  /** Starts to reduce domains using any new assignments which arose since
   * the last call.
   * Keeps track of previously processed assignments and doesn't repeat them.
   * Breaks off early if a new assignment is created, but will subsequently
   * resume automatically from the next unprocessed assignment if this occurs.
   * The reason is that new assignments have cascading effects, often reducing
   * many other domains simultaneously and cheaply. Thus it is better to let
   * the caller perform these cheap reductions first, before carrying out the
   * reductions in this object (which are usually much more expensive).
   * @param accessor An object to get information about the current domains, and
   * alter them if necessary. Also, gets the list of new assignments from this
   * object.
   * @param work_set An object for reuse, to avoid expensive memory
   * reallocations.
   * @return The result of reducing the current node with ALL new assignments
   * not yet processed. Can be resumed later if breaking off early due to new
   * assignments occurring.
   */
  ReductionResult reduce(
      DomainsAccessor& accessor, boost::dynamic_bitset<>& work_set);

 private:
  ReducerInterface& m_reducer;
  std::size_t m_number_of_processed_assignments;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
