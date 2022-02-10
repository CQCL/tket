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
#include <set>
#include <utility>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

/** For a single vertex v in a graph, calculate a filter object
 * based upon the number of vertices at certain distances from v.
 * The objects can then be compared against each other to exclude
 * certain mappings.
 * The trivial theorem here is: if f is a valid monomorphism with
 * f(pv) = tv, and dist(pv,pv2) = d, then dist(f(pv),f(pv2)) <= d.
 * (The target graph may introduce extra edges and thus decrease distances).
 * Thus, we count the number of pv2 at distance d from pv, for d=1,2,3,...
 *
 * It's not necessary to do the whole calculation at once, though;
 * it can be done step-by-step.
 * TODO this is expensive. Make it dynamic, i.e. done throughout the
 * searching on demand, rather than static, i.e. all at initialisation.
 */
class DistanceCounts {
 public:
  /** Perform initial calculations for vertex v.
   * @param neighbours_data an object to extract information about nearby
   * vertices.
   * @param v The vertex attached to this object.
   */
  void initialise(const NeighboursData& neighbours_data, VertexWSM v);

  /** Calculates one more element in "m_counts", the internal data.
   * Returns false if zero has been reached, i.e. the whole graph
   * has been covered, so we should terminate.
   * @param neighbours_data an object to extract information about nearby
   * vertices.
   * @return True if there is still some information to be gained; false if the
   * whole graph has been explored, so we should terminate.
   */
  bool push_back(const NeighboursData& neighbours_data);

  /** Returns the counts.
   * @return A vector x, where x[i] is the number of vertices at distance
   * exactly i+1 from v.
   */
  const std::vector<std::size_t>& get_counts() const;

  /** The size (length) of the counts vector. */
  std::size_t size() const;

  /** Every pattern vertex at distance d can be paired up with a target
   * vertex at distance <= d. Returns TRUE if a matching exists,
   * so it's still possible that the pattern vertex for this object
   * could map to the target vertex used for the other object.
   * Returns FALSE if no matching exists, so it definitely cannot.
   * Notice that it's possible to return FALSE even if the second counts
   * vector has smaller length, so we know that the mapping is impossible
   * WITHOUT having to calculate more target distances.
   * @param other An object calculated for a vertex tv in the target graph.
   * @return False if we are sure that no monomorphism f exists with f(pv)=tv.
   */
  bool test_against_target(const DistanceCounts& other) const;

  /** The function which directly compares the count vectors.
   * @param p_counts The counts, in the pattern graph, for some vertex pv.
   * @param t_counts The counts, in the target graph, for some vertex tv.
   * @return False if the information is sufficient to prove that no
   * monomorphism f exists with f(pv)=tv.
   */
  static bool test_against_target(
      const std::vector<std::size_t>& p_counts,
      const std::vector<std::size_t>& t_counts);

 private:
  /** element [i] gives the number of vertices exactly distance i+1
   * from the initial vertex, which this object was created for.
   */
  std::vector<std::size_t> m_counts;

  std::set<VertexWSM> m_vertices_seen;
  std::set<VertexWSM> m_current_frontier;
  std::set<VertexWSM> m_work_set;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
