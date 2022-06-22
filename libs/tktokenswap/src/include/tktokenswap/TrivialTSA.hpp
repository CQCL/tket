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

#include "PartialTsaInterface.hpp"

namespace tket {
namespace tsa_internal {

/** A full TSA, simple and fast but not giving very good solutions.
 *  This works by decomposing the desired mapping into abstract disjoint
 *  cycles, decomposing the abstract cycles into lists of abstract swaps,
 *  then finally decomposing the abstract swaps into concrete swaps.
 *  ("Abstract" means that the vertices invloved are not necessarily
 *  adjacent, so the actual swaps cannot be calculated without knowing
 *  the graph, and "concrete" swaps are actual swaps beteen adjacent vertices).
 *  Because the ABSTRACT cycles are disjoint, we are free to perform them,
 *  as long as no other vertices are moved when doing so (they may be moved
 *  in intermediate steps, but will be moved back again by the end of each
 *  cycle). Thus we are guaranteed to get a full solution,
 *  although in tests it can easily give 20-30% more swaps than the best TSA.
 */
class TrivialTSA : public PartialTsaInterface {
 public:
  /** Extra options to control behaviour. */
  enum class Options {
    /** Run the algorithm to completion. */
    FULL_TSA,

    /** Start running the calculated swaps,
     *  but terminate as soon as nonzero L decrease occurs
     *  (which thus gives a Partial TSA).
     */
    BREAK_AFTER_PROGRESS
  };

  /** By default, it's a full TSA.
   *  @param options Option to set behaviour; by default, a full TSA.
   */
  explicit TrivialTSA(Options options = Options::FULL_TSA);

  /** Set another option.
   *  @param options The option to be set from now on.
   */
  void set(Options options);

  /** Calculate and append the complete solution (or break off early if
   *  BREAK_AFTER_PROGRESS was set). The point is that this partial TSA
   *  is not so good, but will be combined with other partial TSAs which
   *  are better, so we want to break off ASAP when progress occurs.
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping, will be updated.
   *  @param distances An object to calculate distances between vertices.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices.
   */
  virtual void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& /*not_needed*/,
      RiverFlowPathFinder& path_finder) override;

  /** The same as the standard append_partial_solution interface,
   *  but without needing to pass in a NeighboursInterface.
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping, will be updated.
   *  @param distances An object to calculate distances between vertices.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices.
   */
  void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, RiverFlowPathFinder& path_finder);

 private:
  //  NOTE: the reason this is all a bit more complicated (and so, the word
  //  "trivial" is a bit unfair) is that we have to allow empty vertices.
  //  With full vertices (every vertex having a token), we can find cycles just
  //  by starting anywhere and going forwards until we hit the start again.
  //  But if some vertices can be empty, we may not be able to go forward
  //  once we hit an empty vertex, so we then have to go backwards also
  //  until we cannot anymore, and finally link the empty end with the nonempty
  //  start vertex to make a cycle.
  //  However, it's really just the same algorithm as the full tokens case.

  Options m_options;

  /** This will contain ALL relevant vertices for ALL cycles, but another
   *  object m_cycle_endpoints will store information about where
   *  each cycle starts and ends.
   */
  VectorListHybrid<size_t> m_abstract_cycles_vertices;
  mutable std::set<size_t> m_vertices_seen;

  typedef VectorListHybrid<size_t>::ID ID;

  /** For an abstract cycle: the first is the ID of the start vertex in
   *  "m_abstract_cycles_vertices" (which already has a builtin linked list
   *  structure), the second is the final vertex.
   */
  typedef std::pair<ID, ID> Endpoints;

  /** Information about where each cycle starts and ends,
   *  using the vertices in m_abstract_cycles_vertices.
   */
  std::vector<Endpoints> m_cycle_endpoints;
  std::vector<size_t> m_vertices_work_vector;

  /** Fills m_abstract_cycles_vertices, m_cycle_endpoints with the cycles.
   *  @param vertex_mapping The current desired mapping.
   */
  void fill_disjoint_abstract_cycles(const VertexMapping& vertex_mapping);

  /** Taking the given first element of "endpoints" as the start vertex,
   *  already known to be in "vertex_mapping", follow the arrows forwards
   *  until no more arrows exist, OR it wraps around to the first vertex,
   *  adding the vertices to "m_abstract_cycles_vertices" as we go,
   *  and updating "endpoints". Does NOT change m_vertices_seen.
   *  @param vertex_mapping The current desired mapping.
   *  @param endpoints The IDs of the vertex endpoints of the desired new cycle
   *    (but only the first ID is valid at the start; the second ID will be
   *    updated).
   *  @return TRUE if a cycle is found, FALSE if it ends at an empty vertex.
   */
  bool grow_cycle_forwards(
      const VertexMapping& vertex_mapping, Endpoints& endpoints);

  /** To be called immediately after grow_cycle_forwards,
   *  if the end vertex did NOT wrap around to the start vertex.
   *  So, go backwards from the start vertex until we cannot any more.
   *  (We can't hit the end vertex since it's empty,
   *  so no arrow can come from there).
   *  Update endpoints.first.
   *  Does NOT change m_vertices_seen. Uses m_reversed_vertex_mapping.
   *  @param endpoints The IDs of the partial vertex cycle start and end
   *    vertices, to be updated (the end of the cycle must wrap round
   *    to the start; the start is not yet determined).
   */
  void grow_cycle_backwards(Endpoints& endpoints);

  /** The ordinary vertex mapping is from v1 to v2,
   *  where v2 is the target of the token currently at v1.
   *  For this mapping, the key is v2, the value is v1.
   */
  VertexMapping m_reversed_vertex_mapping;

  /** Checks validity/consistency of the data in m_abstract_cycles_vertices,
   *  m_cycle_endpoints, m_reversed_vertex_mapping and throws if invalid.
   */
  void do_final_checks() const;

  /** Gets the vertices stored in order in m_abstract_cycles_vertices,
   *  given by the Endpoints, and copies them to m_vertices_work_vector.
   *  (Necessary because we need to do random access, which VectorListHybrid
   *  does not have).
   *  @param endpoints The IDs of the complete vertex cycle start and end,
   *    listed in order in m_abstract_cycles_vertices.
   */
  void copy_vertices_to_work_vector(const Endpoints& endpoints);

  /** Once m_abstract_cycles_vertices and m_cycle_endpoints have been filled,
   *  append the complete solution.
   *  (We don't need to find distances any more, we need actual paths).
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping, will be updated.
   *  @param path_finder The object to calculate a shortest path between any
   *    pair of vertices.
   */
  void append_partial_solution_with_all_cycles(
      SwapList& swaps, VertexMapping& vertex_mapping,
      RiverFlowPathFinder& path_finder);

  /** Perform the single abstract cycle, but breaking off as soon as
   *  the overall total home distance (L) decreases.
   *  (Every abstract cycle has strictly positive L-decrease, otherwise
   *  it wouldn't be included at all, so doing the whole thing must decrease L.
   *  But if we're lucky, we'll decrease L earlier).
   *
   *  Note that we ALSO have to do some estimation, not only to choose
   *  which cycle is likely to be cheap, but ALSO to decide where to
   *  start from. (An ABSTRACT cycle [v0, v1, v2, ..., vn] is decomposed into
   *  ABSTRACT swaps (v0, v1).(v1,v2). ... .(v(n-1), vn), which omits the
   *  abstract swap (vn,v0), but we could have chosen any other v(i) to be
   *  the start vertex. Unlike for CONCRETE swaps, abstract swaps have
   *  different costs, so it's important to choose well).
   *
   *  @param endpoints The IDs of the ends of the final cycle
   *    we've decided to use.
   *  @param start_v_index The starting index in the final cycle vertices,
   *    treating it logically as a vector. (The indices wrap round and reduce
   *    modulo the size).
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping, will be updated.
   *  @param distances An object to calculate distances between vertices.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices.
   *  @return the actual L-decrease (will be strictly positive).
   */
  size_t append_partial_solution_with_single_cycle(
      const Endpoints& endpoints, size_t start_v_index,
      // L (the sum of the distances to home) must decrease
      // by at least this amount, to break off early.
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, RiverFlowPathFinder& path_finder);
};

}  // namespace tsa_internal
}  // namespace tket
