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
#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace tket {
namespace WeightedSubgraphMonomorphism {

typedef std::uint_fast32_t VertexWSM;
typedef std::uint64_t WeightWSM;

/** (v1,v2) with v1<v2. */
typedef std::pair<VertexWSM, VertexWSM> EdgeWSM;

/** KEY: Edge (v1,v2) with v1<v2,
 * VALUE: the   weight >= 0  of the edge v1-v2.
 */
typedef std::map<EdgeWSM, WeightWSM> GraphEdgeWeights;

/** KEY: pattern vertex
// VALUE: all possible target vertices it could map to.
 */
typedef std::map<VertexWSM, std::set<VertexWSM>> PossibleAssignments;

/** When you're SURE that Dom(pv) = { tv } has size 1, return tv. */
VertexWSM get_single_value(
    VertexWSM pv, const PossibleAssignments& assignments);

/** A string representation of a graph with edge weights. */
std::string str(const GraphEdgeWeights& gdata);

/** A string representation of a list of edges. */
std::string str(const std::vector<EdgeWSM>&);

/** Ensures consistency. Checks that v1 != v2.
 * @param v1 First vertex.
 * @param v2 Second vertex.
 * @return An edge (a,b), with a<b always.
 */
EdgeWSM get_edge(VertexWSM v1, VertexWSM v2);

/** Get the nonisolated vertices (i.e., vertices appearing in the edges).
 * @param edges_and_weights The graph (with edge weights, which are ignored).
 * @return A list of all nonisolated vertices (i.e., vertices appearing in the
 * edges), unique and sorted.
 */
std::vector<VertexWSM> get_vertices(const GraphEdgeWeights& edges_and_weights);

/** For use in reducing a search node, or individual Domain(pv). */
enum class ReductionResult {
  SUCCESS,

  /** At least one new assignment PV->TV was created;
   * a special case needing different treatment to "SUCCESS".
   */
  NEW_ASSIGNMENTS,

  /** The node is impossible; some domain has become empty. */
  NOGOOD
};

/** KEY: pattern vertex
 * VALUE: target vertex
 */
typedef std::map<VertexWSM, VertexWSM> Assignments;

/** A string representation of some assignments. */
std::string str(const Assignments&);

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
