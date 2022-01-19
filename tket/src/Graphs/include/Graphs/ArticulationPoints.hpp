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

#include <boost/graph/biconnected_components.hpp>
#include <set>

#include "Graphs/ArticulationPoints_impl.hpp"

/**
 * Computation of Articulation Points (APs) on undirected graphs
 *
 * This problem is of interest as we can use APs to maintain connectivity
 * requirements in Architectures and QubitGraphs. The BGL already provides an
 * implementation of this, but it does not quite satisfy our needs -- only a
 * loser definition of connectivity matters to us. The following will give more
 * details.
 *
 * Articulation Points (APs) are vertices in a graph that cannot be removed
 * without breaking the graph connectivity.
 * This concept is closely linked to biconnected components, i.e. connected
 * subgraphs in which removing any vertex will not affect the subgraph
 * connectivity.
 *
 * For a given graph G, we can then define a map `belong_to_comp` that
 * maps any vertex to the set of biconnected components it belongs to. Note that
 * any vertex v will belong to at least one biconnected component (which in the
 * worst case will contain a single vertex). It is a well-known graph
 * theoretical result that APs can equivalently be characterised as the vertices
 * that belong to more than one biconnected component. APs (and biconnected
 * components) can be efficiently computed in linear time O(V + E).
 *
 * For our needs, we replace the connectivity requirements in the definition of
 * APs to only consider subgraph connectivity: given a subgraph, we say that the
 * graph is subgraph-connected iff any two vertices of the subgraph are
 * connected in the graph G. I.e. instead of "vanilla" APs, we are only
 * interested in APs that break subgraph connectivity when removed: we call them
 * subgraph APs. Remark that i) the subgraph APs are necessarily contained in
 * the set of all APs of the graph, as our connectivity requirements are
 * strictly weaker. ii) the subgraph APs will not be elements of the subgraph in
 * general
 *
 * BGL implementation:
 * https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/biconnected_components.html
 */

/** These functions only work for the UndirectedConnGraph graph type (typedef in
 * `ArticulationPoints_impl.hpp`) i.e. undirected graphs of DirectedGraph
 * objects They are not defined for BGL graphs in general
 *
 *  The reason for this is that we need the vertices to have a UnitID property
 *  to identify vertices from the main graph and the subgraph
 */

/**
 *  Implementation specific details (in particular the class
 * detail::BicomponentGraph is in the file `ArticulationPoints_impl.hpp`
 */

namespace tket::graphs {

/** Given a graph and a subgraph, returns the subgraph APs, as defined
 *  in the introduction above
 */
template <typename T>
std::set<T> get_subgraph_aps(
    const UndirectedConnGraph<T>& graph,
    const UndirectedConnGraph<T>& subgraph);

// template explicit instations, with implementations in cpp file
extern template std::set<UnitID> get_subgraph_aps(
    const UndirectedConnGraph<UnitID>& graph,
    const UndirectedConnGraph<UnitID>& subgraph);
extern template std::set<Node> get_subgraph_aps(
    const UndirectedConnGraph<Node>& graph,
    const UndirectedConnGraph<Node>& subgraph);

}  // namespace tket::graphs
