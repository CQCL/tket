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

#include "DAGDefs.hpp"

namespace tket {

/**
 * Check that the DAG satisfies the requirements of the Circuit class.
 *
 * These requirements are described below.
 *
 * Definition: call a DAG _balanced_ if every vertex is either an initial vertex
 * (with out-degree 1), a final vertex (with in-degree 1), or an internal vertex
 * with its inbound edges in a defined bijection with its outbound edges. The
 * _balanced degree_ of an internal vertex of a balanced DAG is the common value
 * of its in-degree and out-degree.
 *
 * The edges of G are partitioned into three types: Quantum, Classical, and
 * Boolean.
 *
 * Let G_Q be the subgraph consisting of Quantum edges and their incident
 * vertices. Let G_C be the subgraph consisting of Classical edges and their
 * incident vertices.
 *
 * From these we define the following classes of vertex:
 *
 * - The _Quantum_ vertices: V(G_Q) ∖ V(G_C).
 * - The _Classical_ vertices: V(G_C) ∖ V(G_Q).
 * - The _Measure_ vertices: V(G_Q) ∩ V(G_C).
 *
 * We check the following properties:
 *
 * - V(G) = V(G_Q) ∪ V(G_C).
 * - G_Q and G_C are balanced DAGs, with the bijections defined by the port
 *   numbers on the edges.
 * - Every Measure vertex has balanced degree 1 in G_Q and in G_C.
 * - A Quantum vertex has no Boolean out-edges.
 * - Every source port number on a Boolean edge matches a source port number on
 *   a Classical edge outgoing from the same vertex.
 * - All port numbers on inbound edges to a vertex are distinct.
 *
 * @param      G DAG to check
 *
 * @return whether the DAG has the required properties
 */
bool is_valid(const DAG &G);

}  // namespace tket
