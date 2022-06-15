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

#include <list>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "OpType/EdgeType.hpp"
#include "Ops/Op.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

/** Description of a node in a circuit, representing some operation */
struct VertexProperties {
  Op_ptr op;                          /**< operation */
  std::optional<std::string> opgroup; /**< operation group identfier */

  VertexProperties(
      Op_ptr op = 0, std::optional<std::string> opgroup = std::nullopt)
      : op(op), opgroup(opgroup) {}
};

/** A specific entry or exit port of a vertex */
typedef unsigned port_t;

/** Whether a vertex port is out-going (source) or in-coming (target) */
enum class PortType { Source, Target };

/** Description of an edge in a circuit, representing a directional wire */
struct EdgeProperties {
  EdgeType type; /**< type of wire */
  std::pair<port_t, port_t> ports;
};

/** Graph representing a circuit, with operations as nodes. */
typedef boost::adjacency_list<
    // OutEdgeList
    boost::listS,

    // VertexList (use listS because we want to be able to remove vertices
    // without invalidating iterators)
    boost::listS,

    // we want access to incoming and outgoing edges
    boost::bidirectionalS,

    // indexing needed for algorithms such as topological sort
    boost::property<boost::vertex_index_t, int, VertexProperties>,

    EdgeProperties>
    DAG;

typedef boost::graph_traits<DAG>::vertex_descriptor Vertex;
typedef boost::graph_traits<DAG>::vertex_iterator V_iterator;
typedef std::unordered_set<Vertex> VertexSet;
typedef std::vector<Vertex> VertexVec;
typedef std::list<Vertex> VertexList;
typedef std::unordered_map<Vertex, unsigned> IndexMap;
typedef boost::adj_list_vertex_property_map<
    DAG, int, int&, boost::vertex_index_t>
    VIndex;

/**
 * A vertex with an index.
 *
 * This can be used instead of a plain @ref Vertex in associative containers
 * where control over the order of iteration is required.
 */
typedef std::pair<unsigned, Vertex> IVertex;

typedef boost::graph_traits<DAG>::edge_descriptor Edge;
typedef boost::graph_traits<DAG>::edge_iterator E_iterator;
typedef DAG::in_edge_iterator E_in_iterator;
typedef DAG::out_edge_iterator E_out_iterator;
typedef std::set<Edge> EdgeSet;
typedef std::vector<Edge> EdgeVec;
typedef std::list<Edge> EdgeList;

typedef std::pair<Vertex, port_t> VertPort;

}  // namespace tket
