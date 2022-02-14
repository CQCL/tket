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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include "Utils/SequencedContainers.hpp"
#include "ZX/ZXGenerator.hpp"

/**
 * This header contains the implementation-specific type definitions for the
 * ZXDiagram class. Currently, Boost is used for the underlying graph structure.
 * The interface can be found in ZXDiagram.hpp.
 *
 * A ZX diagram is implemented as a directed graph with labels for its vertices
 * and edges. Although edges in a ZX diagram are ultimately undirected, some
 * vertices are directed or otherwise non-commutative. That is, the permutation
 * in which the edges are connected to the vertex matters, either to
 * distinguish between inputs and outputs of the vertex or operands with
 * different semantics. We refer to such vertices as "directed".
 *
 * To represent these diagrams with undirected edges but some directed vertices
 * we use a directed underlying graph.
 *
 * - The vertex information (e.g. generator type and properties) is represented
 *   by a `ZXGen_ptr` object (equivalent to `Op_ptr` from circuits) containing
 *   the relevant information. In ZX diagrams, the `ZXGenerator` object pointed
 *   to is guaranteed to be exactly one subtype based on its `ZXType`.
 * - For directed vertices, we have a `DirectedZXGenerator` object which
 *   captures the information about the ports. The mapping of ports to semantic
 *   meaning is dictated by the particular `ZXType`. In general, ports are
 *   optional unsigned integers, taking values for directed vertices and
 *   `std::nullopt` for undirected.
 * - The edges store information including its type (ZXWireType is either Basic
 *   or H), its quantum-ness (QuantumType is either Quantum or Classical), and
 *   the ports it connects to on both ends (as a map from WireEnd::Source/
 *   Target to a the port it connects to on a directed vertex). This
 *   differentiation between source and target vertices is to allow a unique
 *   representation of the edge's connectivity in case of connecting between
 *   two directed vertices.
 */

namespace tket {

namespace zx {

// ZXVert information - each vertex just captures a ZX generator
struct ZXVertProperties {
  ZXGen_ptr op;
};

/**
 * Wire properties:
 * `type`: ZXWireType::Basic or H, whether the wire is an identity / Hadamard
 * `qtype`: QuantumType::Quantum or Classical, whether the wire is doubled or
 * not under the CPM construction `source_port`, `target_port`: the ports the
 * wire connects to on the source and target vertices if they are directed, and
 * `nullopt` if they are undirected
 */
struct WireProperties {
  ZXWireType type;
  QuantumType qtype;
  std::optional<unsigned> source_port;
  std::optional<unsigned> target_port;

  WireProperties();
  WireProperties(
      ZXWireType type, QuantumType qtype,
      std::optional<unsigned> source_port = std::nullopt,
      std::optional<unsigned> target_port = std::nullopt);
  // Use default copy and move constructors

  bool operator==(const WireProperties& other) const;
};

/**
 * A ZX diagram is semantically undirected, but the implementation uses a
 * directed graph in order to represent directed and non-commutative vertices.
 *
 * `listS` is used for both the edge and vertex lists as these maintain
 * descriptor / pointer stability:
 * https://www.boost.org/doc/libs/1_54_0/libs/graph/doc/adjacency_list.html
 *
 * See TKET Wiki entry: "Boost & igraph C++ Graphical Libraries Benchmarks"
 * for earlier BGL tests.
 */
typedef boost::adjacency_list<
    boost::listS, boost::listS, boost::bidirectionalS, ZXVertProperties,
    WireProperties>
    ZXGraph;

typedef boost::graph_traits<ZXGraph>::vertex_descriptor ZXVert;
typedef std::vector<ZXVert> ZXVertVec;
typedef sequence_set_t<ZXVert> ZXVertSeqSet;

typedef boost::graph_traits<ZXGraph>::edge_descriptor Wire;
typedef std::vector<Wire> WireVec;

// (convenience) vertex and edge iterators
typedef boost::graph_traits<ZXGraph>::vertex_iterator ZXVertIterator;
// Iterator over vertex neighbourhoods
typedef boost::graph_traits<ZXGraph>::adjacency_iterator NeighbourIterator;
// Iterator over *undirected* edges
typedef boost::graph_traits<ZXGraph>::edge_iterator WireIterator;
// Iterator over *directed* edges
typedef boost::graph_traits<ZXGraph>::out_edge_iterator OutWireIterator;
typedef boost::graph_traits<ZXGraph>::in_edge_iterator InWireIterator;

}  // namespace zx

}  // namespace tket
