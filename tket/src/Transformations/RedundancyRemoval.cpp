// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "BasicOptimisation.hpp"

#include <optional>
#include <tkassert/Assert.hpp>

#include "Characterisation/ErrorTypes.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Gate/Gate.hpp"
#include "Transform.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket::Transforms {



// A helper struct for holding which vertices have been detached (a.k.a bin)
// the predecessors of those vertices
struct VertexDetachmentInfo {
  VertexList detachedVertices;
  VertexList detachedVertexPredecessors;

  // Merge vertices from temporary VertexDetachmentInfo object to this one
  void append(VertexDetachmentInfo &&detachmentInfoToAppend) {
    detachedVertices.merge(detachmentInfoToAppend.detachedVertices);
    detachedVertexPredecessors.merge(
        detachmentInfoToAppend.detachedVertexPredecessors);
  }
  // Returns VertexDetachmentInfo object with empty lists
  static VertexDetachmentInfo Empty() { return {}; }
};

static VertexList detach_redundant_vertices(Circuit &circuit);
static VertexDetachmentInfo try_detach_vertices(
    Circuit &circuit, const VertexList &vertices);
static VertexDetachmentInfo try_detach_vertex(
    Circuit &circuit, const Vertex &vertex);
static std::optional<VertexDetachmentInfo> try_detach_single_vertex(
    Circuit &circuit, const Vertex &vertex);
static std::optional<VertexDetachmentInfo> try_detach_identity(
    Circuit &circuit, const Vertex &vertex);
static std::optional<VertexDetachmentInfo> try_detach_zbasis_commuting_vertex(
    Circuit &circuit, const Vertex &vertex);
static bool
vertex_is_succeeded_only_by_z_basis_measurements_with_which_it_commutes(
    const Circuit &circuit, const Vertex &vertex);
static VertexDetachmentInfo detach_vertex(
    Circuit &circuit, const Vertex &vertex);
static std::optional<VertexDetachmentInfo> try_detach_vertex_and_successor(
    Circuit &circuit, const Vertex &vertex);

Transform remove_redundancies(){return Transform(redundancy_removal);}
// this method annihilates all primitives next to each other (accounting for
// previous annihilation)
// also removes redundant non-classically controlled Z basis gates before a z
// basis measurement so that eg. -H-X-X-H- always annihilates to -----
bool redundancy_removal(Circuit &circuit) {
  VertexList bin = detach_redundant_vertices(circuit);
  circuit.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return !bin.empty();
}

static VertexList detach_redundant_vertices(Circuit &circuit) {
  VertexList bin;
  VertexList verticesToCheckForRemoval = circuit.vertices_list_in_order();
  while (!verticesToCheckForRemoval.empty()) {
    auto detachmentInfo =
        try_detach_vertices(circuit, verticesToCheckForRemoval);
    bin.merge(detachmentInfo.detachedVertices);
    std::swap(
        verticesToCheckForRemoval, detachmentInfo.detachedVertexPredecessors);
  }
  return bin;
}

static VertexDetachmentInfo try_detach_vertices(
    Circuit &circuit, const VertexList &vertices) {
  VertexDetachmentInfo detachmentInfo;
  for (const auto &vertex : vertices) {
    detachmentInfo.append(try_detach_vertex(circuit, vertex));
  }
  return detachmentInfo;
}

static bool is_apriori_not_detachable(
    const Circuit &circuit, const Vertex &vertex) {
  const OpDesc op_descriptor =
      circuit.get_Op_ptr_from_Vertex(vertex)->get_desc();

  return (not op_descriptor.is_gate()) or
         (circuit.n_out_edges(vertex) == 0 or
          circuit.n_in_edges(vertex) ==
              0);  // vertex is boundary or already detached
}

static VertexDetachmentInfo try_detach_vertex(
    Circuit &circuit, const Vertex &vertex) {
  if (is_apriori_not_detachable(circuit, vertex))
    return VertexDetachmentInfo::Empty();
  if (auto detachmentInfo = try_detach_single_vertex(circuit, vertex))
    return detachmentInfo.value();
  if (auto detachmentInfo = try_detach_vertex_and_successor(circuit, vertex))
    return detachmentInfo.value();

  return VertexDetachmentInfo::Empty();
}

static std::optional<VertexDetachmentInfo> try_detach_single_vertex(
    Circuit &circuit, const Vertex &vertex) {
  if (auto detachmentInfo = try_detach_identity(circuit, vertex))
    return detachmentInfo;
  if (auto detachmentInfo = try_detach_zbasis_commuting_vertex(circuit, vertex))
    return detachmentInfo;
  return std::nullopt;
}

static std::optional<VertexDetachmentInfo> try_detach_identity(
    Circuit &circuit, const Vertex &vertex) {
  auto vertex_operator = circuit.get_Op_ptr_from_Vertex(vertex);
  if (auto phase = vertex_operator->is_identity()) {
    circuit.add_phase(phase.value());
    return std::make_optional(detach_vertex(circuit, vertex));
  }
  return std::nullopt;
}

static std::optional<VertexDetachmentInfo> try_detach_zbasis_commuting_vertex(
    Circuit &circuit, const Vertex &vertex) {
  if (vertex_is_succeeded_only_by_z_basis_measurements_with_which_it_commutes(
          circuit, vertex)) {
    return std::make_optional(detach_vertex(circuit, vertex));
  }
  return std::nullopt;
}

static VertexDetachmentInfo detach_vertex(
    Circuit &circuit, const Vertex &vertex) {
  auto detachInfo =
      VertexDetachmentInfo{{vertex}, circuit.get_predecessors_list(vertex)};
  circuit.remove_vertex(
      vertex, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
  return detachInfo;
}

static VertexDetachmentInfo detach_vertex_and_successor(
    Circuit &circuit, const Vertex &vertex, const Vertex &successor) {
  auto detachInfo = VertexDetachmentInfo{
      {vertex, successor}, circuit.get_predecessors_list(vertex)};
  circuit.remove_vertices(
      VertexList{vertex, successor}, Circuit::GraphRewiring::Yes,
      Circuit::VertexDeletion::No);
  return detachInfo;
}

static bool vertex_is_a_measurement(
    const Circuit &circuit, const Vertex &vertex) {
  return circuit.get_OpType_from_Vertex(vertex) == OpType::Measure;
}

static bool
vertex_is_succeeded_only_by_z_basis_measurements_with_which_it_commutes(
    const Circuit &circuit, const Vertex &vertex) {
  // If vertex has no classical out edges, no need to continue
  if (circuit.n_out_edges_of_type(vertex, EdgeType::Classical) != 0)
    return false;
  auto successors = circuit.get_successors(vertex);
  std::vector<port_t> ports(successors.size());
  std::iota(ports.begin(), ports.end(), 0);

  return std::all_of(ports.cbegin(), ports.cend(), [&](port_t port) {
    return vertex_is_a_measurement(circuit, successors[port]) &&
           circuit.commutes_with_basis(
               vertex, Pauli::Z, PortType::Source, port);
  });
}

static std::optional<VertexDetachmentInfo> try_detach_vertex_and_successor(
    Circuit &circuit, const Vertex &vertex) {
  auto successors = circuit.get_successors(vertex);

  if (successors.size() != 1) return std::nullopt;
  if (circuit.get_predecessors(successors[0]).size() != 1) return std::nullopt;

  auto successor = successors[0];

  for (const Edge &in : circuit.get_in_edges(successor)) {
    if (circuit.get_source_port(in) != circuit.get_target_port(in))
      return std::nullopt;
  }

  if (circuit.n_in_edges_of_type(vertex, EdgeType::Boolean) != 0)
    return std::nullopt;

  const Op_ptr successor_op = circuit.get_Op_ptr_from_Vertex(successor);
  const OpDesc successor_op_descriptor = successor_op->get_desc();

  if (successor_op_descriptor.is_oneway()) return std::nullopt;

  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);

  if (*vertex_op == *successor_op->dagger())
    return std::make_optional(
        detach_vertex_and_successor(circuit, vertex, successor));

  const OpDesc vertex_op_descriptor = vertex_op->get_desc();

  if (not(vertex_op_descriptor.is_rotation() &&
          vertex_op_descriptor.type() == successor_op_descriptor.type()))
    return std::nullopt;

  auto expr1 = vertex_op->get_params()[0];
  auto expr2 = successor_op->get_params()[0];
  Op_ptr op_new = get_op_ptr(
      vertex_op_descriptor.type(), {expr1 + expr2},
      circuit.get_in_edges(successor).size());
  circuit.dag[vertex].op = op_new;
  return std::make_optional(detach_vertex(circuit, successor));
}

}  // namespace tket::Transforms
