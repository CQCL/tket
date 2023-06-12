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

#include <optional>

#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Gate/Gate.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket::Transforms {

// A helper struct for holding which vertices have been detached (a.k.a bin)
// the predecessors of those vertices
struct VertexDetachmentInfo {
  VertexVec detachedVertices;
  VertexVec detachedVertexPredecessors;
};

static void detach_vertex(
    Circuit &circuit, const Vertex &vertex,
    VertexDetachmentInfo &detachmentInfo) {
  detachmentInfo.detachedVertices.emplace_back(vertex);
  auto &vertexPredecessors = detachmentInfo.detachedVertexPredecessors;
  const auto &newVertexPredecessors = circuit.get_predecessors(vertex);
  vertexPredecessors.insert(
      vertexPredecessors.cend(), newVertexPredecessors.cbegin(),
      newVertexPredecessors.cend());
  circuit.remove_vertex(
      vertex, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
}

static void detach_vertex_and_successor(
    Circuit &circuit, const Vertex &vertex, const Vertex &successor,
    VertexDetachmentInfo &detachmentInfo) {
  detachmentInfo.detachedVertices.emplace_back(vertex);
  detachmentInfo.detachedVertices.emplace_back(successor);
  auto &vertexPredecessors = detachmentInfo.detachedVertexPredecessors;
  const auto &newVertexPredecessors = circuit.get_predecessors(vertex);
  vertexPredecessors.insert(
      vertexPredecessors.cend(), newVertexPredecessors.cbegin(),
      newVertexPredecessors.cend());
  circuit.remove_vertices(
      VertexList{vertex, successor}, Circuit::GraphRewiring::Yes,
      Circuit::VertexDeletion::No);
}

static bool try_detach_identity(
    Circuit &circuit, const Vertex &vertex,
    VertexDetachmentInfo &detachmentInfo) {
  auto vertex_operator = circuit.get_Op_ptr_from_Vertex(vertex);
  if (auto phase = vertex_operator->is_identity()) {
    circuit.add_phase(phase.value());
    detach_vertex(circuit, vertex, detachmentInfo);
    return true;
  }
  return false;
}

static bool vertex_is_a_measurement(
    const Circuit &circuit, Vertex const &vertex) {
  return circuit.get_OpType_from_Vertex(vertex) == OpType::Measure;
}

static bool
vertex_is_succeeded_only_by_z_basis_measurements_with_which_it_commutes(
    const Circuit &circuit, const Vertex &vertex) {
  // If vertex has classical out edges, no need to continue
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

static bool try_detach_zbasis_commuting_vertex(
    Circuit &circuit, const Vertex &vertex,
    VertexDetachmentInfo &detachmentInfo) {
  if (vertex_is_succeeded_only_by_z_basis_measurements_with_which_it_commutes(
          circuit, vertex)) {
    detach_vertex(circuit, vertex, detachmentInfo);
    return true;
  }
  return false;
}

static bool port_ordering_is_compatible(
    Circuit &circuit, const Vertex &vertex, const Vertex &successor) {
  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);
  // Vertex port must be symmetrically equivalent to successor port for any edge
  // between them Examples:
  //  CX[0,1]==CX[0,1] passes test because ports line up
  //  CX[0,1]==CX[1,0] fails test because ports don't line up and CX is not
  //  invariant with respect to exchange of qubits CZ[0,1]==CZ[1,0] passes test
  //  because CZ is invariant with respect to exchange of qubits
  for (const Edge &in : circuit.get_in_edges(successor)) {
    auto source_port = circuit.get_source_port(in);
    auto target_port = circuit.get_target_port(in);
    if (not vertex_op->has_symmetry(source_port, target_port)) {
      return false;
    }
  }
  return true;
}

static bool preliminary_vertex_successor_checks_pass(
    Circuit &circuit, const Vertex &vertex) {
  // check that the classical edges match up correctly
  if (circuit.n_in_edges_of_type(vertex, EdgeType::Boolean) != 0) return false;

  // check that both the vertex and its successor have each other and only each
  // other
  auto successors = circuit.get_successors(vertex);
  if (successors.size() != 1) return false;
  auto successor = successors[0];
  if (circuit.get_predecessors(successor).size() != 1) return false;

  // check that successor has adjoint
  if (circuit.get_Op_ptr_from_Vertex(successor)->get_desc().is_oneway()) {
    return false;
  }

  return port_ordering_is_compatible(circuit, vertex, successor);
}

static bool try_detach_both_because_successor_is_adjoint(
    Circuit &circuit, Vertex const &vertex, Vertex const &successor,
    VertexDetachmentInfo &detachmentInfo) {
  const Op_ptr successor_op = circuit.get_Op_ptr_from_Vertex(successor);
  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);
  if (*vertex_op == *successor_op->dagger()) {
    detach_vertex_and_successor(circuit, vertex, successor, detachmentInfo);
    return true;
  }
  return false;
}

static bool try_join_rotations_and_detach_successor(
    Circuit &circuit, Vertex const &vertex, Vertex const &successor,
    VertexDetachmentInfo &detachmentInfo) {
  const Op_ptr successor_op = circuit.get_Op_ptr_from_Vertex(successor);
  const OpDesc successor_op_descriptor = successor_op->get_desc();
  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);
  const OpDesc vertex_op_descriptor = vertex_op->get_desc();

  // check vertex and successor are same rotation type
  if (not(vertex_op_descriptor.is_rotation() &&
          vertex_op_descriptor.type() == successor_op_descriptor.type())) {
    return false;
  }

  // replace vertex with combined rotation
  auto expr1 = vertex_op->get_params()[0];
  auto expr2 = successor_op->get_params()[0];
  Op_ptr op_new = get_op_ptr(
      vertex_op_descriptor.type(), {expr1 + expr2},
      circuit.get_in_edges(successor).size());
  circuit.dag[vertex].op = op_new;

  // detach successor only (this adds vertex to list of vertices to check again,
  // so will be removed later if is identity)
  detach_vertex(circuit, successor, detachmentInfo);
  return true;
}

static bool try_detach_single_vertex(
    Circuit &circuit, const Vertex &vertex,
    VertexDetachmentInfo &detachmentInfo) {
  if (try_detach_identity(circuit, vertex, detachmentInfo)) {
    return true;
  }
  if (try_detach_zbasis_commuting_vertex(circuit, vertex, detachmentInfo)) {
    return true;
  }
  return false;
}

static bool try_detach_vertex_and_successor(
    Circuit &circuit, const Vertex &vertex,
    VertexDetachmentInfo &detachmentInfo) {
  if (not preliminary_vertex_successor_checks_pass(circuit, vertex))
    return false;
  auto successor = circuit.get_successors(vertex)[0];
  if (try_detach_both_because_successor_is_adjoint(
          circuit, vertex, successor, detachmentInfo)) {
    return true;
  }
  if (try_join_rotations_and_detach_successor(
          circuit, vertex, successor, detachmentInfo)) {
    return true;
  }
  return false;
}

static bool is_apriori_not_detachable(
    const Circuit &circuit, Vertex const &vertex) {
  const OpDesc op_descriptor =
      circuit.get_Op_ptr_from_Vertex(vertex)->get_desc();

  return (not op_descriptor.is_gate()) or  // not a gate
         circuit.n_out_edges(vertex) ==
             0 or  // vertex is boundary or already detached
         circuit.n_in_edges(vertex) == 0;  // vertex is boundary
}

static bool try_detach_vertex(
    Circuit &circuit, const Vertex &vertex,
    VertexDetachmentInfo &detachmentInfo) {
  if (is_apriori_not_detachable(circuit, vertex)) {
    return false;
  }
  if (try_detach_single_vertex(circuit, vertex, detachmentInfo)) {
    return true;
  }
  if (try_detach_vertex_and_successor(circuit, vertex, detachmentInfo)) {
    return true;
  }
  return false;
}

static VertexDetachmentInfo detach_vertices_if_redundant(
    Circuit &circuit, const VertexVec &vertices) {
  VertexDetachmentInfo detachmentInfo;
  for (const auto &vertex : vertices) {
    try_detach_vertex(circuit, vertex, detachmentInfo);
  }
  return detachmentInfo;
}

static VertexSet detach_any_redundant_vertices(Circuit &circuit) {
  VertexSet bin;
  VertexVec verticesToCheckForRemoval = circuit.vertices_in_order();
  while (!verticesToCheckForRemoval.empty()) {
    auto detachmentInfo =
        detach_vertices_if_redundant(circuit, verticesToCheckForRemoval);
    bin.insert(
        detachmentInfo.detachedVertices.cbegin(),
        detachmentInfo.detachedVertices.cend());
    std::swap(
        verticesToCheckForRemoval, detachmentInfo.detachedVertexPredecessors);
  }
  return bin;
}

// this method annihilates all primitives next to each other (accounting for
// previous annihilation)
// also removes redundant non-classically controlled Z basis gates before a z
// basis measurement so that eg. -H-X-X-H- always annihilates to -----
bool redundancy_removal(Circuit &circuit) {
  VertexSet bin = detach_any_redundant_vertices(circuit);
  circuit.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return !bin.empty();
}

Transform remove_redundancies() { return Transform(redundancy_removal); }

}  // namespace tket::Transforms
