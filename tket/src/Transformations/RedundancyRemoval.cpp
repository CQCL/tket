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
#include <tkassert/Assert.hpp>

#include "BasicOptimisation.hpp"
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

static std::optional<VertexDetachmentInfo> try_detach_identity(
    Circuit &circuit, const Vertex &vertex) {
  auto vertex_operator = circuit.get_Op_ptr_from_Vertex(vertex);
  if (auto phase = vertex_operator->is_identity()) {
    circuit.add_phase(phase.value());
    return std::make_optional(detach_vertex(circuit, vertex));
  }
  return std::nullopt;
}

bool vertex_is_a_measurement(const Circuit &circuit, Vertex const &vertex) {
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

static std::optional<VertexDetachmentInfo> try_detach_zbasis_commuting_vertex(
    Circuit &circuit, const Vertex &vertex) {
  if (vertex_is_succeeded_only_by_z_basis_measurements_with_which_it_commutes(
      circuit, vertex)) {
    return std::make_optional(detach_vertex(circuit, vertex));
  }
  return std::nullopt;
}

bool preliminary_vertex_successor_checks_pass(Circuit &circuit, const Vertex &vertex){

  // check that the classical edges match up correctly
  if (circuit.n_in_edges_of_type(vertex, EdgeType::Boolean) != 0)
    return false;

  // check that both the vertex and its successor have each other and only each
  // other
  auto successors = circuit.get_successors(vertex);
  if (successors.size() != 1) return false;
  auto successor = successors[0];
  if (circuit.get_predecessors(successor).size() != 1) return false;

  // check that successor has adjoint
  if (circuit.get_Op_ptr_from_Vertex(successor)->get_desc().is_oneway()){
    return false;
  }

  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);
  auto vertex_gate = std::dynamic_pointer_cast<const Gate>(vertex_op);

  // check that the ports respect the (a)symmetry between vertices
  for (const Edge &in : circuit.get_in_edges(successor)) {
    auto source_port = circuit.get_source_port(in);
    auto target_port = circuit.get_target_port(in);
    if (not vertex_gate->port_pair_is_symmetric(source_port, target_port)){
      return false;
    }
  }

  return true;
}

std::optional<VertexDetachmentInfo>
try_detach_both_because_successor_is_adjoint(
    Circuit &circuit, Vertex const &vertex, Vertex const &successor) {

  const Op_ptr successor_op = circuit.get_Op_ptr_from_Vertex(successor);
  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);

  if (*vertex_op == *successor_op->dagger())
    return std::make_optional(
        detach_vertex_and_successor(circuit, vertex, successor));

  return std::nullopt;
}

std::optional<VertexDetachmentInfo> try_join_rotations_and_detach_successor(
    Circuit &circuit, Vertex const &vertex, Vertex const &successor) {

  const Op_ptr successor_op = circuit.get_Op_ptr_from_Vertex(successor);
  const OpDesc successor_op_descriptor = successor_op->get_desc();
  const Op_ptr vertex_op = circuit.get_Op_ptr_from_Vertex(vertex);
  const OpDesc vertex_op_descriptor = vertex_op->get_desc();

  // check vertex and successor are same rotation type
  if (not(vertex_op_descriptor.is_rotation() &&
          vertex_op_descriptor.type() == successor_op_descriptor.type()))
    return std::nullopt;

  // replace vertex with combined rotation
  auto expr1 = vertex_op->get_params()[0];
  auto expr2 = successor_op->get_params()[0];
  Op_ptr op_new = get_op_ptr(
      vertex_op_descriptor.type(), {expr1 + expr2},
      circuit.get_in_edges(successor).size());
  circuit.dag[vertex].op = op_new;

  // detach successor only (this adds vertex to list of vertices to check again,
  // so will be removed later if is identity)
  return std::make_optional(detach_vertex(circuit, successor));
}

static std::optional<VertexDetachmentInfo> try_detach_single_vertex(
    Circuit &circuit, const Vertex &vertex) {
  if (auto detachmentInfo = try_detach_identity(circuit, vertex))
    return detachmentInfo;
  if (auto detachmentInfo = try_detach_zbasis_commuting_vertex(circuit, vertex))
    return detachmentInfo;
  return std::nullopt;
}

static std::optional<VertexDetachmentInfo> try_detach_vertex_and_successor(
    Circuit &circuit, const Vertex &vertex) {

  if(not preliminary_vertex_successor_checks_pass(circuit, vertex)) return std::nullopt;
  auto successor = circuit.get_successors(vertex)[0];
  if(auto detachmentInfo = try_detach_both_because_successor_is_adjoint(circuit, vertex, successor)){
    return detachmentInfo;
  }
  if(auto detachmentInfo = try_join_rotations_and_detach_successor(circuit, vertex, successor)){
    return detachmentInfo;
  }
  return std::nullopt;
}

bool is_apriori_not_detachable(const Circuit &circuit, Vertex const &vertex) {
  const OpDesc op_descriptor =
      circuit.get_Op_ptr_from_Vertex(vertex)->get_desc();

  return (not op_descriptor.is_gate()) or  // not a gate
         circuit.n_out_edges(vertex) ==
         0 or  // vertex is boundary or already detached
         circuit.n_in_edges(vertex) == 0;  // vertex is boundary
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

static VertexDetachmentInfo try_detach_vertices(
    Circuit &circuit, const VertexList &vertices) {
  VertexDetachmentInfo detachmentInfo;
  for (const auto &vertex : vertices) {
    detachmentInfo.append(try_detach_vertex(circuit, vertex));
  }
  return detachmentInfo;
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

Transform remove_redundancies() { return Transform(redundancy_removal); }

}  // namespace tket::Transforms
