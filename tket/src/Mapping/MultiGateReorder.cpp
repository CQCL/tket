// Copyright Quantinuum
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

#include "tket/Mapping/MultiGateReorder.hpp"

#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Mapping/MappingFrontier.hpp"

namespace tket {

MultiGateReorder::MultiGateReorder(
    const ArchitecturePtr &_architecture,
    MappingFrontier_ptr &_mapping_frontier)
    : architecture_(_architecture), mapping_frontier_(_mapping_frontier) {
  // This needs to be updated every time the frontier changes
  this->u_frontier_edges_ =
      convert_u_frontier_to_edges(*frontier_convert_vertport_to_edge(
          _mapping_frontier->circuit_, _mapping_frontier->linear_boundary));
}

// Traverse the DAG to the quantum frontier
// to find the UnitID associated with an VertPort
static UnitID get_unitid_from_vertex_port(
    const MappingFrontier_ptr &frontier, const VertPort &vert_port) {
  VertPort current_vert_port = vert_port;
  while (true) {
    auto it =
        frontier->linear_boundary->get<TagValue>().find(current_vert_port);
    if (it != frontier->linear_boundary->get<TagValue>().end()) {
      return it->first;
    }
    Edge current_e = frontier->circuit_.get_nth_out_edge(
        current_vert_port.first, current_vert_port.second);
    Vertex prev_vert;
    Edge prev_e;
    std::tie(prev_vert, prev_e) =
        frontier->circuit_.get_prev_pair(current_vert_port.first, current_e);
    current_vert_port = {prev_vert, frontier->circuit_.get_source_port(prev_e)};
  }
}

static bool is_multiq_quantum_gate(const Circuit &circ, const Vertex &vert) {
  Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
  return (
      op->get_desc().is_gate() && circ.n_in_edges(vert) > 1 &&
      circ.n_in_edges_of_type(vert, EdgeType::Quantum) ==
          circ.n_in_edges(vert) &&
      circ.n_out_edges_of_type(vert, EdgeType::Quantum) ==
          circ.n_out_edges(vert));
}

// This method will try to commute a vertex to the quantum frontier
// If successful, returns the current in-edges and the target in-edges to
// rewire the vertex
static std::optional<std::pair<EdgeVec, EdgeVec>> try_find_commute_edges(
    const Circuit &circ, const EdgeVec &frontier_edges, const Vertex &vert) {
  // Initialize to be the in_edges for the given vertex
  EdgeVec current_edges = circ.get_in_edges(vert);
  EdgeVec initial_edges(current_edges.begin(), current_edges.end());

  // Record the colour of each port of the vertex.
  std::vector<std::optional<Pauli>> colours;
  for (const Edge &edge : current_edges) {
    port_t target_port = circ.get_target_port(edge);
    std::optional<Pauli> colour = circ.commuting_basis(vert, target_port);
    colours.push_back(colour);
  }
  // Stores all edges which the vertex can be commuted to
  EdgeVec dest_edges;
  while (true) {
    // The vertex can be commuted to the front
    bool success = true;
    for (unsigned i = 0; i < current_edges.size(); ++i) {
      // Check if the edge is already in the quantum frontier
      if (std::find(
              frontier_edges.begin(), frontier_edges.end(), current_edges[i]) !=
          frontier_edges.end()) {
        dest_edges.push_back(current_edges[i]);
        continue;
      }
      // Check prev_op is a gate
      Vertex prev_vert = circ.source(current_edges[i]);
      Op_ptr prev_op = circ.get_Op_ptr_from_Vertex(prev_vert);
      if (!prev_op->get_desc().is_gate()) {
        // not commute
        return std::nullopt;
      }

      // Check commute
      port_t source_port = circ.get_source_port(current_edges[i]);
      if (!circ.commutes_with_basis(prev_vert, colours[i], source_port)) {
        // not commute
        return std::nullopt;
      } else {
        // Update dest_edges
        Vertex prev_prev_v;
        Edge prev_e;
        std::tie(prev_prev_v, prev_e) =
            circ.get_prev_pair(prev_vert, current_edges[i]);
        dest_edges.push_back(prev_e);
      }
      // Only true if all edges are in frontier
      success = false;
    }
    if (success) {
      std::pair<EdgeVec, EdgeVec> p(initial_edges, dest_edges);
      return p;
    } else {
      current_edges = dest_edges;
      dest_edges = {};
    }
  }
}

static void partial_rewire(
    const Vertex &vert, Circuit &circ, EdgeVec &src_edges,
    EdgeVec &dest_edges) {
  // move the vertex to the frontier
  // Notice that if one of the vertex's in edge is already a destination
  // edge then the circuit::remove_vertex will delete the destination edge
  // hence circuit::rewire would result in an error due to the missing edge.
  // We need a partial rewire for that reason.
  // Example:
  // Moving the second vertex (CX gate) to the front we only need to rewire
  // the "x" part.
  // --o-----
  //   |
  // --x--x--
  //      |
  // -----o--

  for (unsigned i = 0; i < dest_edges.size(); i++) {
    Edge &dest_in_edge = dest_edges[i];
    Edge &curr_in_edge = src_edges[i];
    // If the vertex is already connected to an edge in the frontier, do
    // nothing.
    if (dest_in_edge != curr_in_edge) {
      // Add first edge
      Vertex dest_prev_vert = circ.source(dest_in_edge);
      circ.add_edge(
          {dest_prev_vert, circ.get_source_port(dest_in_edge)},
          {vert, circ.get_target_port(curr_in_edge)}, EdgeType::Quantum);
      // Add second edge
      Vertex curr_next_vert;
      Edge curr_out_edge;
      Vertex dest_next_vert = circ.target(dest_in_edge);
      std::tie(curr_next_vert, curr_out_edge) =
          circ.get_next_pair(vert, curr_in_edge);
      circ.add_edge(
          {vert, circ.get_source_port(curr_out_edge)},
          {dest_next_vert, circ.get_target_port(dest_in_edge)},
          EdgeType::Quantum);
      // Add third edge
      Vertex curr_prev_vert = circ.source(curr_in_edge);
      circ.add_edge(
          {curr_prev_vert, circ.get_source_port(curr_in_edge)},
          {curr_next_vert, circ.get_target_port(curr_out_edge)},
          EdgeType::Quantum);
      // Remove edges
      circ.remove_edge(dest_in_edge);
      circ.remove_edge(curr_in_edge);
      circ.remove_edge(curr_out_edge);
    }
  }
}

bool MultiGateReorder::solve(unsigned max_depth, unsigned max_size) {
  // Assume the frontier has been advanced

  // store a copy of the original this->mapping_frontier_->linear_boundary
  // this object will be updated and reset throughout the procedure
  // so need to return it to original setting at end.
  unit_vertport_frontier_t copy;
  for (const std::pair<UnitID, VertPort> &pair :
       this->mapping_frontier_->linear_boundary->get<TagKey>()) {
    copy.insert({pair.first, pair.second});
  }
  // get subcircuit vertices in topological order
  CutFrontier current_cut = this->mapping_frontier_->circuit_.next_cut(
      frontier_convert_vertport_to_edge(
          this->mapping_frontier_->circuit_,
          this->mapping_frontier_->linear_boundary),
      this->mapping_frontier_->boolean_boundary);
  unsigned subcircuit_depth = 1;
  std::vector<Vertex> subcircuit_vertices(
      current_cut.slice->begin(), current_cut.slice->end());
  // add cuts of vertices to subcircuit_vertices until constraints met, or end
  // of circuit reached
  while (subcircuit_depth < max_depth &&
         unsigned(subcircuit_vertices.size()) < max_size &&
         current_cut.slice->size() > 0) {
    current_cut = this->mapping_frontier_->circuit_.next_cut(
        current_cut.u_frontier, current_cut.b_frontier);
    subcircuit_depth++;
    subcircuit_vertices.insert(
        subcircuit_vertices.end(), current_cut.slice->begin(),
        current_cut.slice->end());
  }
  // For each of the multi-q vertices in the subcircuit, find its unitids
  std::vector<std::pair<Vertex, std::vector<Node>>> vertex_nodes;
  for (const Vertex &vert : subcircuit_vertices) {
    if (!is_multiq_quantum_gate(this->mapping_frontier_->circuit_, vert)) {
      continue;
    }
    std::vector<Node> nodes;
    for (port_t port = 0;
         port < this->mapping_frontier_->circuit_.n_ports(vert); ++port) {
      nodes.push_back(Node(
          get_unitid_from_vertex_port(this->mapping_frontier_, {vert, port})));
    }
    vertex_nodes.push_back({vert, nodes});
  }

  bool modification_made = false;
  for (auto it = vertex_nodes.begin(); it != vertex_nodes.end(); it++) {
    const Vertex &vert = it->first;
    const std::vector<Node> &nodes = it->second;
    // determine whether a vertex has been advanced by the
    // advance_frontier_boundary call in this loop by traversing to the right by
    // max_depth steps to see if it hits a frontier vertex.
    // checking one port is enough
    bool advanced = false;
    VertPort current_vert_port(vert, 0);
    for (unsigned i = 0; i < max_depth; i++) {
      if (this->mapping_frontier_->linear_boundary->get<TagValue>().find(
              current_vert_port) !=
          this->mapping_frontier_->linear_boundary->get<TagValue>().end()) {
        advanced = true;
        break;
      }
      if (this->mapping_frontier_->circuit_.detect_boundary_Op(
              current_vert_port.first)) {
        break;
      }
      Edge current_in_edge = this->mapping_frontier_->circuit_.get_nth_in_edge(
          current_vert_port.first, current_vert_port.second);
      Vertex next_vert;
      Edge next_e;
      std::tie(next_vert, next_e) =
          this->mapping_frontier_->circuit_.get_next_pair(
              current_vert_port.first, current_in_edge);
      current_vert_port = {
          next_vert, this->mapping_frontier_->circuit_.get_target_port(next_e)};
    }
    if (advanced) {
      continue;
    }

    // try to commute a vertex to the frontier if it's physically permitted
    if (this->mapping_frontier_->valid_boundary_operation(
            this->architecture_,
            this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(vert),
            nodes)) {
      std::optional<std::pair<EdgeVec, EdgeVec>> commute_pairs =
          try_find_commute_edges(
              this->mapping_frontier_->circuit_, this->u_frontier_edges_, vert);

      if (commute_pairs != std::nullopt) {
        modification_made = true;
        partial_rewire(
            vert, this->mapping_frontier_->circuit_, (*commute_pairs).first,
            (*commute_pairs).second);
        // Update the frontier
        this->mapping_frontier_->advance_frontier_boundary(this->architecture_);
        this->u_frontier_edges_ =
            convert_u_frontier_to_edges(*frontier_convert_vertport_to_edge(
                this->mapping_frontier_->circuit_,
                this->mapping_frontier_->linear_boundary));
      }
    }
  }
  // Return the quantum boundary to its original setting
  this->mapping_frontier_->set_linear_boundary(copy);
  return modification_made;
}

MultiGateReorderRoutingMethod::MultiGateReorderRoutingMethod(
    unsigned _max_depth, unsigned _max_size)
    : max_depth_(_max_depth), max_size_(_max_size) {}

std::pair<bool, unit_map_t> MultiGateReorderRoutingMethod::routing_method(
    MappingFrontier_ptr &mapping_frontier,
    const ArchitecturePtr &architecture) const {
  MultiGateReorder mr(architecture, mapping_frontier);
  return {mr.solve(this->max_depth_, this->max_size_), {}};
}

unsigned MultiGateReorderRoutingMethod::get_max_depth() const {
  return this->max_depth_;
}

unsigned MultiGateReorderRoutingMethod::get_max_size() const {
  return this->max_size_;
}

nlohmann::json MultiGateReorderRoutingMethod::serialize() const {
  nlohmann::json j;
  j["depth"] = this->max_depth_;
  j["size"] = this->max_size_;
  j["name"] = "MultiGateReorderRoutingMethod";
  return j;
}

MultiGateReorderRoutingMethod MultiGateReorderRoutingMethod::deserialize(
    const nlohmann::json &j) {
  return MultiGateReorderRoutingMethod(
      j.at("depth").get<unsigned>(), j.at("size").get<unsigned>());
}

}  // namespace tket
