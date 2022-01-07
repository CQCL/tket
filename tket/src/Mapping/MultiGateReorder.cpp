#include "Mapping/MultiGateReorder.hpp"

#include "Mapping/MappingFrontier.hpp"

namespace tket {

MultiGateReorder::MultiGateReorder(
    const ArchitecturePtr &_architecture,
    std::shared_ptr<MappingFrontier> &_mapping_frontier)
    : architecture_(_architecture), mapping_frontier_(_mapping_frontier) {
  // This needs to be updated every time the frontier changes
  this->u_frontier_edges_ =
      convert_u_frontier_to_edges(*frontier_convert_vertport_to_edge(
          _mapping_frontier->circuit_, _mapping_frontier->quantum_boundary));
}

// Traverse the DAG to find the UnitID associated with an edge
Vertex get_input_from_vertex_edge(
    Circuit &circ, const Vertex &current_vertex, const Edge &current_outedge) {
  Edge current_e = current_outedge;
  Vertex current_v = current_vertex;
  while (true) {
    if (is_initial_q_type(circ.get_OpType_from_Vertex(current_v))) {
      return current_v;
    }
    std::tie(current_v, current_e) = circ.get_prev_pair(current_v, current_e);
  }
}

// This method will try to commute a vertex to the quantum frontier
bool MultiGateReorder::try_commute_multi_to_front(const Vertex &vert) {
  // Initialize to be the in_edges for the given vertex
  EdgeVec current_edges = this->mapping_frontier_->circuit_.get_in_edges(vert);
  // Stores all edges which the vertex can be commuted to
  EdgeVec dest_edges;
  while (true) {
    // The vertex can be commuted to the front
    bool success = true;
    for (const Edge &in_edge : current_edges) {
      // Check if the edge is already in the quantum frontier
      if (std::find(
              this->u_frontier_edges_.begin(), this->u_frontier_edges_.end(),
              in_edge) != this->u_frontier_edges_.end()) {
        dest_edges.push_back(in_edge);
        continue;
      }
      // Check prev_op is a gate
      Vertex prev_vert = this->mapping_frontier_->circuit_.source(in_edge);
      Op_ptr prev_op =
          this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(prev_vert);
      if (!prev_op->get_desc().is_gate()) {
        // not commute
        return false;
      }

      // Check commute
      port_t source_port =
          this->mapping_frontier_->circuit_.get_source_port(in_edge);
      port_t target_port =
          this->mapping_frontier_->circuit_.get_target_port(in_edge);
      Op_ptr current_op =
          this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(vert);
      std::optional<Pauli> current_colour =
          current_op->commuting_basis(target_port);
      if (!prev_op->commutes_with_basis(current_colour, source_port)) {
        // not commute
        return false;
      } else {
        // Update dest_edges
        Vertex prev_prev_v;
        Edge prev_e;
        std::tie(prev_prev_v, prev_e) =
            this->mapping_frontier_->circuit_.get_prev_pair(prev_vert, in_edge);
        dest_edges.push_back(prev_e);
      }
      // Only true if all edges are in frontier
      success = false;
    }
    if (success) {
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

      EdgeVec initial_in_edges =
          this->mapping_frontier_->circuit_.get_in_edges(vert);
      for (unsigned i = 0; i < dest_edges.size(); i++) {
        Edge &dest_in_edge = dest_edges[i];
        Edge &curr_in_edge = initial_in_edges[i];
        // If the vertex is already connected to an edge in the frontier, do
        // nothing.
        if (dest_in_edge != curr_in_edge) {
          // Add first edge
          Vertex dest_prev_vert =
              this->mapping_frontier_->circuit_.source(dest_in_edge);
          this->mapping_frontier_->circuit_.add_edge(
              {dest_prev_vert,
               this->mapping_frontier_->circuit_.get_source_port(dest_in_edge)},
              {vert,
               this->mapping_frontier_->circuit_.get_target_port(curr_in_edge)},
              EdgeType::Quantum);
          // Add second edge
          Vertex curr_next_vert;
          Edge curr_out_edge;
          Vertex dest_next_vert =
              this->mapping_frontier_->circuit_.target(dest_in_edge);
          std::tie(curr_next_vert, curr_out_edge) =
              this->mapping_frontier_->circuit_.get_next_pair(
                  vert, curr_in_edge);
          this->mapping_frontier_->circuit_.add_edge(
              {vert, this->mapping_frontier_->circuit_.get_source_port(
                         curr_out_edge)},
              {dest_next_vert,
               this->mapping_frontier_->circuit_.get_target_port(dest_in_edge)},
              EdgeType::Quantum);
          // Add third edge
          Vertex curr_prev_vert =
              this->mapping_frontier_->circuit_.source(curr_in_edge);
          this->mapping_frontier_->circuit_.add_edge(
              {curr_prev_vert,
               this->mapping_frontier_->circuit_.get_source_port(curr_in_edge)},
              {curr_next_vert,
               this->mapping_frontier_->circuit_.get_target_port(
                   curr_out_edge)},
              EdgeType::Quantum);
          // Remove edges
          this->mapping_frontier_->circuit_.remove_edge(dest_in_edge);
          this->mapping_frontier_->circuit_.remove_edge(curr_in_edge);
          this->mapping_frontier_->circuit_.remove_edge(curr_out_edge);
        }
      }
      return true;
    } else {
      current_edges = dest_edges;
      dest_edges = {};
    }
  }
}

void MultiGateReorder::solve(unsigned max_depth, unsigned max_size) {
  // Assume the frontier has been advanced

  // store a copy of the original this->mapping_frontier_->quantum_boundray
  // this object will be updated and reset throughout the procedure
  // so need to return it to original setting at end
  unit_vertport_frontier_t copy;
  for (const std::pair<UnitID, VertPort> &pair :
       this->mapping_frontier_->quantum_boundary->get<TagKey>()) {
    copy.insert({pair.first, pair.second});
  }
  // Get a subcircuit only for iterating vertices
  Subcircuit circ =
      this->mapping_frontier_->get_frontier_subcircuit(max_depth, max_size);
  for (const Vertex &vert : circ.verts) {
    // Check if the vertex is:
    //  1. physically permitted
    //  2. is a multi qubit quantum operation without classical controls
    Op_ptr op = this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(vert);
    std::vector<Node> nodes;
    for (const Edge &edge :
         this->mapping_frontier_->circuit_.get_all_out_edges(vert)) {
      Vertex input = get_input_from_vertex_edge(
          this->mapping_frontier_->circuit_, vert, edge);
      nodes.push_back(
          Node(this->mapping_frontier_->circuit_.get_id_from_in(input)));
    }
    if (op->get_desc().is_gate() &&
        this->mapping_frontier_->circuit_.n_in_edges(vert) > 1 &&
        this->mapping_frontier_->circuit_.n_in_edges_of_type(
            vert, EdgeType::Quantum) ==
            this->mapping_frontier_->circuit_.n_in_edges(vert) &&
        this->mapping_frontier_->circuit_.n_out_edges_of_type(
            vert, EdgeType::Quantum) ==
            this->mapping_frontier_->circuit_.n_out_edges(vert) &&
        this->architecture_->valid_operation(nodes)) {
      if (this->try_commute_multi_to_front(vert)) {
        // Update the frontier
        this->mapping_frontier_->advance_frontier_boundary(this->architecture_);
        this->u_frontier_edges_ =
            convert_u_frontier_to_edges(*frontier_convert_vertport_to_edge(
                this->mapping_frontier_->circuit_,
                this->mapping_frontier_->quantum_boundary));
      }
    }
  }
  // Return the quantum boundary to its original setting
  this->mapping_frontier_->set_quantum_boundary(copy);
}

MultiGateReorderRoutingMethod::MultiGateReorderRoutingMethod(
    unsigned _max_depth, unsigned _max_size)
    : max_depth_(_max_depth), max_size_(_max_size){};

bool MultiGateReorderRoutingMethod::check_method(
    const std::shared_ptr<MappingFrontier> & /*mapping_frontier*/,
    const ArchitecturePtr & /*architecture*/) const {
  return true;
}

unit_map_t MultiGateReorderRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier> &mapping_frontier,
    const ArchitecturePtr &architecture) const {
  MultiGateReorder mr(architecture, mapping_frontier);
  mr.solve(this->max_depth_, this->max_size_);
  return {};
}

}  // namespace tket
