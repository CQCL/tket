#include "Mapping/MappingFrontier.hpp"

#include "Circuit/Circuit.hpp"
namespace tket {

/**
 * unit_vertport_frontier_t is <UnitID, VertPort>, helper function returned
 * UnitID given Edge
 */
UnitID get_unitid_from_unit_frontier(
    const std::shared_ptr<unit_vertport_frontier_t>& u_frontier,
    const VertPort& vp) {
  for (auto it = u_frontier->get<TagKey>().begin();
       it != u_frontier->get<TagKey>().end(); ++it) {
    if (it->second == vp) {
      return it->first;
    }
  }
  throw MappingFrontierError(
      std::string("Edge provided not in unit_frontier_t object."));
}

/**
 * quantum_boundary stored as vertport so that correct edge can be recovered
 * after subcircuit substitution method uses Vertex and port_t and
 * Circuit::get_nth_out_edge to generate unit_frontier_t object
 */
std::shared_ptr<unit_frontier_t> frontier_convert_vertport_to_edge(
    const Circuit& circuit,
    const std::shared_ptr<unit_vertport_frontier_t>& u_frontier) {
  // make empty unit_frontier_t object
  std::shared_ptr<unit_frontier_t> output_frontier =
      std::make_shared<unit_frontier_t>();
  // iterate through u_frontier, convert VertPort to Edge and insert
  for (const std::pair<UnitID, VertPort>& pair : u_frontier->get<TagKey>()) {
    output_frontier->insert(
        {pair.first,
         circuit.get_nth_out_edge(pair.second.first, pair.second.second)});
  }
  return output_frontier;
}

/**
 * Initialise quantum_boundary and classical_boundary from
 * out edges of Input vertices
 */
MappingFrontier::MappingFrontier(Circuit& _circuit) : circuit_(_circuit) {
  this->quantum_boundary = std::make_shared<unit_vertport_frontier_t>();
  this->classical_boundary = std::make_shared<b_frontier_t>();
  // Set up {UnitID, VertPort} objects for quantum and classical boundaries
  for (const Qubit& qb : this->circuit_.all_qubits()) {
    this->quantum_boundary->insert({qb, {this->circuit_.get_in(qb), 0}});
  }
  for (const Bit& bit : this->circuit_.all_bits()) {
    this->classical_boundary->insert(
        {bit,
         this->circuit_.get_nth_b_out_bundle(this->circuit_.get_in(bit), 0)});
  }
}

void MappingFrontier::advance_next_2qb_slice(unsigned max_advance) {
  bool boundary_updated = false;
  unsigned loop = 0;
  std::shared_ptr<unit_frontier_t> current_frontier =
      frontier_convert_vertport_to_edge(this->circuit_, this->quantum_boundary);
  // Get all vertices in first cut
  VertexVec immediate_cut_vertices_v =
      *(this->circuit_
            .next_cut(current_frontier, std::make_shared<b_frontier_t>())
            .slice);
  do {
    // each do section first finds the next set of edges after the held set
    // for edges with target vertices with all their edges presented in the
    // first set
    loop++;
    boundary_updated = false;
    // produce next frontier object
    std::shared_ptr<unit_frontier_t> next_frontier =
        std::make_shared<unit_frontier_t>();

    for (const std::pair<UnitID, Edge>& pair :
         current_frontier->get<TagKey>()) {
      // if target_v not in immediate_cut_vertices, then do not pass it
      Vertex target_v = this->circuit_.target(pair.second);

      EdgeVec in_edges =
          this->circuit_.get_in_edges_of_type(target_v, EdgeType::Quantum);

      bool in_slice =
          std::find(
              immediate_cut_vertices_v.begin(), immediate_cut_vertices_v.end(),
              target_v) != immediate_cut_vertices_v.end();

      if (((!in_slice && in_edges.size() > 1) ||
           this->circuit_.get_OpType_from_Vertex(target_v) == OpType::Output) &&
          this->circuit_.get_OpType_from_Vertex(target_v) != OpType::Barrier) {
        // Vertex either not allowed to pass, or is output vertex => update
        // nothing
        next_frontier->insert({pair.first, pair.second});
      } else {
        // vertex can be surpassed, so update quantum_boundary and next_frontier
        // with next edge
        Edge next_edge = this->circuit_.get_next_edge(target_v, pair.second);
        this->quantum_boundary->replace(
            this->quantum_boundary->get<TagKey>().find(pair.first),
            {pair.first,
             {target_v, this->circuit_.get_source_port(next_edge)}});
        next_frontier->insert({pair.first, next_edge});
      }
    }
    // Given new frontier, find the actual next cut
    CutFrontier next_cut = this->circuit_.next_cut(
        next_frontier, std::make_shared<b_frontier_t>());
    // For each vertex in a slice, if its physically permitted, update
    // quantum_boundary with quantum out edges from vertex (i.e.
    // next_cut.u_frontier)
    for (const Vertex& vert : *next_cut.slice) {
      // Output means we don't want to pass, so just leave
      if (this->circuit_.get_OpType_from_Vertex(vert) == OpType::Output) {
        continue;
      }
      EdgeVec in_edges =
          this->circuit_.get_in_edges_of_type(vert, EdgeType::Quantum);
      // More than 1 edge means we want to keep edges, so continue
      if (in_edges.size() > 1) {
        continue;
      }
      // can guarantee that we update now as non-updating cases have been
      // continued
      boundary_updated = true;
      // push edge past single qubit vertex, repeat
      UnitID uid = get_unitid_from_unit_frontier(
          this->quantum_boundary,
          {this->circuit_.source(in_edges[0]),
           this->circuit_.get_source_port(in_edges[0])});

      Edge replacement_edge =
          next_cut.u_frontier->get<TagKey>().find(uid)->second;

      Vertex source_vertex = this->circuit_.source(replacement_edge);
      port_t source_port = this->circuit_.get_source_port(replacement_edge);

      this->quantum_boundary->replace(
          this->quantum_boundary->get<TagKey>().find(uid),
          {uid, {source_vertex, source_port}});
    }
    current_frontier = next_frontier;
  } while (boundary_updated && loop <= max_advance);
  return;
}

/**
 * advance_frontier_boundary
 * terminates when next_cut returns a "slice" where
 * no vertices are physically permitted by the architecture
 * quantum_boundary and classical_boundary updated to reflect this
 */
void MappingFrontier::advance_frontier_boundary(
    const ArchitecturePtr& architecture) {
  bool boundary_updated = false;
  do {
    // next_cut.slice vertices in_edges from this->quantum_boundary
    // TODO: add optional skip function later to skip vertices that don't have
    // physical requirements
    boundary_updated = false;

    std::shared_ptr<unit_frontier_t> frontier_edges =
        frontier_convert_vertport_to_edge(
            this->circuit_, this->quantum_boundary);
    // Add all classical edges that share the same target
    unsigned dummy_bit_index = 0;
    for (const std::pair<UnitID, Edge>& pair : frontier_edges->get<TagKey>()) {
      Vertex vert = this->circuit_.target(pair.second);
      for (const Edge& e :
           this->circuit_.get_in_edges_of_type(vert, EdgeType::Classical)) {
        frontier_edges->insert({Bit(dummy_bit_index), e});
        dummy_bit_index++;
      }
    }

    CutFrontier next_cut = this->circuit_.next_cut(
        frontier_edges, std::make_shared<b_frontier_t>());

    // For each vertex in a slice, if its physically permitted, update
    // quantum_boundary with quantum out edges from vertex (i.e.
    // next_cut.u_frontier)
    for (const Vertex& vert : *next_cut.slice) {
      std::vector<UnitID> uids;
      for (const Edge& e :
           this->circuit_.get_in_edges_of_type(vert, EdgeType::Quantum)) {
        uids.push_back(get_unitid_from_unit_frontier(
            this->quantum_boundary,
            {this->circuit_.source(e), this->circuit_.get_source_port(e)}));
      }

      // TODO: update architecture valid operation to reflect devices supporting
      // different multi qubit operations also, like, think about how best that
      // should actually be done?
      std::vector<Node> nodes;
      for (const UnitID& uid : uids) {
        nodes.push_back(Node(uid));
      }
      if (architecture->valid_operation(
              /* this->circuit_.get_OpType_from_Vertex(vert), */
              nodes) ||
          this->circuit_.get_OpType_from_Vertex(vert) == OpType::Barrier) {
        // if no valid operation, boundary not updated and while loop terminates
        boundary_updated = true;
        for (const UnitID& uid : uids) {
          Edge replacement_edge =
              next_cut.u_frontier->get<TagKey>().find(uid)->second;
          Vertex source_vertex = this->circuit_.source(replacement_edge);
          port_t source_port = this->circuit_.get_source_port(replacement_edge);
          this->quantum_boundary->replace(
              this->quantum_boundary->get<TagKey>().find(uid),
              {uid, {source_vertex, source_port}});
        }
      }
    }
  } while (boundary_updated);
  return;
}

/**
 * convert_u_frontier_to_edges
 * Subcircuit requires EdgeVec, not unit_frontier_t as boundary information
 * Helper Functions to convert types
 * TODO: also probably another way of doing this? EdgeVec required for
 * subcircuit. Double check with someone who knows better than I...
 */
EdgeVec convert_u_frontier_to_edges(const unit_frontier_t& u_frontier) {
  EdgeVec edges;
  for (const std::pair<UnitID, Edge>& pair : u_frontier.get<TagKey>()) {
    edges.push_back(pair.second);
  }
  return edges;
}

Subcircuit MappingFrontier::get_frontier_subcircuit(
    unsigned _max_subcircuit_depth, unsigned _max_subcircuit_size) const {
  CutFrontier current_cut = this->circuit_.next_cut(
      frontier_convert_vertport_to_edge(this->circuit_, this->quantum_boundary),
      this->classical_boundary);

  unsigned subcircuit_depth = 1;
  VertexSet subcircuit_vertices(
      current_cut.slice->begin(), current_cut.slice->end());
  // add cuts of vertices to subcircuit_vertices until constraints met, or end
  // of circuit reached
  while (subcircuit_depth < _max_subcircuit_depth &&
         unsigned(subcircuit_vertices.size()) < _max_subcircuit_size &&
         current_cut.slice->size() > 0) {
    current_cut =
        this->circuit_.next_cut(current_cut.u_frontier, current_cut.b_frontier);
    subcircuit_depth++;
    subcircuit_vertices.insert(
        current_cut.slice->begin(), current_cut.slice->end());
  }
  if (subcircuit_vertices.size() == 0) {
    throw MappingFrontierError("Subcircuit being produced with no gates.");
  }
  return Subcircuit(
      convert_u_frontier_to_edges(*frontier_convert_vertport_to_edge(
          this->circuit_, this->quantum_boundary)),
      convert_u_frontier_to_edges(*current_cut.u_frontier),
      subcircuit_vertices);
}

// TODO: Update to support ancillas
void MappingFrontier::update_quantum_boundary_uids(
    const unit_map_t& relabelled_uids) {
  for (const std::pair<const UnitID, UnitID>& label : relabelled_uids) {
    // implies new labelling
    if (label.first != label.second) {
      // by type, label.first already assumed in circuit
      // this condition means label.second also in circuit
      // implies that a merging is done -> remove first qubit
      if (this->quantum_boundary->get<TagKey>().find(label.second) !=
          this->quantum_boundary->get<TagKey>().end()) {
        // erase, assume updated already
        this->quantum_boundary->erase(label.first);
      } else {
        auto current_label_it =
            this->quantum_boundary->get<TagKey>().find(label.first);
        // relabel "label.first" with "label.second"
        this->quantum_boundary->replace(
            current_label_it, {label.second, current_label_it->second});
        unit_map_t relabel = {label};
        this->circuit_.rename_units(relabel);
      }
    }
  }
}

// TODO: expects every qubit is present in permutation, even if unmoved
void MappingFrontier::permute_subcircuit_q_out_hole(
    const unit_map_t& final_permutation, Subcircuit& subcircuit) {
  EdgeVec new_q_out_hole;
  int i = 0;
  // Change to iterate through final permutation first?
  if (this->quantum_boundary->size() != final_permutation.size()) {
    throw MappingFrontierError(
        "Number of Qubits in mapping permutation does not match number of "
        "Qubits in MappingFrontier boundary, for permuting Qubits as with "
        "routed Subcircuit.");
  }
  for (const std::pair<UnitID, VertPort>& pair :
       this->quantum_boundary->get<TagKey>()) {
    // other iteration avoids this...
    // TODO: change this when making route different subcircuits
    auto it = final_permutation.find(pair.first);
    if (it == final_permutation.end()) {
      throw MappingFrontierError("Qubit in boundary not in permutation.");
    }
    std::pair<UnitID, UnitID> uid_pair = *it;
    if (uid_pair.first == uid_pair.second) {
      new_q_out_hole.push_back(subcircuit.q_out_hole[i]);
    } else {
      int j = 0;
      for (auto it = this->quantum_boundary->get<TagKey>().begin();
           it != this->quantum_boundary->get<TagKey>().end(); ++it) {
        if (it->first == uid_pair.second) {
          new_q_out_hole.push_back(subcircuit.q_out_hole[j]);
          break;
        }
        j++;
      }
    }
    i++;
  }
  subcircuit.q_out_hole = new_q_out_hole;
}

/**
 * MappingFrontier::get_u_frontier_default_unit_map
 * Map from default qubit register qubits to UnitIDs in quantum_boundary
 */
unit_map_t MappingFrontier::get_default_to_quantum_boundary_unit_map() const {
  unsigned i = 0;
  unit_map_t default_to_u_frontier_map;
  for (const std::pair<UnitID, VertPort>& pair :
       this->quantum_boundary->get<TagKey>()) {
    default_to_u_frontier_map.insert({Qubit(i), pair.first});
    i++;
  }
  return default_to_u_frontier_map;
}

void MappingFrontier::set_quantum_boundary(
    const unit_vertport_frontier_t& new_boundary) {
  this->quantum_boundary = std::make_shared<unit_vertport_frontier_t>();
  for (const std::pair<UnitID, VertPort>& pair : new_boundary.get<TagKey>()) {
    this->quantum_boundary->insert(pair);
  }
}

/**
 * add_qubit
 * Adds given UnitID as a qubit to held circuit.
 * Updates boundary.
 */
void MappingFrontier::add_qubit(const UnitID& uid) {
  Qubit qb(uid);
  this->circuit_.add_qubit(qb);
  this->quantum_boundary->insert({qb, {this->circuit_.get_in(qb), 0}});
}

/**
 * add_swap
 * Inserts an OpType::SWAP gate into the uid_0 and uid_1 edges held in
 * quantum_boundary This directly modifies circuit_ Updates quantum_boundary to
 * reflect new edges
 */
void MappingFrontier::add_swap(const UnitID& uid_0, const UnitID& uid_1) {
  // get iterators to quantum_boundary uids
  auto uid0_in_it = this->quantum_boundary->find(uid_0);
  auto uid1_in_it = this->quantum_boundary->find(uid_1);

  // Add Qubit if not in MappingFrontier boundary (i.e. not in circuit)
  if (uid0_in_it == this->quantum_boundary->end()) {
    this->add_qubit(uid_0);
    uid0_in_it = this->quantum_boundary->find(uid_0);
  }
  if (uid1_in_it == this->quantum_boundary->end()) {
    this->add_qubit(uid_1);
    uid1_in_it = this->quantum_boundary->find(uid_1);
  }

  // update held ancillas
  Node n0 = Node(uid_0);
  Node n1 = Node(uid_1);

  bool uid0_ancilla =
      this->ancilla_nodes_.find(n0) != this->ancilla_nodes_.end();
  bool uid1_ancilla =
      this->ancilla_nodes_.find(n1) != this->ancilla_nodes_.end();

  if (uid0_ancilla && !uid1_ancilla) {
    this->ancilla_nodes_.erase(n0);
    this->ancilla_nodes_.insert(n1);
  }
  if (!uid0_ancilla && uid1_ancilla) {
    this->ancilla_nodes_.erase(n1);
    this->ancilla_nodes_.insert(n0);
  }

  // Get predecessor edges to SWAP insert location
  VertPort vp0 = uid0_in_it->second;
  VertPort vp1 = uid1_in_it->second;
  EdgeVec predecessors = {
      this->circuit_.get_nth_out_edge(vp0.first, vp0.second),
      this->circuit_.get_nth_out_edge(vp1.first, vp1.second)};

  // add SWAP vertex to circuit_ and rewire into predecessor
  Vertex swap_v = this->circuit_.add_vertex(OpType::SWAP);
  this->circuit_.rewire(
      swap_v, predecessors, {EdgeType::Quantum, EdgeType::Quantum});

  // Update boundary to reflect new edges
  EdgeVec successors = this->circuit_.get_all_out_edges(swap_v);
  this->circuit_.dag[successors[0]].ports.first = 1;
  this->circuit_.dag[successors[1]].ports.first = 0;

  this->quantum_boundary->replace(
      uid0_in_it, {uid_0, {this->circuit_.source(successors[1]), 0}});
  this->quantum_boundary->replace(
      uid1_in_it, {uid_1, {this->circuit_.source(successors[0]), 1}});

  // update output vertices of quantum boundary of circuit to reflect changing
  // qubit paths
  auto uid0_circuit_boundary_it =
      this->circuit_.boundary.get<TagID>().find(uid_0);
  auto uid1_circuit_boundary_it =
      this->circuit_.boundary.get<TagID>().find(uid_1);

  Vertex uid0_out = uid0_circuit_boundary_it->out_;
  Vertex uid1_out = uid1_circuit_boundary_it->out_;
  Vertex uid0_in = uid0_circuit_boundary_it->in_;
  Vertex uid1_in = uid1_circuit_boundary_it->in_;

  this->circuit_.boundary.get<TagID>().erase(uid_0);
  this->circuit_.boundary.get<TagID>().erase(uid_1);

  this->circuit_.boundary.get<TagID>().insert({uid_0, uid0_in, uid1_out});
  this->circuit_.boundary.get<TagID>().insert({uid_1, uid1_in, uid0_out});
}

void MappingFrontier::add_bridge(
    const UnitID& control, const UnitID& central, const UnitID& target) {
  // get predecessors
  auto control_in_it = this->quantum_boundary->find(control);
  auto central_in_it = this->quantum_boundary->find(central);
  auto target_in_it = this->quantum_boundary->find(target);

  // by virtue of method, control and target qubit will always be in BRIDGE.
  // However, distances used to check BRIDGE and find PATH may use
  // central qubit that is unallocated, in which add it.
  if (central_in_it == this->quantum_boundary->end()) {
    this->add_qubit(central);
    central_in_it = this->quantum_boundary->find(central);
  }

  VertPort vp_control = control_in_it->second;
  VertPort vp_central = central_in_it->second;
  VertPort vp_target = target_in_it->second;

  EdgeVec predecessors = {
      this->circuit_.get_nth_out_edge(vp_control.first, vp_control.second),
      this->circuit_.get_nth_out_edge(vp_central.first, vp_central.second),
      this->circuit_.get_nth_out_edge(vp_target.first, vp_target.second),
  };  // get cx vertex
      // this should be guaranteeds by pre-checks
  Vertex cx_v = this->circuit_.target(predecessors[0]);
  // add bridge
  Vertex bridge_v = this->circuit_.add_vertex(OpType::BRIDGE);
  // add bridge vertex to circuit
  this->circuit_.rewire(
      bridge_v, predecessors,
      {EdgeType::Quantum, EdgeType::Quantum, EdgeType::Quantum});
  // remove old cx vertex
  this->circuit_.remove_vertex(
      cx_v, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);
}

}  // namespace tket
