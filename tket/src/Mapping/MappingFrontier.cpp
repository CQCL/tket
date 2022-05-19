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

#include "Mapping/MappingFrontier.hpp"

#include "Circuit/Circuit.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

/**
 * unit_vertport_frontier_t is <UnitID, VertPort>, helper function returns
 * UnitID corresponding to given VertPort
 */
UnitID get_unitid_from_unit_frontier(
    const std::shared_ptr<unit_vertport_frontier_t>& u_frontier,
    const VertPort& vp) {
  auto it = u_frontier->get<TagValue>().find(vp);
  if (it != u_frontier->get<TagValue>().end()) return it->first;
  throw MappingFrontierError(
      std::string("Edge provided not in unit_frontier_t object."));
}

/**
 * bit_frontier_t is <Bit, EdgeVec>, helper function returns
 * Bit corresponding to given Edge
 */
static Bit get_bit_from_bool_frontier(
    const std::shared_ptr<b_frontier_t>& b_frontier, const EdgeVec& ev) {
  TKET_ASSERT(ev.size() > 0);
  // condition is that if one Edge in EdgeVector ev is in a
  // held bundle, then all are so return true
  for (auto it = b_frontier->get<TagKey>().begin();
       it != b_frontier->get<TagKey>().end(); ++it) {
    for (const Edge& e0 : ev) {
      for (const Edge& e1 : it->second) {
        if (e0 == e1) {
          return it->first;
        }
      }
    }
  }
  throw MappingFrontierError(
      std::string("EdgeVec provided not in b_frontier_t object."));
}

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
 * Initialise linear_boundary and boolean_boundary from
 * out edges of Input vertices
 */
MappingFrontier::MappingFrontier(Circuit& _circuit) : circuit_(_circuit) {
  this->linear_boundary = std::make_shared<unit_vertport_frontier_t>();
  this->boolean_boundary = std::make_shared<b_frontier_t>();
  this->bimaps_ = std::make_shared<unit_bimaps_t>();

  // Set up {UnitID, VertPort} objects for quantum and classical boundaries
  for (const Qubit& qb : this->circuit_.all_qubits()) {
    this->linear_boundary->insert({qb, {this->circuit_.get_in(qb), 0}});
    this->bimaps_->initial.insert({qb, qb});
    this->bimaps_->final.insert({qb, qb});
  }
  for (const Bit& bit : this->circuit_.all_bits()) {
    Vertex bit_input = this->circuit_.get_in(bit);
    EdgeVec bool_bundle = this->circuit_.get_nth_b_out_bundle(bit_input, 0);
    if (bool_bundle.size() != 0) {
      this->boolean_boundary->insert({bit, bool_bundle});
    }
    if (this->circuit_.n_out_edges_of_type(bit_input, EdgeType::Classical) >
        0) {
      this->linear_boundary->insert({bit, {bit_input, 0}});
    }
  }
}

/**
 * Initialise linear_boundary and boolean_boundary from
 * out edges of Input vertices
 */
MappingFrontier::MappingFrontier(
    Circuit& _circuit, std::shared_ptr<unit_bimaps_t> _bimaps)
    : circuit_(_circuit), bimaps_(_bimaps) {
  // Check that the maps are valid
  for (const Qubit& q : _circuit.all_qubits()) {
    if (_bimaps->initial.right.find(q) == _bimaps->initial.right.end()) {
      throw MappingFrontierError(
          "Uid " + q.repr() + " not found in initial map.");
    }
    if (_bimaps->final.right.find(q) == _bimaps->final.right.end()) {
      throw MappingFrontierError(
          "Uid " + q.repr() + " not found in final map.");
    }
  }

  this->linear_boundary = std::make_shared<unit_vertport_frontier_t>();
  this->boolean_boundary = std::make_shared<b_frontier_t>();

  // Set up {UnitID, VertPort} objects for quantum and classical boundaries
  for (const Qubit& qb : this->circuit_.all_qubits()) {
    this->linear_boundary->insert({qb, {this->circuit_.get_in(qb), 0}});
  }
  for (const Bit& bit : this->circuit_.all_bits()) {
    Vertex bit_input = this->circuit_.get_in(bit);
    EdgeVec bool_bundle = this->circuit_.get_nth_b_out_bundle(bit_input, 0);
    if (bool_bundle.size() != 0) {
      this->boolean_boundary->insert({bit, bool_bundle});
    } else {
      this->linear_boundary->insert({bit, {bit_input, 0}});
    }
  }
}

MappingFrontier::MappingFrontier(const MappingFrontier& mapping_frontier)
    : circuit_(mapping_frontier.circuit_), bimaps_(mapping_frontier.bimaps_) {
  this->linear_boundary = std::make_shared<unit_vertport_frontier_t>();
  this->boolean_boundary = std::make_shared<b_frontier_t>();

  for (const std::pair<UnitID, VertPort>& pair :
       mapping_frontier.linear_boundary->get<TagKey>()) {
    this->linear_boundary->insert({pair.first, pair.second});
  }
  for (const std::pair<Bit, EdgeVec>& pair :
       mapping_frontier.boolean_boundary->get<TagKey>()) {
    EdgeVec edges;
    for (const Edge& edge : pair.second) {
      edges.push_back(edge);
    }
    this->boolean_boundary->insert({pair.first, edges});
  }
  for (const Node& node : mapping_frontier.ancilla_nodes_) {
    this->ancilla_nodes_.insert(node);
  }
}

void MappingFrontier::advance_next_2qb_slice(unsigned max_advance) {
  bool boundary_updated = false;
  unsigned loop = 0;
  std::shared_ptr<unit_frontier_t> current_frontier =
      frontier_convert_vertport_to_edge(this->circuit_, this->linear_boundary);

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
      OpType ot = this->circuit_.get_OpType_from_Vertex(target_v);
      if (((!in_slice && in_edges.size() > 1) || ot == OpType::Output ||
           ot == OpType::ClOutput) &&
          this->circuit_.get_OpType_from_Vertex(target_v) != OpType::Barrier) {
        // Vertex either not allowed to pass, or is output vertex => update
        // nothing
        next_frontier->insert({pair.first, pair.second});
      } else {
        // vertex can be surpassed, so update linear_boundary and
        // next_frontier with next edge
        Edge next_edge = this->circuit_.get_next_edge(target_v, pair.second);
        this->linear_boundary->replace(
            this->linear_boundary->get<TagKey>().find(pair.first),
            {pair.first,
             {target_v, this->circuit_.get_source_port(next_edge)}});
        next_frontier->insert({pair.first, next_edge});
      }
    }
    // Given new frontier, find the actual next cut
    CutFrontier next_cut = this->circuit_.next_cut(
        next_frontier, std::make_shared<b_frontier_t>());
    // For each vertex in a slice, if its physically permitted, update
    // linear_boundary with quantum out edges from vertex (i.e.
    // next_cut.u_frontier)
    for (const Vertex& vert : *next_cut.slice) {
      // Output means we don't want to pass, so just leave
      OpType ot = this->circuit_.get_OpType_from_Vertex(vert);
      if (ot == OpType::Output || ot == OpType::ClOutput) {
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
          this->linear_boundary, {this->circuit_.source(in_edges[0]),
                                  this->circuit_.get_source_port(in_edges[0])});

      Edge replacement_edge =
          next_cut.u_frontier->get<TagKey>().find(uid)->second;

      Vertex source_vertex = this->circuit_.source(replacement_edge);
      port_t source_port = this->circuit_.get_source_port(replacement_edge);

      this->linear_boundary->replace(
          this->linear_boundary->get<TagKey>().find(uid),
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
 * linear_boundary and boolean_boundary updated to reflect this
 */
void MappingFrontier::advance_frontier_boundary(
    const ArchitecturePtr& architecture) {
  bool boundary_updated = false;
  do {
    // next_cut.slice vertices in_edges from this->linear_boundary
    boundary_updated = false;
    std::shared_ptr<unit_frontier_t> l_frontier_edges =
        frontier_convert_vertport_to_edge(
            this->circuit_, this->linear_boundary);

    CutFrontier next_cut =
        this->circuit_.next_cut(l_frontier_edges, this->boolean_boundary);
    // For each vertex in a slice, if its physically permitted, update
    // linear_boundary with quantum out edges from vertex (i.e.
    // next_cut.u_frontier)
    // update boolean_boundary in line
    for (const Vertex& vert : *next_cut.slice) {
      // for each boolean edge into vertex, collect associated Bit and port
      // number n.b. a single Bit may have multiple "in bundles" to different
      // vertices in the same cut
      std::map<Bit, port_t> bool_uid_port_set;
      std::vector<EdgeVec> b_in_bundles = this->circuit_.get_b_in_bundles(vert);
      for (unsigned i = 0; i < b_in_bundles.size(); i++) {
        EdgeVec ev = b_in_bundles[i];
        if (ev.size() > 0) {
          bool_uid_port_set.insert(
              {get_bit_from_bool_frontier(this->boolean_boundary, ev), i});
        }
      }

      // for each quantum edge into vertex, collect associated Qubit/Node
      // don't collect port as this->linear_boundary holds this
      // each UnitID will only have one quantum edge active
      std::vector<UnitID> l_uids;  // linear unit id
      std::vector<Node> nodes;     // quantum/node only
      for (const Edge& e :
           this->circuit_.get_in_edges_of_type(vert, EdgeType::Quantum)) {
        UnitID uid = get_unitid_from_unit_frontier(
            this->linear_boundary,
            {this->circuit_.source(e), this->circuit_.get_source_port(e)});
        l_uids.push_back(uid);
        nodes.push_back(Node(uid));
      }

      // for each classical edge store related UnitID in l_uids
      // each Bit will only have one classical edge active
      // n.b. some vertices may introduce new Bit to the boolean_boundary
      // e.g. A Measurement result may be passed to a conditional as
      // therefore, all Bit not in the Boolean boundary are also stored
      // in case the operation does this and the this->boolean_boundary
      // needs to be updated
      std::map<Bit, port_t> extra_bool_uid_port_set;
      for (const Edge& e :
           this->circuit_.get_in_edges_of_type(vert, EdgeType::Classical)) {
        // for updating linear boundary
        port_t port_source = this->circuit_.get_source_port(e);
        UnitID uid = get_unitid_from_unit_frontier(
            this->linear_boundary, {this->circuit_.source(e), port_source});
        l_uids.push_back(uid);

        // for potentially adding new Bit to boolean boundary
        // port_target makes it possible to track which "out bundle" corresponds
        // to this Bit
        port_t port_target = this->circuit_.get_target_port(e);
        Bit bit = Bit(uid);
        if (bool_uid_port_set.find(bit) == bool_uid_port_set.end()) {
          extra_bool_uid_port_set.insert({bit, port_target});
        }
      }

      if (nodes.size() == 0 ||
          this->valid_boundary_operation(
              architecture, this->circuit_.get_Op_ptr_from_Vertex(vert),
              nodes)) {
        // if no valid operation, boundary not updated and while loop terminates
        boundary_updated = true;
        // update linear UnitID (Qubits&Quantum edges, Bits&Classical edges)
        for (const UnitID& uid : l_uids) {
          Edge replacement_edge =
              next_cut.u_frontier->get<TagKey>().find(uid)->second;
          Vertex source_vertex = this->circuit_.source(replacement_edge);
          port_t source_port = this->circuit_.get_source_port(replacement_edge);
          this->linear_boundary->replace(
              this->linear_boundary->get<TagKey>().find(uid),
              {uid, {source_vertex, source_port}});
        }
        // update booleans
        // n.b. its possible a boolean path terminates with an operation
        // however, the port should be preserved so we can track the correct Bit
        // {Bit, port_t}
        for (auto it = bool_uid_port_set.begin(); it != bool_uid_port_set.end();
             ++it) {
          std::vector<EdgeVec> out_bundles =
              this->circuit_.get_b_out_bundles(vert);

          port_t port = it->second;
          TKET_ASSERT(out_bundles.size() > port);  // safe port indexing
          // However, this Bit may have boolean values in other Vertices in
          // slice therefore, we remove every edge from the vertex in_bundle for
          // this port from the boolean_boundary and then insert these new edges
          std::vector<EdgeVec> in_bundles =
              this->circuit_.get_b_in_bundles(vert);

          TKET_ASSERT(in_bundles.size() > port);  // safe port indexing
          EdgeVec in_bundle = in_bundles[port];
          // Bit should be in boolean_boundary
          Bit bit = it->first;
          auto jt = this->boolean_boundary->get<TagKey>().find(bit);
          TKET_ASSERT(jt != this->boolean_boundary->get<TagKey>().end());
          // construct a new EdgeVec object with replaced Edge and persisting
          // Edge

          EdgeVec new_boolean_edges;

          for (const Edge& e : jt->second) {
            // => edge isn't being replaced
            if (std::find(in_bundle.begin(), in_bundle.end(), e) ==
                in_bundle.end()) {
              new_boolean_edges.push_back(e);
            }
          }
          // add all new edges from out bundle to the boundary
          for (const Edge& e : out_bundles[port]) {
            new_boolean_edges.push_back(e);
          }

          // boolean no longer needed
          if (new_boolean_edges.size() == 0) {
            this->boolean_boundary->erase(jt);
          } else {
            // replace boolean boundary
            this->boolean_boundary->replace(jt, {bit, new_boolean_edges});
          }
        }
        // Some operations may spawn a Boolean wire not held in boolean_boundary
        // this checks for any new wires and if true, adds to boolean_boundary
        // {Bit, port_t}
        for (auto it = extra_bool_uid_port_set.begin();
             it != extra_bool_uid_port_set.end(); ++it) {
          std::vector<EdgeVec> source_out =
              this->circuit_.get_b_out_bundles(vert);
          // If source_out has more bundles than port value, then we know
          // it's been spawned (multiple could be spawned at same vertex)
          port_t port = it->second;
          if (source_out.size() > port) {
            EdgeVec new_boolean_wire = source_out[port];
            // add new edges to boolean_boundary
            // note that a boolean cannot be spawned in multiple vertices
            // as the incoming Bit wire is linear
            // Measure always create a boolean, even if empty of edges
            // => check size before adding
            if (new_boolean_wire.size() > 0) {
              this->boolean_boundary->insert({it->first, new_boolean_wire});
            }
          }
        }
      }
    }
  } while (boundary_updated);
  return;
}

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
      frontier_convert_vertport_to_edge(this->circuit_, this->linear_boundary),
      this->boolean_boundary);

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
  TKET_ASSERT(subcircuit_vertices.size() != 0);
  return Subcircuit(
      convert_u_frontier_to_edges(*frontier_convert_vertport_to_edge(
          this->circuit_, this->linear_boundary)),
      convert_u_frontier_to_edges(*current_cut.u_frontier),
      subcircuit_vertices);
}

UnitID MappingFrontier::get_qubit_from_circuit_uid(const UnitID& uid) {
  auto it = this->bimaps_->initial.right.find(uid);
  if (it == this->bimaps_->initial.right.end()) {
    throw MappingFrontierError("UnitID not found in initial map.");
  }
  return it->second;
}

void MappingFrontier::update_bimaps(UnitID qubit, UnitID node) {
  // Update initial map
  auto init_it = this->bimaps_->initial.left.find(qubit);
  if (init_it == this->bimaps_->initial.left.end())
    throw MappingFrontierError("Qubit not found in initial map.");
  this->bimaps_->initial.left.erase(init_it);
  this->bimaps_->initial.left.insert({qubit, node});
  // Update final map
  auto final_it = this->bimaps_->final.left.find(qubit);
  if (final_it == this->bimaps_->final.left.end())
    throw MappingFrontierError("Qubit not found in final map.");
  this->bimaps_->final.left.erase(final_it);
  this->bimaps_->final.left.insert({qubit, node});
}

void MappingFrontier::update_linear_boundary_uids(
    const unit_map_t& relabelled_uids) {
  for (const std::pair<const UnitID, UnitID>& label : relabelled_uids) {
    // implies new labelling
    if (label.first != label.second) {
      // by type, label.first already assumed in circuit
      // this condition means label.second also in circuit
      // implies that a merging is done -> remove first qubit

      if (this->linear_boundary->get<TagKey>().find(label.second) !=
          this->linear_boundary->get<TagKey>().end()) {
        // erase, assume updated already
        this->linear_boundary->erase(label.first);
      } else {
        auto current_label_it =
            this->linear_boundary->get<TagKey>().find(label.first);
        // relabel "label.first" with "label.second"
        this->linear_boundary->replace(
            current_label_it, {label.second, current_label_it->second});
        unit_map_t relabel = {label};
        this->circuit_.rename_units(relabel);
      }
    }
  }
}

void MappingFrontier::permute_subcircuit_q_out_hole(
    const unit_map_t& final_permutation, Subcircuit& subcircuit) {
  EdgeVec new_q_out_hole;
  int i = 0;
  // Change to iterate through final permutation first?
  if (this->linear_boundary->size() != final_permutation.size()) {
    throw MappingFrontierError(
        "Number of Qubits in mapping permutation does not match number of "
        "Qubits in MappingFrontier boundary, for permuting Qubits as with "
        "routed Subcircuit.");
  }
  for (const std::pair<UnitID, VertPort>& pair :
       this->linear_boundary->get<TagKey>()) {
    auto it = final_permutation.find(pair.first);
    if (it == final_permutation.end()) {
      throw MappingFrontierError("Qubit in boundary not in permutation.");
    }
    std::pair<UnitID, UnitID> uid_pair = *it;
    if (uid_pair.first == uid_pair.second) {
      new_q_out_hole.push_back(subcircuit.q_out_hole[i]);
    } else {
      int j = 0;
      for (auto it = this->linear_boundary->get<TagKey>().begin();
           it != this->linear_boundary->get<TagKey>().end(); ++it) {
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
 * Map from default qubit register qubits to UnitIDs in linear_boundary
 */
unit_map_t MappingFrontier::get_default_to_linear_boundary_unit_map() const {
  unsigned i = 0;
  unit_map_t default_to_u_frontier_map;
  for (const std::pair<UnitID, VertPort>& pair :
       this->linear_boundary->get<TagKey>()) {
    default_to_u_frontier_map.insert({Qubit(i), pair.first});
    i++;
  }
  return default_to_u_frontier_map;
}

void MappingFrontier::set_linear_boundary(
    const unit_vertport_frontier_t& new_boundary) {
  this->linear_boundary = std::make_shared<unit_vertport_frontier_t>();
  for (const std::pair<UnitID, VertPort>& pair : new_boundary.get<TagKey>()) {
    this->linear_boundary->insert(pair);
  }
}

/**
 * add_swap
 * Inserts an OpType::SWAP gate into the uid_0 and uid_1 edges held in
 * linear_boundary This directly modifies circuit_ Updates linear_boundary to
 * reflect new edges
 */
bool MappingFrontier::add_swap(const UnitID& uid_0, const UnitID& uid_1) {
  // get iterators to linear_boundary uids
  auto uid0_in_it = this->linear_boundary->find(uid_0);
  auto uid1_in_it = this->linear_boundary->find(uid_1);

  // Add Qubit if not in MappingFrontier boundary (i.e. not in circuit)
  if (uid0_in_it == this->linear_boundary->end()) {
    this->add_ancilla(uid_0);
    uid0_in_it = this->linear_boundary->find(uid_0);
  }
  if (uid1_in_it == this->linear_boundary->end()) {
    this->add_ancilla(uid_1);
    uid1_in_it = this->linear_boundary->find(uid_1);
  }

  // update held ancillas
  // the location/id of the "ancilla node" changes when a SWAP occurs
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

  // If the immediate vertex before the new vertex is the same SWAP, this
  // implies a rut has been hit In which case return false that SWAP can't be
  // added Can safely do this check here after relabelling as adding/updating
  // ancillas => fresh SWAP
  Vertex source0 = this->circuit_.source(predecessors[0]);
  Vertex source1 = this->circuit_.source(predecessors[1]);
  if (source0 == source1) {
    if (this->circuit_.get_OpType_from_Vertex(source0) == OpType::SWAP) {
      return false;
    }
  }

  // add SWAP vertex to circuit_ and rewire into predecessor
  Vertex swap_v = this->circuit_.add_vertex(OpType::SWAP);
  this->circuit_.rewire(
      swap_v, predecessors, {EdgeType::Quantum, EdgeType::Quantum});

  // Update boundary to reflect new edges
  EdgeVec successors = this->circuit_.get_all_out_edges(swap_v);
  this->circuit_.dag[successors[0]].ports.first = 1;
  this->circuit_.dag[successors[1]].ports.first = 0;

  this->linear_boundary->replace(
      uid0_in_it, {uid_0, {this->circuit_.source(successors[1]), 0}});
  this->linear_boundary->replace(
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

  std::map<Node, Node> final_map = {{n0, n1}, {n1, n0}};

  update_maps(this->bimaps_, {}, final_map);
  return true;
}

void MappingFrontier::add_bridge(
    const UnitID& control, const UnitID& central, const UnitID& target) {
  // get predecessors
  auto control_in_it = this->linear_boundary->find(control);
  auto central_in_it = this->linear_boundary->find(central);
  auto target_in_it = this->linear_boundary->find(target);

  // by virtue of method, control and target qubit will always be in BRIDGE.
  // However, distances used to check BRIDGE and find PATH may use
  // central qubit that is unallocated, in which add it.
  if (central_in_it == this->linear_boundary->end()) {
    this->add_ancilla(central);
    central_in_it = this->linear_boundary->find(central);
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

void MappingFrontier::add_ancilla(const UnitID& ancilla) {
  Qubit qb(ancilla);
  this->circuit_.add_qubit(qb);
  this->linear_boundary->insert({qb, {this->circuit_.get_in(qb), 0}});

  this->bimaps_->initial.insert({qb, qb});
  this->bimaps_->final.insert({qb, qb});
  this->ancilla_nodes_.insert(Node(ancilla));
  UnitID uid_ancilla(ancilla);

  unit_map_t update_map;
  update_map.insert({uid_ancilla, uid_ancilla});

  update_maps(this->bimaps_, update_map, update_map);
}

void MappingFrontier::merge_ancilla(
    const UnitID& merge, const UnitID& ancilla) {
  // get output and input vertices
  Vertex merge_v_in = this->circuit_.get_in(merge);
  Vertex merge_v_out = this->circuit_.get_out(merge);
  Vertex ancilla_v_out = this->circuit_.get_out(ancilla);
  // find source vertex & port of merge_v_out
  // output vertex, so can assume single edge
  Edge merge_out_edge = this->circuit_.get_nth_out_edge(merge_v_in, 0);
  Edge ancilla_in_edge = this->circuit_.get_nth_in_edge(ancilla_v_out, 0);
  // Find port number
  port_t merge_target_port = this->circuit_.get_target_port(merge_out_edge);
  port_t ancilla_source_port = this->circuit_.get_source_port(ancilla_in_edge);
  // Find vertices
  Vertex merge_v_target = this->circuit_.target(merge_out_edge);
  Vertex ancilla_v_source = this->circuit_.source(ancilla_in_edge);

  // remove and replace edges
  this->circuit_.remove_edge(merge_out_edge);
  this->circuit_.remove_edge(ancilla_in_edge);
  this->circuit_.add_edge(
      {ancilla_v_source, ancilla_source_port},
      {merge_v_target, merge_target_port}, EdgeType::Quantum);

  // instead of manually updating all boundaries, we change which output
  // vertex the qubit paths to
  Edge merge_in_edge = this->circuit_.get_nth_in_edge(merge_v_out, 0);
  port_t merge_source_port = this->circuit_.get_source_port(merge_in_edge);
  Vertex merge_v_source = this->circuit_.source(merge_in_edge);

  this->circuit_.remove_edge(merge_in_edge);
  this->circuit_.add_edge(
      {merge_v_source, merge_source_port}, {ancilla_v_out, 0},
      EdgeType::Quantum);

  // remove empty vertex wire, relabel dag vertices
  this->circuit_.dag[merge_v_in].op = get_op_ptr(OpType::noop);
  this->circuit_.dag[merge_v_out].op = get_op_ptr(OpType::noop);
  this->circuit_.remove_vertex(
      merge_v_in, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  this->circuit_.remove_vertex(
      merge_v_out, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

  // Can now just erase "merge" qubit from the circuit
  this->circuit_.boundary.get<TagID>().erase(merge);

  // Update the qubit mappings
  // let's call the arguments ancilla_node and merge_node
  // e.g. before merge:
  //  initial := {ancilla_q:node_x, merge_q:some_uid}
  //  final := {ancilla_q:ancilla_node, merge_q:merge_node}
  // e.g. after merge:
  //  initial := {merge_q:node_x}
  //  final := {merge_q:ancilla_node}
  // Basically, in both qubit maps, erase the entry with qubit merge_q
  // then replace the entry ancilla_q -> x with the merge_q -> x

  auto merge_it = this->bimaps_->initial.right.find(merge);
  TKET_ASSERT(merge_it != this->bimaps_->initial.right.end());
  UnitID merge_q = merge_it->second;
  this->bimaps_->initial.right.erase(merge_it);
  this->bimaps_->final.left.erase(merge_q);
  // Find ancilla_q
  auto final_it = this->bimaps_->final.right.find(ancilla);
  UnitID ancilla_q = final_it->second;
  // Replace in final map
  this->bimaps_->final.right.erase(final_it);
  this->bimaps_->final.left.insert({merge_q, ancilla});
  // Replace in initial map
  auto init_it = this->bimaps_->initial.left.find(ancilla_q);
  UnitID init_ancilla_node = init_it->second;
  this->bimaps_->initial.left.erase(init_it);
  this->bimaps_->initial.left.insert({merge_q, init_ancilla_node});
}

bool MappingFrontier::valid_boundary_operation(
    const ArchitecturePtr& architecture, const Op_ptr& op,
    const std::vector<Node>& uids) const {
  // boxes are never allowed
  OpType ot = op->get_type();
  if (is_box_type(ot)) {
    return false;
  }

  if (ot == OpType::Conditional) {
    Op_ptr cond_op_ptr = static_cast<const Conditional&>(*op).get_op();
    // conditional boxes are never allowed, too
    OpType ot = cond_op_ptr->get_type();
    while (ot == OpType::Conditional) {
      cond_op_ptr = static_cast<const Conditional&>(*op).get_op();
      ot = cond_op_ptr->get_type();
      if (is_box_type(ot)) {
        return false;
      }
    }
  }

  // Barriers are allways allowed
  if (ot == OpType::Barrier) {
    return true;
  }

  // this currently allows unplaced single qubits gates
  // this should be changes in the future
  if (uids.size() == 1) {
    return true;
  }

  // allow two qubit gates only for placed and connected nodes
  if (uids.size() == 2) {
    bool n0 = architecture->node_exists(uids[0]);
    bool n1 = architecture->node_exists(uids[1]);
    if (n0 && n1) {
      bool bde = architecture->bidirectional_edge_exists(uids[0], uids[1]);
      if (bde) {
        return true;
      }
    }
  } else if (uids.size() == 3 && ot == OpType::BRIDGE) {
    bool con_0_exists =
        architecture->bidirectional_edge_exists(uids[0], uids[1]);
    bool con_1_exists =
        architecture->bidirectional_edge_exists(uids[2], uids[1]);
    if (architecture->node_exists(uids[0]) &&
        architecture->node_exists(uids[1]) &&
        architecture->node_exists(uids[2]) && con_0_exists && con_1_exists) {
      return true;
    }
  }

  return false;
}

}  // namespace tket
