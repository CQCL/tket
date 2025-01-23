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

////////////////////////////////////////////////////
// ALL METHODS TO OBTAIN COMPLEX GRAPH INFORMATIION//
////////////////////////////////////////////////////

#include <tklog/TketLog.hpp>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Circuit/Slices.hpp"
#include "tket/OpType/EdgeType.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/Ops/OpPtr.hpp"

namespace tket {

bool Circuit::is_simple() const {
  if (!default_regs_ok()) return false;
  for (const BoundaryElement& el : boundary.get<TagID>()) {
    std::string reg = el.id_.reg_name();
    if (!(reg == q_default_reg() || reg == c_default_reg())) return false;
  }
  return true;
}

bool Circuit::default_regs_ok() const {
  opt_reg_info_t q_info = get_reg_info(q_default_reg());
  register_info_t correct_q_info = {UnitType::Qubit, 1};
  if (q_info && q_info.value() != correct_q_info) return false;
  opt_reg_info_t c_info = get_reg_info(c_default_reg());
  register_info_t correct_c_info = {UnitType::Bit, 1};
  if (c_info && c_info.value() != correct_c_info) return false;
  return true;
}

std::map<OpType, unsigned> Circuit::op_counts() const {
  std::map<OpType, unsigned> counts;
  BGL_FORALL_VERTICES(v, dag, DAG) { counts[get_OpType_from_Vertex(v)]++; }
  return counts;
}

unsigned Circuit::count_gates(
    const OpType& op_type, const bool include_conditional) const {
  unsigned counter = 0;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    if (get_OpType_from_Vertex(v) == op_type) {
      ++counter;
    } else if (
        include_conditional &&
        (get_OpType_from_Vertex(v) == OpType::Conditional) &&
        (static_cast<const Conditional&>(*get_Op_ptr_from_Vertex(v))
             .get_op()
             ->get_type() == op_type)) {
      ++counter;
    }
  }
  return counter;
}

VertexSet Circuit::get_gates_of_type(const OpType& op_type) const {
  VertexSet vset;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    if (get_OpType_from_Vertex(v) == op_type) {
      vset.insert(v);
    }
  }
  return vset;
}

std::list<Command> Circuit::get_commands_of_type(OpType op_type) const {
  std::list<Command> coms;
  for (const Command& cmd : *this) {
    if (cmd.get_op_ptr()->get_type() == op_type) {
      coms.push_back(cmd);
    }
  }
  return coms;
}
unsigned Circuit::count_n_qubit_gates(unsigned size) const {
  unsigned counter = 0;
  if (size == 0) return counter;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    if (n_in_edges_of_type(v, EdgeType::Quantum) == size) {
      const Op_ptr op_ptr = get_Op_ptr_from_Vertex(v);
      switch (op_ptr->get_type()) {
        case OpType::Input:
        case OpType::Create:
        case OpType::Output:
        case OpType::Discard:
        case OpType::Reset:
        case OpType::Measure:
        case OpType::Barrier: {
          break;
        }
        default: {
          counter++;
        }
      }
    }
  }
  return counter;
}

Circuit Circuit::subcircuit(const Subcircuit& sc) const {
  Circuit sub;
  vertex_map_t vmap;
  VertexVec q_ins;
  VertexVec q_outs;
  VertexVec c_ins;
  VertexVec c_outs;
  std::map<Edge, Vertex> in_boundary_map;
  std::map<Edge, Vertex> out_boundary_map;
  for (const Edge& e : sc.q_in_hole) {
    Vertex added = sub.add_vertex(OpType::Input);
    vmap[source(e)] = added;
    q_ins.push_back(added);
    in_boundary_map.insert({e, added});
  }
  for (const Edge& e : sc.q_out_hole) {
    Vertex added = sub.add_vertex(OpType::Output);
    vmap[target(e)] = added;
    q_outs.push_back(added);
    out_boundary_map.insert({e, added});
  }
  for (const Edge& e : sc.c_in_hole) {
    Vertex added = sub.add_vertex(OpType::ClInput);
    vmap[source(e)] = added;
    c_ins.push_back(added);
    in_boundary_map.insert({e, added});
  }
  for (const Edge& e : sc.c_out_hole) {
    Vertex added = sub.add_vertex(OpType::ClOutput);
    vmap[target(e)] = added;
    c_outs.push_back(added);
    out_boundary_map.insert({e, added});
  }
  for (unsigned i = 0; i < q_ins.size(); i++) {
    sub.boundary.insert({Qubit(i), q_ins[i], q_outs[i]});
  }
  for (unsigned i = 0; i < c_ins.size(); i++) {
    sub.boundary.insert({Bit(i), c_ins[i], c_outs[i]});
  }
  for (const Vertex& v : sc.verts) {
    Vertex added = sub.add_vertex(get_Op_ptr_from_Vertex(v));
    vmap[v] = added;
  }
  for (const Vertex& v : sc.verts) {
    // iterate through original circuit as order of the set varies between Mac
    // and Docker
    BGL_FORALL_INEDGES(v, e, dag, DAG) {
      Vertex source = this->source(e);
      Vertex sub_source = vmap.at(source);
      port_t in_port = get_source_port(e);
      OpType type = sub.get_OpType_from_Vertex(sub_source);
      if (is_initial_q_type(type) || type == OpType::ClInput) {
        Edge boundary_edge = get_linear_edge(e);
        // Multiple inputs might be mapped to the same source
        // so need to distinguish them.
        sub_source = in_boundary_map.at(boundary_edge);
        in_port = 0;
      }
      sub.add_edge(
          {sub_source, in_port}, {vmap.at(v), get_target_port(e)}, dag[e].type);
    }
  }
  for (const Edge& e : sc.q_out_hole) {
    // Multiple outputs might be mapped to the same target
    // so need to distinguish them.
    Vertex out = out_boundary_map[e];
    Vertex sub_source = vmap.at(source(e));
    port_t in_port = get_source_port(e);
    std::map<Edge, Vertex>::iterator found = in_boundary_map.find(e);
    if (found != in_boundary_map.end()) {
      sub_source = found->second;
      in_port = 0;
    }
    sub.add_edge({sub_source, in_port}, {out, 0}, EdgeType::Quantum);
  }
  for (const Edge& e : sc.c_out_hole) {
    Vertex out = out_boundary_map[e];
    Vertex sub_source = vmap.at(source(e));
    port_t in_port = get_source_port(e);
    std::map<Edge, Vertex>::iterator found = in_boundary_map.find(e);
    if (found != in_boundary_map.end()) {
      sub_source = found->second;
      in_port = 0;
    }
    sub.add_edge({sub_source, in_port}, {out, 0}, EdgeType::Classical);
  }
  return sub;
}

Subcircuit Circuit::singleton_subcircuit(const Vertex& v) const {
  return {
      get_in_edges_of_type(v, EdgeType::Quantum),
      get_out_edges_of_type(v, EdgeType::Quantum),
      get_in_edges_of_type(v, EdgeType::Classical),
      get_out_edges_of_type(v, EdgeType::Classical),
      get_out_edges_of_type(v, EdgeType::Boolean),
      {v}};
}

// returns qubit path via vertices & inhabited port in vertices
// used to construct a routing grid
QPathDetailed Circuit::unit_path(const UnitID& unit) const {
  Vertex current_v = get_in(unit);

  QPathDetailed path = {{current_v, 0}};
  Edge betweenEdge = get_nth_out_edge(current_v, 0);
  current_v = target(betweenEdge);

  while (detect_final_Op(current_v) == false) {
    if (n_out_edges(current_v) == 0) {
      throw CircuitInvalidity("A path ends before reaching an output vertex.");
    }
    port_t n = get_target_port(betweenEdge);
    VertPort v_and_port = {current_v, n};
    path.push_back(v_and_port);
    betweenEdge = get_nth_out_edge(current_v, n);
    current_v = target(betweenEdge);
  }
  path.push_back({current_v, 0});
  return path;
}

// returns a vector of each qubits path via qubit_path
// this is all the information required to make a circuit
std::vector<QPathDetailed> Circuit::all_qubit_paths() const {
  std::vector<QPathDetailed> new_list_of_paths;
  for (const Qubit& q : all_qubits()) {
    new_list_of_paths.push_back(unit_path(q));
  }
  return new_list_of_paths;
}

std::map<UnitID, QPathDetailed> Circuit::all_unit_paths() const {
  std::map<UnitID, QPathDetailed> new_list_of_paths;
  for (const Qubit& q : all_qubits()) {
    new_list_of_paths.insert({q, unit_path(q)});
  }
  for (const Bit& b : all_bits()) {
    new_list_of_paths.insert({b, unit_path(b)});
  }
  return new_list_of_paths;
}

// safely return boundary iterator to UnitID
boundary_t::iterator boundary_elem(const Circuit& circ, const UnitID& unit) {
  boundary_t::iterator found = circ.boundary.get<TagID>().find(unit);

  if (found == circ.boundary.get<TagID>().end()) {
    throw CircuitInvalidity("Unit not found in circuit: " + unit.repr());
  }
  return found;
}

/*
    Permute output boundary of circuit according to qubit map
    Assumes all circuit Qubits are mapped
*/
void Circuit::permute_boundary_output(const qubit_map_t& qm) {
  std::map<UnitID, BoundaryElement> new_entries;

  for (const auto& qb_pair : qm) {
    boundary_t::iterator input = boundary_elem(*this, qb_pair.first);
    boundary_t::iterator output = boundary_elem(*this, qb_pair.second);
    new_entries.insert({output->id_, {output->id_, output->in_, input->out_}});
  }

  for (const auto& pair : new_entries) {
    boundary.erase(boundary.get<TagID>().find(pair.first));
  }

  for (const auto& pair : new_entries) {
    boundary.insert(pair.second);
  }
}

qubit_map_t Circuit::implicit_qubit_permutation() const {
  qubit_map_t perm;
  for (const QPathDetailed& path : all_qubit_paths()) {
    Qubit in_q(get_id_from_in(path.front().first));
    Qubit out_q(get_id_from_out(path.back().first));
    perm.insert({in_q, out_q});
  }
  return perm;
}

bool Circuit::has_implicit_wireswaps() const {
  qubit_map_t perm = implicit_qubit_permutation();
  for (const std::pair<const Qubit, Qubit>& pair : perm) {
    if (pair.first != pair.second) return true;
  }
  return false;
}

// returns a basic qubit path consisting of just vertices
// this was used in early methods for printing
// it is now only used for MatrixAnalysis methods
// TODO:: remove from circuit
VertexVec Circuit::qubit_path_vertices(const Qubit& qubit) const {
  QPathDetailed path = unit_path(qubit);
  VertexVec follow_q;
  for (const VertPort& pair : path) {
    follow_q.push_back(pair.first);
  }
  return follow_q;
}

// returns 'slices' of 'parallel' actions in dag as a vector encompassing all
// vertices
// requires the boundaries to be correct and the circuit to be fully connected
SliceVec Circuit::get_slices() const {
  SliceVec slices;
  for (SliceIterator sit = this->slice_begin(); sit != slice_end(); ++sit) {
    slices.push_back(*sit);
  }
  return slices;
}

Edge Circuit::skip_irrelevant_edges(Edge current) const {
  Vertex try_next_v = target(current);
  while (n_out_edges_of_type(try_next_v, EdgeType::Quantum) == 1) {
    std::tie(try_next_v, current) = get_next_pair(try_next_v, current);
  }
  return current;
}

SliceVec Circuit::get_reverse_slices() const {
  vertex_map_t mapping;
  vertex_map_t rev_mapping;
  Circuit rev;
  for (const BoundaryElement& el : boundary.get<TagID>()) {
    Vertex new_in, new_out;
    if (el.type() == UnitType::Qubit) {
      new_in = rev.add_vertex(OpType::Input);
      new_out = rev.add_vertex(OpType::Output);
    } else {
      new_in = rev.add_vertex(OpType::ClInput);
      new_out = rev.add_vertex(OpType::ClOutput);
    }
    mapping[el.in_] = new_out;
    rev_mapping[new_out] = el.in_;
    mapping[el.out_] = new_in;
    rev_mapping[new_in] = el.out_;
    rev.boundary.insert({el.id_, new_in, new_out});
  }
  BGL_FORALL_VERTICES(v, dag, DAG) {
    const Op_ptr op_ptr = get_Op_ptr_from_Vertex(v);
    switch (op_ptr->get_type()) {
      case OpType::Input:
      case OpType::Create:
      case OpType::Output:
      case OpType::ClInput:
      case OpType::ClOutput:
      case OpType::Discard:
      case OpType::WASMInput:
      case OpType::WASMOutput: {
        break;
      }
      default: {
        Vertex v0 = rev.add_vertex(op_ptr);
        mapping[v] = v0;
        rev_mapping[v0] = v;
        break;
      }
    }
  }
  BGL_FORALL_EDGES(e, dag, DAG) {
    Vertex s = source(e);
    port_t sp = get_source_port(e);
    Vertex t = target(e);
    port_t tp = get_target_port(e);
    EdgeType type = dag[e].type;
    if (type == EdgeType::Boolean) {
      // Move Boolean to read from bit wire just before next write
      Edge bit_wire = get_nth_out_edge(s, sp);
      Vertex next_on_bit = target(bit_wire);
      port_t next_p = get_target_port(bit_wire);
      rev.add_edge({mapping[next_on_bit], next_p}, {mapping[t], tp}, type);
    } else {
      rev.add_edge({mapping[t], tp}, {mapping[s], sp}, type);
    }
  }
  SliceVec slices_of_rev = rev.get_slices();
  SliceVec rev_slices;
  for (const Slice& s : slices_of_rev) {
    Slice sl;
    for (const Vertex& v : s) {
      sl.push_back(rev_mapping[v]);
    }
    rev_slices.push_back(sl);
  }
  return rev_slices;
}

unsigned Circuit::depth() const {
  unsigned count = 0;
  std::function<bool(Op_ptr)> skip_func = [&](Op_ptr op) {
    return (op->get_type() == OpType::Barrier);
  };
  SliceIterator slice_iter(*this, skip_func);
  if (!(*slice_iter).empty()) count++;
  while (!slice_iter.finished()) {
    slice_iter.cut_ = this->next_cut(
        slice_iter.cut_.u_frontier, slice_iter.cut_.b_frontier, skip_func);
    if (!(*slice_iter).empty()) count++;
  }
  return count;
}

unsigned Circuit::depth_by_type(OpType _type) const {
  unsigned count = 0;
  std::function<bool(Op_ptr)> skip_func = [&](Op_ptr op) {
    return (op->get_type() != _type);
  };
  SliceIterator slice_iter(*this, skip_func);
  if (!(*slice_iter).empty()) count++;
  while (!slice_iter.finished()) {
    slice_iter.cut_ = this->next_cut(
        slice_iter.cut_.u_frontier, slice_iter.cut_.b_frontier, skip_func);
    if (!(*slice_iter).empty()) count++;
  }
  return count;
}

unsigned Circuit::depth_by_types(const OpTypeSet& _types) const {
  unsigned count = 0;
  std::function<bool(Op_ptr)> skip_func = [&](Op_ptr op) {
    return (_types.find(op->get_type()) == _types.end());
  };
  SliceIterator slice_iter(*this, skip_func);
  if (!(*slice_iter).empty()) count++;
  while (!slice_iter.finished()) {
    slice_iter.cut_ = this->next_cut(
        slice_iter.cut_.u_frontier, slice_iter.cut_.b_frontier, skip_func);
    if (!(*slice_iter).empty()) count++;
  }
  return count;
}

unsigned Circuit::depth_2q() const {
  unsigned count = 0;
  std::function<bool(Op_ptr)> skip_func = [&](Op_ptr op) {
    return (op->n_qubits() != 2 || op->get_type() == OpType::Barrier);
  };
  SliceIterator slice_iter(*this, skip_func);
  if (!(*slice_iter).empty()) count++;
  while (!slice_iter.finished()) {
    slice_iter.cut_ = this->next_cut(
        slice_iter.cut_.u_frontier, slice_iter.cut_.b_frontier, skip_func);
    if (!(*slice_iter).empty()) count++;
  }
  return count;
}

std::map<Vertex, unit_set_t> Circuit::vertex_unit_map() const {
  std::map<Vertex, unit_set_t> map;
  BGL_FORALL_VERTICES(v, dag, DAG) { map[v] = {}; }
  for (const std::pair<const UnitID, QPathDetailed>& path : all_unit_paths()) {
    for (const VertPort& vp : path.second) {
      map[vp.first].insert(path.first);
    }
  }
  return map;
}

std::map<Vertex, unsigned> Circuit::vertex_depth_map() const {
  std::map<Vertex, unsigned> map;
  unsigned i = 0;
  for (SliceIterator it = slice_begin(); it != slice_end(); ++it) {
    for (const Vertex& v : *it) {
      map[v] = i;
    }
    i++;
  }
  for (const BoundaryElement& el : boundary) {
    map[el.in_] = 0;
    map[el.out_] = i;
  }
  return map;
}

std::map<Vertex, unsigned> Circuit::vertex_rev_depth_map() const {
  std::map<Vertex, unsigned> map;
  SliceVec slices = get_reverse_slices();
  for (unsigned i = 0; i < slices.size(); i++) {
    for (const Vertex& v : slices[i]) {
      map[v] = i;
    }
  }
  for (const BoundaryElement& el : boundary) {
    map[el.in_] = slices.size();
    map[el.out_] = 0;
  }
  return map;
}

std::map<Edge, UnitID> Circuit::edge_unit_map() const {
  std::map<Edge, UnitID> map;
  for (const std::pair<const UnitID, QPathDetailed>& path : all_unit_paths()) {
    QPathDetailed::const_iterator it = path.second.begin();
    for (++it; it != path.second.end(); ++it) {
      map.insert({get_nth_in_edge(it->first, it->second), path.first});
    }
  }
  return map;
}

Circuit::CommandIterator::CommandIterator(const Circuit& circ)
    : current_slice_iterator_(circ.slice_begin()),
      current_index_(0),
      circ_(&circ) {
  if ((*current_slice_iterator_).size() == 0) {
    *this = circ.end();
  } else {
    current_vertex_ = (*current_slice_iterator_)[0];
    current_command_ = circ.command_from_vertex(
        current_vertex_, current_slice_iterator_.get_u_frontier(),
        current_slice_iterator_.get_prev_b_frontier());
  }
}

Circuit::CommandIterator Circuit::begin() const {
  return CommandIterator(*this);
}

const Circuit::CommandIterator Circuit::end() const { return nullcit; }

const Circuit::CommandIterator Circuit::nullcit = CommandIterator();

Circuit::CommandIterator::Commandholder Circuit::CommandIterator::operator++(
    int) {
  Commandholder ret(current_command_);
  ++*this;
  return ret;
}

Circuit::CommandIterator& Circuit::CommandIterator::operator++() {
  if (*this == circ_->end()) return *this;
  if (current_index_ == (*current_slice_iterator_).size() - 1) {
    if (current_slice_iterator_.finished()) {
      *this = circ_->end();
      return *this;
    }
    ++current_slice_iterator_;
    current_index_ = 0;
  } else {
    ++current_index_;
  }
  if (current_index_ == (*current_slice_iterator_).size()) {
    TKET_ASSERT(!"slice is empty");
  }
  current_vertex_ = (*current_slice_iterator_)[current_index_];
  current_command_ = circ_->command_from_vertex(
      current_vertex_, current_slice_iterator_.get_u_frontier(),
      current_slice_iterator_.get_prev_b_frontier());
  return *this;
}

unit_vector_t Circuit::args_from_frontier(
    const Vertex& vert, std::shared_ptr<const unit_frontier_t> u_frontier,
    std::shared_ptr<const b_frontier_t> prev_b_frontier) const {
  EdgeVec ins = get_in_edges(vert);
  unit_vector_t args;
  for (port_t p = 0; p < ins.size(); ++p) {
    switch (get_edgetype(ins[p])) {
      case EdgeType::WASM: {
        Edge out = get_next_edge(vert, ins[p]);
        bool found = false;
        for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
          if (pair.second == out) {
            args.push_back(pair.first);
            found = true;
            break;
          }
        }
        TKET_ASSERT(found);  // Vertex edges not found in frontier.
        break;
      }
      case EdgeType::Boolean: {
        bool found = false;
        for (const std::pair<Bit, EdgeVec>& pair :
             prev_b_frontier->get<TagKey>()) {
          for (const Edge& edge : pair.second) {
            if (edge == ins[p]) {
              args.push_back(pair.first);
              found = true;
              break;
            }
          }
          if (found) {
            break;
          }
        }
        TKET_ASSERT(found);  // Vertex edges not found in Boolean frontier.
        break;
      }
      case EdgeType::Classical:
      case EdgeType::Quantum: {
        Edge out = get_next_edge(vert, ins[p]);
        bool found = false;
        for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
          if (pair.second == out) {
            args.push_back(pair.first);
            found = true;
            break;
          }
        }
        if (!found)
          throw CircuitInvalidity(
              "Vertex edges not found in frontier. Edge: " +
              get_Op_ptr_from_Vertex(source(out))->get_name() + " -> " +
              get_Op_ptr_from_Vertex(target(out))->get_name());
        break;
      }
      default: {
        TKET_ASSERT(!"args_from_frontier found invalid edge type in signature");
      }
    }
  }
  return args;
}

Command Circuit::command_from_vertex(
    const Vertex& vert, std::shared_ptr<const unit_frontier_t> u_frontier,
    std::shared_ptr<const b_frontier_t> prev_b_frontier) const {
  unit_vector_t args = args_from_frontier(vert, u_frontier, prev_b_frontier);
  return Command(
      get_Op_ptr_from_Vertex(vert), args, get_opgroup_from_Vertex(vert), vert);
}

}  // namespace tket
