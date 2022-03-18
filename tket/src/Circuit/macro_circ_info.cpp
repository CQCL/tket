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

////////////////////////////////////////////////////
// ALL METHODS TO OBTAIN COMPLEX GRAPH INFORMATIION//
////////////////////////////////////////////////////

#include "Circuit.hpp"
#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/TketLog.hpp"

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

unsigned Circuit::count_gates(const OpType& op_type) const {
  unsigned counter = 0;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    if (get_OpType_from_Vertex(v) == op_type) {
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
  std::function<bool(Op_ptr)> skip_func = [=](Op_ptr op) {
    return (op->get_type() != op_type);
  };
  Circuit::SliceIterator slice_iter(*this, skip_func);
  for (const Vertex& v : *slice_iter) {
    coms.push_back(command_from_vertex(
        v, slice_iter.get_u_frontier(), slice_iter.get_prev_b_frontier()));
  }
  while (!slice_iter.finished()) {
    slice_iter.cut_ = next_cut(
        slice_iter.cut_.u_frontier, slice_iter.cut_.b_frontier, skip_func);
    for (const Vertex& v : *slice_iter) {
      coms.push_back(command_from_vertex(
          v, slice_iter.get_u_frontier(), slice_iter.get_prev_b_frontier()));
    }
  }
  return coms;
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
        // For Quantum and Classical edges, boundary_edge == *it;
        // for Boolean, this gives the corresponding Classical
        Edge boundary_edge = this->get_nth_out_edge(source, in_port);
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

static std::shared_ptr<unit_frontier_t> get_next_u_frontier(
    const Circuit& circ, std::shared_ptr<const unit_frontier_t> u_frontier,
    const VertexSet& next_slice_lookup) {
  std::shared_ptr<unit_frontier_t> next_frontier =
      std::make_shared<unit_frontier_t>();
  for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
    Vertex next_v = circ.target(pair.second);
    if (next_slice_lookup.find(next_v) == next_slice_lookup.end()) {
      next_frontier->insert(pair);
    } else {
      next_frontier->insert(
          {pair.first, circ.get_next_edge(next_v, pair.second)});
    }
  }
  return next_frontier;
}

static std::shared_ptr<b_frontier_t> get_next_b_frontier(
    const Circuit& circ, std::shared_ptr<const b_frontier_t> b_frontier,
    std::shared_ptr<const unit_frontier_t> u_frontier,
    const VertexSet& next_slice_lookup) {
  std::shared_ptr<b_frontier_t> next_b_frontier =
      std::make_shared<b_frontier_t>();
  // Copy any remaining edges
  for (const std::pair<Bit, EdgeVec>& pair : b_frontier->get<TagKey>()) {
    EdgeVec remaining;
    for (const Edge& e : pair.second) {
      Vertex targ = circ.target(e);
      if (next_slice_lookup.find(targ) == next_slice_lookup.end()) {
        remaining.push_back(e);
      }
    }
    if (!remaining.empty()) {
      next_b_frontier->insert({pair.first, remaining});
    }
  }
  // Add any new bits introduced in this slice
  for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
    if (circ.get_edgetype(pair.second) == EdgeType::Quantum) continue;
    Vertex next_v = circ.target(pair.second);
    if (next_slice_lookup.find(next_v) == next_slice_lookup.end()) continue;
    if (next_b_frontier->get<TagKey>().find(Bit(pair.first)) !=
        next_b_frontier->end()) {
      throw CircuitInvalidity("RAW hazard created in slicing");
    }
    port_t p = circ.get_target_port(pair.second);
    EdgeVec reads = circ.get_nth_b_out_bundle(next_v, p);
    if (!reads.empty()) next_b_frontier->insert({Bit(pair.first), reads});
  }
  return next_b_frontier;
}

CutFrontier Circuit::next_cut(
    std::shared_ptr<const unit_frontier_t> u_frontier,
    std::shared_ptr<const b_frontier_t> b_frontier) const {
  auto next_slice = std::make_shared<Slice>();
  VertexSet next_slice_lookup;
  VertexSet bad_vertices;
  std::list<Edge> all_edges;
  EdgeSet edge_lookup;
  for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
    if (pair.first.type() == UnitType::Bit) {
      Vertex targ = target(pair.second);
      b_frontier_t::const_iterator found =
          b_frontier->get<TagKey>().find(Bit(pair.first));
      if (found != b_frontier->get<TagKey>().end()) {
        bool still_live = false;
        for (const Edge& e : found->second) {
          if (target(e) != targ) {
            still_live = true;
            break;
          }
        }
        if (still_live) continue;
      }
    }
    all_edges.push_back(pair.second);
    edge_lookup.insert(pair.second);
  }
  for (const std::pair<Bit, EdgeVec>& pair : b_frontier->get<TagKey>()) {
    for (const Edge& e : pair.second) {
      all_edges.push_back(e);
      edge_lookup.insert(e);
    }
  }
  // find the next slice first
  for (const Edge& e : all_edges) {
    Vertex try_v = target(e);
    if (detect_final_Op(try_v)) continue;
    if (next_slice_lookup.find(try_v) != next_slice_lookup.end())
      continue;  // already going to be in next slice
    bool good_vertex = bad_vertices.find(try_v) == bad_vertices.end();
    if (!good_vertex) continue;
    EdgeVec ins = get_in_edges(try_v);
    for (const Edge& in : ins) {
      if (edge_lookup.find(in) == edge_lookup.end()) {
        good_vertex = false;
        bad_vertices.insert(try_v);
        break;
      }
    }
    if (good_vertex) {
      next_slice_lookup.insert(try_v);
      next_slice->push_back(try_v);
    }
  }

  return {
      next_slice, get_next_u_frontier(*this, u_frontier, next_slice_lookup),
      get_next_b_frontier(*this, b_frontier, u_frontier, next_slice_lookup)};
}

CutFrontier Circuit::next_cut(
    std::shared_ptr<const unit_frontier_t> u_frontier,
    std::shared_ptr<const b_frontier_t> b_frontier,
    const std::function<bool(Op_ptr)>& skip_func) const {
  VertexSet bad_vertices;
  std::list<Edge> all_edges;
  EdgeSet edge_lookup;
  for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
    Edge e = pair.second;
    all_edges.push_back(e);
    edge_lookup.insert(e);
  }
  for (const std::pair<Bit, EdgeVec>& pair : b_frontier->get<TagKey>()) {
    for (const Edge& edge : pair.second) {
      Edge e = edge;
      all_edges.push_back(e);
      edge_lookup.insert(e);
    }
  }
  // advance through skippable
  bool can_skip;
  do {
    can_skip = false;
    VertexSet skip_slice_lookup;
    for (const Edge& e : all_edges) {
      Vertex try_v = target(e);
      if (detect_final_Op(try_v) || (!skip_func(get_Op_ptr_from_Vertex(try_v))))
        continue;
      if (skip_slice_lookup.find(try_v) != skip_slice_lookup.end()) continue;
      bool good_vertex = bad_vertices.find(try_v) == bad_vertices.end();
      if (!good_vertex) continue;
      const EdgeVec ins = get_in_edges(try_v);
      for (const Edge& in : ins) {
        if (edge_lookup.find(in) == edge_lookup.end()) {
          good_vertex = false;
          break;
        }
      }
      if (!good_vertex) {
        bad_vertices.insert(try_v);
        continue;
      }
      skip_slice_lookup.insert(try_v);
    }
    if (!skip_slice_lookup.empty()) {
      b_frontier =
          get_next_b_frontier(*this, b_frontier, u_frontier, skip_slice_lookup);
      u_frontier = get_next_u_frontier(*this, u_frontier, skip_slice_lookup);
      bad_vertices = {};
      all_edges = {};
      edge_lookup = {};

      for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
        Edge e = pair.second;
        all_edges.push_back(e);
        edge_lookup.insert(e);
      }
      for (const std::pair<Bit, EdgeVec>& pair : b_frontier->get<TagKey>()) {
        for (const Edge& edge : pair.second) {
          Edge e = edge;
          all_edges.push_back(e);
          edge_lookup.insert(e);
        }
      }
      can_skip = true;
    }
  } while (can_skip);

  // find the next slice first
  auto next_slice = std::make_shared<Slice>();
  VertexSet next_slice_lookup;

  for (const Edge& e : all_edges) {
    Vertex try_v = target(e);
    if (detect_final_Op(try_v)) continue;
    if (next_slice_lookup.find(try_v) != next_slice_lookup.end())
      continue;  // already going to be in next slice
    bool good_vertex = bad_vertices.find(try_v) == bad_vertices.end();
    if (!good_vertex) continue;
    const EdgeVec ins = get_in_edges(try_v);
    for (const Edge& in : ins) {
      if (edge_lookup.find(in) == edge_lookup.end()) {
        good_vertex = false;
        bad_vertices.insert(try_v);
        break;
      }
    }
    if (good_vertex) {
      next_slice_lookup.insert(try_v);
      next_slice->push_back(try_v);
    }
  }

  return {
      next_slice, get_next_u_frontier(*this, u_frontier, next_slice_lookup),
      get_next_b_frontier(*this, b_frontier, u_frontier, next_slice_lookup)};
}

CutFrontier Circuit::next_q_cut(
    std::shared_ptr<const unit_frontier_t> u_frontier) const {
  auto next_slice = std::make_shared<Slice>();
  VertexSet next_slice_lookup;
  VertexSet bad_vertices;
  EdgeSet edge_lookup;
  for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
    edge_lookup.insert(pair.second);
  }

  // find the next slice first
  for (const std::pair<UnitID, Edge>& pair : u_frontier->get<TagKey>()) {
    Vertex try_v = target(pair.second);
    if (detect_final_Op(try_v)) continue;
    if (next_slice_lookup.contains(try_v))
      continue;  // already going to be in next slice
    bool good_vertex = !bad_vertices.contains(try_v);
    if (!good_vertex) continue;
    EdgeVec ins = get_in_edges(try_v);
    for (const Edge& in : ins) {
      if (!edge_lookup.contains(in) && get_edgetype(in) == EdgeType::Quantum) {
        good_vertex = false;
        bad_vertices.insert(try_v);
        break;
      }
    }
    if (good_vertex) {
      next_slice_lookup.insert(try_v);
      next_slice->push_back(try_v);
    }
  }

  return {
      next_slice, get_next_u_frontier(*this, u_frontier, next_slice_lookup),
      std::make_shared<b_frontier_t>()};
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
      case OpType::Discard:
      case OpType::ClInput:
      case OpType::ClOutput: {
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
  Circuit::SliceIterator slice_iter(*this, skip_func);
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
  Circuit::SliceIterator slice_iter(*this, skip_func);
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
  Circuit::SliceIterator slice_iter(*this, skip_func);
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

/*SliceIterator related methods*/
Circuit::SliceIterator::SliceIterator(const Circuit& circ)
    : cut_(), circ_(&circ) {
  cut_.init();
  for (const Qubit& q : circ.all_qubits()) {
    Vertex in = circ.get_in(q);
    cut_.slice->push_back(in);
    cut_.u_frontier->insert({q, circ.get_nth_out_edge(in, 0)});
  }
  for (const Bit& b : circ.all_bits()) {
    Vertex in = circ.get_in(b);
    cut_.slice->push_back(in);
    cut_.b_frontier->insert({b, circ.get_nth_b_out_bundle(in, 0)});
    cut_.u_frontier->insert({b, circ.get_nth_out_edge(in, 0)});
  }
  prev_b_frontier_ = cut_.b_frontier;
  cut_ = circ.next_cut(cut_.u_frontier, cut_.b_frontier);
}

Circuit::SliceIterator::SliceIterator(
    const Circuit& circ, const std::function<bool(Op_ptr)>& skip_func)
    : cut_(), circ_(&circ) {
  cut_.init();
  for (const Qubit& q : circ.all_qubits()) {
    Vertex in = circ.get_in(q);
    cut_.u_frontier->insert({q, circ.get_nth_out_edge(in, 0)});
  }
  for (const Bit& b : circ.all_bits()) {
    Vertex in = circ.get_in(b);
    cut_.b_frontier->insert({b, circ.get_nth_b_out_bundle(in, 0)});
    cut_.u_frontier->insert({b, circ.get_nth_out_edge(in, 0)});
  }
  prev_b_frontier_ = cut_.b_frontier;
  cut_ = circ.next_cut(cut_.u_frontier, cut_.b_frontier, skip_func);
}

Circuit::SliceIterator Circuit::slice_begin() const {
  return Circuit::SliceIterator(*this);
}

Circuit::SliceIterator Circuit::slice_end() { return nullsit; }
const Circuit::SliceIterator Circuit::nullsit = Circuit::SliceIterator();

Circuit::SliceIterator::Sliceholder Circuit::SliceIterator::operator++(int) {
  Sliceholder ret(*cut_.slice);
  ++*this;
  return ret;
}

Circuit::SliceIterator& Circuit::SliceIterator::operator++() {
  if (this->finished()) {
    *this = circ_->slice_end();
    return *this;
  }
  prev_b_frontier_ = cut_.b_frontier;
  cut_ = circ_->next_cut(cut_.u_frontier, cut_.b_frontier);
  return *this;
}

bool Circuit::SliceIterator::finished() const {
  for (const std::pair<UnitID, Edge>& pair :
       this->cut_.u_frontier->get<TagKey>()) {
    if (!circ_->detect_final_Op(circ_->target(pair.second))) return false;
  }
  for (const std::pair<Bit, EdgeVec>& pair :
       this->cut_.b_frontier->get<TagKey>()) {
    if (!pair.second.empty()) return false;
  }
  return true;
}

Circuit::CommandIterator::CommandIterator(const Circuit& circ)
    : current_slice_iterator_(circ.slice_begin()),
      current_index_(0),
      circ_(&circ) {
  if ((*current_slice_iterator_).size() == 0)
    *this = circ.end();
  else {
    current_vertex_ = (*current_slice_iterator_)[0];
    current_command_ = circ.command_from_vertex(
        current_vertex_, current_slice_iterator_.get_u_frontier(),
        current_slice_iterator_.get_prev_b_frontier());
  }
}

const Circuit::CommandIterator Circuit::begin() const {
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
  } else
    ++current_index_;
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
    if (get_edgetype(ins[p]) == EdgeType::Boolean) {
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
        if (found) break;
      }
      if (!found)
        throw CircuitInvalidity(
            "Vertex edges not found in CRead frontier. Edge: " +
            get_Op_ptr_from_Vertex(source(ins[p]))->get_name() + " -> " +
            get_Op_ptr_from_Vertex(target(ins[p]))->get_name());
    } else {
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
