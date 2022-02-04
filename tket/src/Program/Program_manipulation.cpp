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

#include "Program.hpp"

namespace tket {

Program::Program() {
  entry_ = add_vertex(Circuit());
  exit_ = add_vertex(Circuit());
  add_edge(entry_, exit_);
}

Program::Program(const Program &to_copy) {
  std::map<FGVert, FGVert> isomap = copy_graph(to_copy);
  entry_ = isomap.at(to_copy.entry_);
  exit_ = isomap.at(to_copy.exit_);
}

Program::Program(unsigned qubits, unsigned bits) : Program() {
  add_q_register(q_default_reg(), qubits);
  add_c_register(c_default_reg(), bits);
}

FGVert Program::add_vertex(
    const Circuit &circ, std::optional<Bit> branch_condition,
    const std::optional<std::string> &label) {
  for (const Qubit &qb : circ.all_qubits()) {
    add_qubit(qb, false);
  }
  for (const Bit &b : circ.all_bits()) {
    add_bit(b, false);
  }
  FGVert new_v = boost::add_vertex(flow_);
  flow_[new_v] = {circ, branch_condition, label};
  return new_v;
}

void Program::remove_vertex(const FGVert &vert) {
  boost::clear_vertex(vert, flow_);
  boost::remove_vertex(vert, flow_);
}

FGEdge Program::add_edge(
    const FGVert &source, const FGVert &target, bool branch) {
  std::pair<FGEdge, bool> edge_pair = boost::add_edge(source, target, flow_);
  if (!edge_pair.second) throw ProgramError("Could not add edge to flow graph");
  FGEdge new_e = edge_pair.first;
  flow_[new_e] = {branch};
  return new_e;
}

void Program::remove_edge(const FGEdge &edge) {
  boost::remove_edge(edge, flow_);
}

std::map<FGVert, FGVert> Program::copy_graph(const Program &to_copy) {
  std::map<FGVert, FGVert> isomap;
  if (&to_copy == this) {
    throw ProgramError("Cannot copy a program into itself");
  }
  for (const Qubit &qb : to_copy.all_qubits()) {
    add_qubit(qb, false);
  }
  for (const Bit &b : to_copy.all_bits()) {
    add_bit(b, false);
  }
  BGL_FORALL_VERTICES(v, to_copy.flow_, FlowGraph) {
    FGVert new_v = boost::add_vertex(this->flow_);
    this->flow_[new_v] = to_copy.flow_[v];
    isomap.insert({v, new_v});
  }
  BGL_FORALL_EDGES(e, to_copy.flow_, FlowGraph) {
    FGVert source = isomap.at(to_copy.get_source(e));
    FGVert target = isomap.at(to_copy.get_target(e));
    bool branch = to_copy.get_branch(e);
    add_edge(source, target, branch);
  }
  return isomap;
}

FGVert Program::add_block(const Circuit &circ) {
  FGVert block = add_vertex(circ);
  FGEdgeVec ins = get_in_edges(exit_);
  for (const FGEdge &e : ins) {
    add_edge(get_source(e), block, get_branch(e));
    remove_edge(e);
  }
  add_edge(block, exit_);
  return block;
}

void Program::add_op(const Op_ptr &op, const std::vector<unsigned> &args) {
  unit_vector_t arg_ids;
  op_signature_t sig = op->get_signature();
  for (unsigned i = 0; i < args.size(); ++i) {
    if (sig.at(i) == EdgeType::Quantum) {
      arg_ids.push_back(Qubit(args[i]));
    } else {
      arg_ids.push_back(Bit(args[i]));
    }
  }
  return add_op(op, arg_ids);
}

void Program::add_op(const Op_ptr &op, const unit_vector_t &args) {
  FGVertVec lasts = get_predecessors(exit_);
  FGVert block;
  if (lasts.size() == 1 && lasts.front() != entry_ &&
      !get_condition(lasts.front())) {
    block = lasts.front();
  } else {
    block = add_block({});
  }
  Circuit &circ = flow_[block].circ;
  op_signature_t sig = op->get_signature();
  for (unsigned i = 0; i < args.size(); ++i) {
    if (sig.at(i) == EdgeType::Quantum) {
      circ.add_qubit(Qubit(args[i]), false);
    } else {
      circ.add_bit(Bit(args[i]), false);
    }
  }
  circ.add_op(op, args);
}

void Program::append(const Program &to_append) {
  std::map<FGVert, FGVert> isomap = copy_graph(to_append);
  FGEdgeVec ins = get_in_edges(exit_);
  FGVert added_entry = isomap.at(to_append.entry_);
  FGVert target = get_branch_successor(added_entry, false);
  for (const FGEdge &e : ins) {
    FGVert source = get_source(e);
    bool branch = get_branch(e);
    add_edge(source, target, branch);
  }
  remove_vertex(added_entry);
  remove_vertex(exit_);
  exit_ = isomap.at(to_append.exit_);
}

void Program::append_if(const Bit &condition_bit, const Program &body) {
  std::map<FGVert, FGVert> isomap = copy_graph(body);
  FGVert added_entry = isomap.at(body.entry_);
  FGVert added_exit = isomap.at(body.exit_);
  FGVert target = get_branch_successor(added_entry, false);
  flow_[exit_].branch_condition = condition_bit;
  add_edge(exit_, target, true);
  add_edge(exit_, added_exit, false);
  remove_vertex(added_entry);
  exit_ = added_exit;
}

void Program::append_if_else(
    const Bit &condition_bit, const Program &if_body,
    const Program &else_body) {
  std::map<FGVert, FGVert> if_map = copy_graph(if_body);
  FGVert if_entry = if_map.at(if_body.entry_);
  FGVert if_exit = if_map.at(if_body.exit_);
  FGVert if_target = get_branch_successor(if_entry, false);
  std::map<FGVert, FGVert> else_map = copy_graph(else_body);
  FGVert else_entry = else_map.at(else_body.entry_);
  FGVert else_exit = else_map.at(else_body.exit_);
  FGVert else_target = get_branch_successor(else_entry, false);
  flow_[exit_].branch_condition = condition_bit;
  add_edge(exit_, if_target, true);
  add_edge(exit_, else_target, false);
  remove_vertex(if_entry);
  remove_vertex(else_entry);
  add_edge(if_exit, else_exit, false);
  exit_ = else_exit;
}

void Program::append_while(const Bit &condition_bit, const Program &body) {
  std::map<FGVert, FGVert> isomap = copy_graph(body);
  FGVert added_entry = isomap.at(body.entry_);
  FGVert added_exit = isomap.at(body.exit_);
  FGVert target = get_branch_successor(added_entry, false);
  FGVert new_exit = add_vertex(Circuit());
  flow_[added_exit].branch_condition = condition_bit;
  add_edge(added_exit, target, true);
  add_edge(added_exit, new_exit, false);
  add_edge(exit_, added_exit, false);
  remove_vertex(added_entry);
  exit_ = new_exit;
}

Program operator>>(const Program &p1, const Program &p2) {
  Program new_prog = p1;
  new_prog.append(p2);
  return new_prog;
}

}  // namespace tket
