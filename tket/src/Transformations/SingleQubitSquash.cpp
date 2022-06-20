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

#include "SingleQubitSquash.hpp"

#include "Circuit/DAGDefs.hpp"
#include "Gate/Gate.hpp"

namespace tket {

SingleQubitSquash::SingleQubitSquash(const SingleQubitSquash &other)
    : squasher_(other.squasher_->clone()),
      circ_(other.circ_),
      reversed_(other.reversed_) {}

SingleQubitSquash &SingleQubitSquash::operator=(
    const SingleQubitSquash &other) {
  squasher_ = other.squasher_->clone();
  circ_ = other.circ_;
  reversed_ = other.reversed_;
  return *this;
}

SingleQubitSquash::SingleQubitSquash(SingleQubitSquash &&other)
    : squasher_(std::move(other.squasher_)),
      circ_(other.circ_),
      reversed_(other.reversed_) {}

SingleQubitSquash &SingleQubitSquash::operator=(SingleQubitSquash &&other) {
  squasher_ = std::move(other.squasher_);
  circ_ = other.circ_;
  reversed_ = other.reversed_;
  return *this;
}

bool SingleQubitSquash::squash() {
  bool success = false;

  VertexVec inputs = circ_.q_inputs();
  VertexVec outputs = circ_.q_outputs();
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    Edge in = circ_.get_nth_out_edge(inputs[i], 0);
    Edge out = circ_.get_nth_in_edge(outputs[i], 0);
    if (reversed_) {
      success |= squash_between(out, in);
    } else {
      success |= squash_between(in, out);
    }
  }

  return success;
}

bool SingleQubitSquash::squash_between(const Edge &in, const Edge &out) {
  squasher_->clear();
  Edge e = in;
  Vertex v = next_vertex(e);
  std::vector<Gate_ptr> single_chain;
  VertexVec bin;
  bool success = false;
  Condition condition = std::nullopt;
  while (true) {
    Op_ptr v_op = circ_.get_Op_ptr_from_Vertex(v);
    OpType v_type = v_op->get_type();
    bool move_to_next_vertex = false;
    bool reset_search = false;
    Condition this_condition = std::nullopt;

    if (v_type == OpType::Conditional) {
      // => deal with conditional case
      this_condition = get_condition(v);
      v_op = static_cast<const Conditional &>(*v_op).get_op();
      v_type = v_op->get_type();

      if (single_chain.empty()) {
        condition = this_condition;
      }
    }

    bool is_squashable = circ_.n_in_edges_of_type(v, EdgeType::Quantum) == 1 &&
                         is_gate_type(v_type) &&
                         squasher_->accepts(as_gate_ptr(v_op));

    if (e != out && condition == this_condition && is_squashable) {
      // => add gate to current squash
      squasher_->append(as_gate_ptr(reversed_ ? v_op->dagger() : v_op));
      move_to_next_vertex = true;
    } else {
      // => squash and reset
      reset_search = true;
      if (single_chain.empty()) {
        // => nothing to do, move on
        move_to_next_vertex = true;
      } else {
        Circuit sub;
        std::optional<Pauli> commutation_colour = std::nullopt;
        if (is_gate_type(v_type) && v_op->n_qubits() > 1) {
          commutation_colour =
              circ_.commuting_basis(v, PortType::Target, next_port(e));
          move_to_next_vertex = true;
        }
        auto pair = squasher_->flush(commutation_colour);
        sub = pair.first;
        Gate_ptr left_over_gate = pair.second;
        if (left_over_gate != nullptr) {
          // => commute leftover through before squashing
          insert_left_over_gate(left_over_gate, next_edge(v, e), condition);
          left_over_gate = nullptr;
        }
        if (reversed_) {
          sub = sub.dagger();
        }

        // we squash if the replacement is at least as good as the original
        // (and it's not a no-op)
        if (sub_is_better(sub, single_chain)) {
          substitute(sub, bin, e, condition);
          success = true;
        }
      }
    }
    if (e == out || is_last_optype(v_type)) {
      squasher_->clear();
      break;
    }
    if (move_to_next_vertex) {
      if (is_gate_type(v_type)) {
        bin.push_back(v);
        single_chain.push_back(as_gate_ptr(v_op));
      }
      e = next_edge(v, e);
      v = next_vertex(e);
    }
    if (reset_search) {
      bin.clear();
      single_chain.clear();
      squasher_->clear();
      condition = std::nullopt;
    }
  }
  return success;
}

void SingleQubitSquash::substitute(
    const Circuit &sub, const VertexVec &single_chain, Edge &e,
    const Condition &condition) {
  // backup edge
  VertPort backup = {next_vertex(e), next_port(e)};

  if (condition) {
    circ_.substitute_conditional(
        sub, single_chain.front(), Circuit::VertexDeletion::No);
  } else {
    circ_.substitute(sub, single_chain.front(), Circuit::VertexDeletion::No);
  }
  circ_.remove_vertices(
      VertexSet{single_chain.begin(), single_chain.end()},
      Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

  // restore backup
  e = prev_edge(backup);
}

void SingleQubitSquash::insert_left_over_gate(
    Op_ptr left_over, const Edge &e, const Condition &condition) {
  if (reversed_) {
    left_over = left_over->dagger();
  }
  EdgeVec preds;
  op_signature_t sigs;
  if (condition) {
    left_over = std::make_shared<Conditional>(
        left_over, (unsigned)condition->first.size(), condition->second);
  }
  Vertex new_v = circ_.add_vertex(left_over);
  if (condition) {
    for (const VertPort &vp : condition->first) {
      preds.push_back(circ_.get_nth_out_edge(vp.first, vp.second));
      sigs.push_back(EdgeType::Boolean);
    }
  }
  preds.push_back(e);
  sigs.push_back(EdgeType::Quantum);
  circ_.rewire(new_v, preds, sigs);
}

bool SingleQubitSquash::sub_is_better(
    const Circuit &sub, const std::vector<Gate_ptr> chain) const {
  const unsigned n_gates = sub.n_gates();
  return n_gates < chain.size() ||
         (n_gates == chain.size() && !is_equal(sub, chain, reversed_));
}

// returns a description of the condition of current vertex
SingleQubitSquash::Condition SingleQubitSquash::get_condition(Vertex v) const {
  Op_ptr v_op = circ_.get_Op_ptr_from_Vertex(v);
  OpType v_type = v_op->get_type();
  if (v_type != OpType::Conditional) {
    throw NotValid("Cannot get condition from non-conditional OpType");
  }
  const Conditional &cond_op = static_cast<const Conditional &>(*v_op);
  EdgeVec ins = circ_.get_in_edges(v);
  Condition cond = std::pair<std::list<VertPort>, unsigned>();
  for (port_t p = 0; p < cond_op.get_width(); ++p) {
    Edge in_p = ins.at(p);
    VertPort vp = {circ_.source(in_p), circ_.get_source_port(in_p)};
    cond->first.push_back(vp);
  }
  cond->second = cond_op.get_value();
  return cond;
}

// simple utils respecting reversed boolean
Vertex SingleQubitSquash::next_vertex(const Edge &e) const {
  return reversed_ ? circ_.source(e) : circ_.target(e);
}

port_t SingleQubitSquash::next_port(const Edge &e) const {
  return reversed_ ? circ_.get_source_port(e) : circ_.get_target_port(e);
}

Edge SingleQubitSquash::prev_edge(const VertPort &pair) const {
  return reversed_ ? circ_.get_nth_out_edge(pair.first, pair.second)
                   : circ_.get_nth_in_edge(pair.first, pair.second);
}

Edge SingleQubitSquash::next_edge(const Vertex &v, const Edge &e) const {
  return reversed_ ? circ_.get_last_edge(v, e) : circ_.get_next_edge(v, e);
}

bool SingleQubitSquash::is_last_optype(OpType type) const {
  return (reversed_ && is_initial_q_type(type)) ||
         (!reversed_ && is_final_q_type(type));
}

bool SingleQubitSquash::is_equal(
    const Circuit &circ, const std::vector<Gate_ptr> &gates, bool reversed) {
  if (reversed) {
    return is_equal(circ, {gates.rbegin(), gates.rend()});
  }
  if (circ.n_qubits() != 1) {
    throw NotValid("Only circuits with one qubit are supported");
  }

  auto it1 = circ.begin();
  auto it2 = gates.cbegin();

  while (it1 != circ.end() && it2 != gates.end()) {
    const Gate_ptr op1 = as_gate_ptr(it1->get_op_ptr());
    const Gate_ptr op2 = as_gate_ptr(*it2);
    if (!(*op1 == *op2)) {
      return false;
    }
    ++it1;
    ++it2;
  }

  if (it1 != circ.end() || it2 != gates.cend()) {
    return false;
  }
  return true;
}

}  // namespace tket
