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

#pragma once

#include <optional>

#include "Gate/Gate.hpp"
#include "Gate/GatePtr.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Expression.hpp"

namespace tket {

/**
 * @brief Squashes single qubit gates using given Squasher
 *
 * The template Squasher must define the following interface:
 *   - bool accepts(OpType): whether the OpType can be added to current squash
 *   - void append(Gate_ptr): add gate to current squash
 *   - std::pair<Circuit, Gate_ptr> flush(std::optional<Pauli>):
 *     return a squashed circuit, optionally use the commutation colour of the
 *     next gate to return an additional Gate_ptr to be commuted through
 *     (otherwise, set Gate_ptr = nullptr)
 *   - void clear(): reset the current squash
 */
template <typename Squasher>
class SingleQubitSquash {
 private:
  using Condition = std::optional<std::pair<std::list<VertPort>, unsigned>>;

 public:
  SingleQubitSquash(const Squasher &_squasher, bool _reversed = false)
      : squasher(_squasher), reversed(_reversed), circ_ptr(nullptr) {}

  // squash a circuit
  // main method to be called by user
  bool squash(Circuit &circ) {
    bool success = false;
    circ_ptr = &circ;

    VertexVec starts(reversed ? circ.q_outputs() : circ.q_inputs());
    for (Vertex v : starts) {
      Edge e(
          reversed ? circ.get_nth_in_edge(v, 0) : circ.get_nth_out_edge(v, 0));
      success |= squash_wire(e);
    }

    circ_ptr = nullptr;
    return success;
  }

  // squash a single wire, starting at e
  bool squash_wire(Edge e) {
    squasher.clear();
    Vertex v = next_vertex(e);
    std::vector<Gate_ptr> single_chain;
    VertexVec bin;
    bool success = false;
    Condition condition = std::nullopt;
    while (true) {
      Op_ptr v_op = circ_ptr->get_Op_ptr_from_Vertex(v);
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

      if (condition == this_condition && is_squashable(v, v_type)) {
        // => add gate to current squash
        squasher.append(as_gate_ptr(reversed ? v_op->dagger() : v_op));
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
            commutation_colour = v_op->commuting_basis(next_port(e));
            move_to_next_vertex = true;
          }
          auto pair = squasher.flush(commutation_colour);
          sub = pair.first;
          Gate_ptr left_over_gate = pair.second;
          if (left_over_gate != nullptr) {
            // => commute leftover through before squashing
            insert_left_over_gate(left_over_gate, next_edge(v, e), condition);
            left_over_gate = nullptr;
          }
          if (reversed) {
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
      if (is_last_optype(v_type)) {
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
        squasher.clear();
        condition = std::nullopt;
      }
    }
    return success;
  }

 private:
  Squasher squasher;
  bool reversed;
  // points to the current circuit during squashing
  Circuit *circ_ptr;

  // substitute chain by a sub circuit, handling conditions
  // and backing up + restoring current edge
  void substitute(
      const Circuit &sub, const VertexVec &single_chain, Edge &e,
      const Condition &condition) {
    // backup edge
    VertPort backup = {next_vertex(e), next_port(e)};

    if (condition) {
      circ_ptr->substitute_conditional(
          sub, single_chain.front(), Circuit::VertexDeletion::No);
    } else {
      circ_ptr->substitute(
          sub, single_chain.front(), Circuit::VertexDeletion::No);
    }
    circ_ptr->remove_vertices(
        VertexSet{single_chain.begin(), single_chain.end()},
        Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

    // restore backup
    e = prev_edge(backup);
  }

  // insert a gate at the given edge, respecting condition
  void insert_left_over_gate(
      Op_ptr left_over, Edge e, const Condition &condition) {
    if (reversed) {
      left_over = left_over->dagger();
    }
    EdgeVec preds;
    op_signature_t sigs;
    if (condition) {
      left_over = std::make_shared<Conditional>(
          left_over, condition->first.size(), condition->second);
    }
    Vertex new_v = circ_ptr->add_vertex(left_over);
    if (condition) {
      for (const VertPort &vp : condition->first) {
        preds.push_back(circ_ptr->get_nth_out_edge(vp.first, vp.second));
        sigs.push_back(EdgeType::Boolean);
      }
    }
    preds.push_back(e);
    sigs.push_back(EdgeType::Quantum);
    circ_ptr->rewire(new_v, preds, sigs);
  }

  // whether a vertex can be squashed with the previous vertices
  bool is_squashable(Vertex v, OpType v_type) {
    return circ_ptr->n_in_edges_of_type(v, EdgeType::Quantum) == 1 &&
           is_gate_type(v_type) && squasher.accepts(v_type);
  }

  // whether the sub circuit is shorter than chain
  bool sub_is_better(const Circuit &sub, const std::vector<Gate_ptr> chain) {
    return sub.n_gates() < chain.size() ||
           (sub.n_gates() == chain.size() && !is_equal(sub, chain, reversed));
  }

  // returns a description of the condition of current vertex
  Condition get_condition(Vertex v) {
    Op_ptr v_op = circ_ptr->get_Op_ptr_from_Vertex(v);
    OpType v_type = v_op->get_type();
    if (v_type != OpType::Conditional) {
      throw NotValid("Cannot get condition from non-conditional OpType");
    }
    const Conditional &cond_op = static_cast<const Conditional &>(*v_op);
    EdgeVec ins = circ_ptr->get_in_edges(v);
    Condition cond = std::pair<std::list<VertPort>, unsigned>();
    for (port_t p = 0; p < cond_op.get_width(); ++p) {
      Edge in_p = ins.at(p);
      VertPort vp = {circ_ptr->source(in_p), circ_ptr->get_source_port(in_p)};
      cond->first.push_back(vp);
    }
    cond->second = cond_op.get_value();
    return cond;
  }

  // simple utils respecting reversed boolean
  Vertex next_vertex(const Edge &e) const {
    return reversed ? circ_ptr->source(e) : circ_ptr->target(e);
  }
  port_t next_port(const Edge &e) const {
    return reversed ? circ_ptr->get_source_port(e)
                    : circ_ptr->get_target_port(e);
  }
  Edge prev_edge(const VertPort &pair) const {
    return reversed ? circ_ptr->get_nth_out_edge(pair.first, pair.second)
                    : circ_ptr->get_nth_in_edge(pair.first, pair.second);
  }
  Edge next_edge(const Vertex &v, const Edge &e) const {
    return reversed ? circ_ptr->get_last_edge(v, e)
                    : circ_ptr->get_next_edge(v, e);
  }
  bool is_last_optype(OpType type) const {
    return (reversed && is_initial_q_type(type)) ||
           (!reversed && is_final_q_type(type));
  }

  // checks whether a 1-qb circuit is equal to a chain of single-qb gates
  // it checks (i) that the sequence of OpTypes is identical
  //          (ii) that the parameters match for each gate
  static bool is_equal(
      const Circuit &circ, const std::vector<Gate_ptr> &gates,
      bool reversed = false) {
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
      const OpType ot1 = op1->get_type();
      const OpType ot2 = op2->get_type();
      if (ot1 != ot2) {
        return false;
      }
      const std::vector<Expr> params1 = op1->get_params();
      const std::vector<Expr> params2 = op2->get_params();
      auto p1 = params1.cbegin();
      auto p2 = params2.cbegin();
      unsigned i = 0;
      while (p1 != params1.cend() && p2 != params2.cend()) {
        unsigned mod = op1->get_desc().param_mod(i);
        if (!equiv_expr(*p1, *p2, mod)) {
          return false;
        }
        ++p1;
        ++p2;
        ++i;
      }
      if (p1 != params1.cend() || p2 != params2.cend()) {
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
};

}  // namespace tket
