// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <algorithm>

#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/PauliGraph/PauliGraph.hpp"
#include "tket/Transformations/CliffordOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisationLookupTables.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

// convert a Clifford tableau to a vector of PauliNode_ptr
static std::vector<PauliNode_ptr> get_nodes_from_tableau(
    const UnitaryRevTableau& tab, unsigned n_qubits) {
  std::vector<PauliNode_ptr> rows;
  for (unsigned i = 0; i < n_qubits; i++) {
    Qubit q(i);
    SpPauliStabiliser z_stab = tab.get_zrow(q);
    SpPauliStabiliser x_stab = tab.get_xrow(q);
    bool z_sign = cast_coeff<quarter_turns_t, Complex>(z_stab.coeff) == 1.;
    bool x_sign = cast_coeff<quarter_turns_t, Complex>(x_stab.coeff) == 1.;
    TKET_ASSERT(z_stab.string.size() == n_qubits);
    std::vector<Pauli> z_string;
    std::vector<Pauli> x_string;
    for (unsigned j = 0; j < n_qubits; j++) {
      z_string.push_back(z_stab.string.at(Qubit(j)));
      x_string.push_back(x_stab.string.at(Qubit(j)));
    }
    rows.push_back(std::make_shared<PauliPropagation>(
        z_string, x_string, z_sign, x_sign, i));
  }
  return rows;
}

std::tuple<std::vector<PauliNode_ptr>, std::vector<PauliNode_ptr>>
gpg_from_unordered_set(const std::vector<SymPauliTensor>& unordered_set) {
  std::vector<PauliNode_ptr> rotation_set;
  unsigned n_qubits = unordered_set.at(0).string.size();
  for (auto& pauli : unordered_set) {
    TKET_ASSERT(pauli.string.size() == n_qubits);
    rotation_set.push_back(
        std::make_shared<PauliRotation>(pauli.string, pauli.coeff));
  }
  UnitaryRevTableau tab(n_qubits);
  std::vector<PauliNode_ptr> rows = get_nodes_from_tableau(tab, n_qubits);
  return {rotation_set, rows};
}

// given a stabiliser Pauli string and an angle, return a dense string and an
// angle
static std::pair<std::vector<Pauli>, Expr> dense_pauli(
    const SpPauliStabiliser& pauli, const unsigned& n_qubits,
    const Expr& angle) {
  bool sign = cast_coeff<quarter_turns_t, Complex>(pauli.coeff) == 1.;
  std::vector<Pauli> string(n_qubits, Pauli::I);
  for (const auto& pair : pauli.string) {
    string[pair.first.index().at(0)] = pair.second;
  }
  return {string, sign ? angle : -angle};
}

static bool strings_commute(
    const std::vector<Pauli>& s1, const std::vector<Pauli>& s2) {
  unsigned n_conflicts = 0;
  TKET_ASSERT(s1.size() == s2.size());
  for (unsigned i = 0; i < s1.size(); ++i) {
    Pauli p = s1[i];
    Pauli p2 = s2[i];
    if (p != Pauli::I && p2 != Pauli::I && p != p2) n_conflicts++;
  }
  return (n_conflicts % 2) == 0;
}

static bool nodes_commute(const PauliNode_ptr& n1, const PauliNode_ptr& n2) {
  CommuteInfo c1 = n1->get_commute_info();
  CommuteInfo c2 = n2->get_commute_info();
  // check if every string in n1 commutes with all strings in n2
  for (const std::vector<Pauli>& p1 : c1.paulis) {
    for (const std::vector<Pauli>& p2 : c2.paulis) {
      if (!strings_commute(p1, p2)) return false;
    }
  }
  // check if the bits commute
  for (const std::pair<UnitID, BitType>& b1 : c1.bits_info) {
    for (const std::pair<UnitID, BitType>& b2 : c2.bits_info) {
      if (b1.first == b2.first) {
        // if two nodes read the same bit it's OK
        if (b1.second == BitType::READ && b2.second == BitType::READ) {
          break;
        }
        return false;
      }
    }
  }
  return true;
}

GPGraph::GPGraph(const Circuit& circ)
    : n_qubits_(circ.n_qubits()), n_bits_(circ.n_bits()) {
  qubit_vector_t qubits = circ.all_qubits();
  bit_vector_t bits = circ.all_bits();
  for (const Qubit& q : qubits) {
    TKET_ASSERT(q.reg_name() == q_default_reg());
    TKET_ASSERT(q.index().at(0) < qubits.size());
  }
  for (const Bit& b : bits) {
    TKET_ASSERT(b.reg_name() == c_default_reg());
    TKET_ASSERT(b.index().at(0) < bits.size());
  }
  cliff_ = UnitaryRevTableau(n_qubits_);
  for (const Command& cmd : circ.get_commands()) {
    apply_gate_at_end(cmd);
  }
}

GPVertSet GPGraph::get_successors(const GPVert& vert) const {
  GPVertSet succs;
  for (auto iter = boost::adjacent_vertices(vert, graph_);
       iter.first != iter.second; iter.first++) {
    succs.insert(*iter.first);
  }
  return succs;
}

GPVertSet GPGraph::get_predecessors(const GPVert& vert) const {
  GPVertSet preds;
  for (auto iter = boost::inv_adjacent_vertices(vert, graph_);
       iter.first != iter.second; iter.first++) {
    preds.insert(*iter.first);
  }
  return preds;
}

void GPGraph::apply_node_at_end(PauliNode_ptr& node) {
  GPVertSet to_search = end_line_;
  GPVertSet commuted;
  GPVert new_vert = boost::add_vertex(graph_);
  graph_[new_vert] = node;
  while (!to_search.empty()) {
    // Get next candidate parent
    GPVert to_compare = *to_search.begin();
    to_search.erase(to_search.begin());
    // Check that we have already commuted past all of its children
    bool ready = true;
    for (const PauliVert& child : get_successors(to_compare)) {
      if (commuted.get<TagKey>().find(child) == commuted.get<TagKey>().end()) {
        ready = false;
        break;
      }
    }
    if (!ready) continue;
    // Check if we can commute past it
    PauliNode_ptr compare_node = graph_[to_compare];
    if (nodes_commute(node, compare_node)) {
      // Check if two pauli rotations can be merged
      if (node->get_type() == PauliNodeType::Rotation &&
          compare_node->get_type() == PauliNodeType::Rotation) {
        const PauliRotation& rot1 = dynamic_cast<const PauliRotation&>(*node);
        const PauliRotation& rot2 =
            dynamic_cast<const PauliRotation&>(*compare_node);
        if (rot1.string() == rot2.string()) {
          boost::clear_vertex(new_vert, graph_);
          boost::remove_vertex(new_vert, graph_);
          Expr merged_angle = rot1.angle() + rot2.angle();
          std::optional<unsigned> cl_ang = equiv_Clifford(merged_angle);
          if (cl_ang) {
            cliff_.apply_pauli_at_front(
                SpPauliStabiliser(rot1.string()), *cl_ang);
            start_line_.erase(to_compare);
            for (const GPVert& v : get_predecessors(to_compare)) {
              if (boost::out_degree(v, graph_) == 1) {
                end_line_.insert(v);
              }
            }
            end_line_.erase(to_compare);
            boost::clear_vertex(to_compare, graph_);
            boost::remove_vertex(to_compare, graph_);
          } else {
            graph_[to_compare] =
                std::make_shared<PauliRotation>(rot1.string(), merged_angle);
          }
          return;
        }
      }
      // Commute and continue searching
      PauliVertSet preds = get_predecessors(to_compare);
      to_search.insert(preds.begin(), preds.end());
      commuted.insert(to_compare);
    } else {
      // Does not commute - add dependency edge
      boost::add_edge(to_compare, new_vert, graph_);
      end_line_.erase(to_compare);
    }
  }
  end_line_.insert(new_vert);
  if (get_predecessors(new_vert).empty()) start_line_.insert(new_vert);
}

void GPGraph::apply_pauli_at_end(
    const std::vector<Pauli>& paulis, const Expr& angle,
    const qubit_vector_t& qbs, bool conditional,
    std::vector<unsigned> cond_bits, unsigned cond_value) {
  // Note that global phase is ignored
  if (static_cast<std::size_t>(std::count(
          paulis.begin(), paulis.end(), Pauli::I)) == paulis.size()) {
    return;
  }
  std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
  if (cliff_angle && cliff_angle.value() == 0) {
    return;
  }
  QubitPauliMap qpm;
  for (unsigned i = 0; i != qbs.size(); ++i)
    qpm.insert({Qubit(qbs[i]), paulis[i]});
  if (cliff_angle && !conditional) {
    cliff_.apply_pauli_at_end(SpPauliStabiliser(qpm), *cliff_angle);
    return;
  }
  SpPauliStabiliser qpt = cliff_.get_row_product(SpPauliStabiliser(qpm));
  auto [pauli_dense, theta] = dense_pauli(qpt, n_qubits_, angle);
  PauliNode_ptr node;
  if (conditional) {
    node = std::make_shared<ConditionalPauliRotation>(
        pauli_dense, theta, cond_bits, cond_value);
  } else {
    node = std::make_shared<PauliRotation>(pauli_dense, theta);
  }
  apply_node_at_end(node);
}

void GPGraph::apply_gate_at_end(
    const Command& cmd, bool conditional, std::vector<unsigned> cond_bits,
    unsigned cond_value) {
  const Op_ptr op = cmd.get_op_ptr();
  unit_vector_t args = cmd.get_args();
  qubit_vector_t qbs = cmd.get_qubits();
  OpType type = op->get_type();

  for (const UnitID& arg : args) {
    if (arg.type() == UnitType::Qubit) {
      if (measures_.left.find(arg.index().at(0)) != measures_.left.end()) {
        throw MidCircuitMeasurementNotAllowed(
            "PauliGraph does not support mid-circuit measurements "
            "- cannot add gate after measure on qubit " +
            arg.repr());
      }
    } else if (
        measures_.right.find(arg.index().at(0)) != measures_.right.end()) {
      throw MidCircuitMeasurementNotAllowed(
          "PauliGraph does not support mid-circuit measurements - "
          "cannot add gate after measure to bit " +
          arg.repr());
    }
  }

  std::vector<std::pair<std::vector<Pauli>, Expr>> pauli_rots;
  switch (type) {
    case OpType::Conditional: {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      std::vector<unsigned> cond_bits;
      unit_vector_t inner_args;
      for (unsigned i = 0; i < cond.get_width(); ++i)
        cond_bits.push_back(Bit(args.at(i)).index().at(0));
      for (unsigned i = cond.get_width(); i < args.size(); ++i)
        inner_args.push_back(args.at(i));
      apply_gate_at_end(
          Command(cond.get_op(), inner_args), true, cond_bits,
          cond.get_value());
      return;
    }
    case OpType::Measure: {
      measures_.insert({args.at(0).index().at(0), args.at(1).index().at(0)});
      return;
    }
    case OpType::Z: {
      pauli_rots.push_back({{Pauli::Z}, 1});
      break;
    }
    case OpType::X: {
      pauli_rots.push_back({{Pauli::X}, 1});
      break;
    }
    case OpType::Y: {
      pauli_rots.push_back({{Pauli::Y}, 1});
      break;
    }
    case OpType::S: {
      pauli_rots.push_back({{Pauli::Z}, 0.5});
      break;
    }
    case OpType::V: {
      pauli_rots.push_back({{Pauli::X}, 0.5});
      break;
    }
    case OpType::Sdg: {
      pauli_rots.push_back({{Pauli::Z}, 1.5});
      break;
    }
    case OpType::Vdg: {
      pauli_rots.push_back({{Pauli::X}, 1.5});
      break;
    }
    case OpType::H: {
      pauli_rots.push_back({{Pauli::Y}, 0.5});
      pauli_rots.push_back({{Pauli::X}, 1});
      break;
    }
    case OpType::CX:
    case OpType::CY:
    case OpType::CZ: {
      Pauli t = (type == OpType::CZ)   ? Pauli::Z
                : (type == OpType::CX) ? Pauli::X
                                       : Pauli::Y;
      pauli_rots.push_back({{Pauli::Z, Pauli::I}, 1.5});
      pauli_rots.push_back({{Pauli::I, t}, 1.5});
      pauli_rots.push_back({{Pauli::Z, t}, 0.5});
      break;
    }
    case OpType::SWAP: {
      pauli_rots.push_back({{Pauli::Z, Pauli::Z}, 0.5});
      pauli_rots.push_back({{Pauli::X, Pauli::X}, 0.5});
      pauli_rots.push_back({{Pauli::Y, Pauli::Y}, 0.5});
      break;
    }
    case OpType::noop:
    case OpType::Phase: {
      // ignore global phase
      return;
    }
    case OpType::Rz: {
      pauli_rots.push_back({{Pauli::Z}, op->get_params().at(0)});
      break;
    }
    case OpType::Rx: {
      pauli_rots.push_back({{Pauli::X}, op->get_params().at(0)});
      break;
    }
    case OpType::Ry: {
      pauli_rots.push_back({{Pauli::Y}, op->get_params().at(0)});
      break;
    }
    case OpType::PhasedX: {
      Expr alpha = op->get_params().at(0);
      Expr beta = op->get_params().at(1);
      pauli_rots.push_back({{Pauli::Z}, -beta});
      pauli_rots.push_back({{Pauli::X}, alpha});
      pauli_rots.push_back({{Pauli::Z}, beta});
      break;
    }
    case OpType::T: {
      pauli_rots.push_back({{Pauli::Z}, 0.25});
      break;
    }
    case OpType::Tdg: {
      pauli_rots.push_back({{Pauli::Z}, -0.25});
      break;
    }
    case OpType::ZZMax: {
      pauli_rots.push_back({{Pauli::Z, Pauli::Z}, 0.5});
      break;
    }
    case OpType::PhaseGadget:
    case OpType::ZZPhase: {
      Expr angle = op->get_params().at(0);
      std::vector<Pauli> paulis(qbs.size(), Pauli::Z);
      pauli_rots.push_back({paulis, angle});
      break;
    }
    case OpType::XXPhase: {
      Expr angle = op->get_params().at(0);
      pauli_rots.push_back({{Pauli::X, Pauli::X}, angle});
      break;
    }
    case OpType::YYPhase: {
      Expr angle = op->get_params().at(0);
      pauli_rots.push_back({{Pauli::Y, Pauli::Y}, angle});
      break;
    }
    case OpType::PauliExpBox: {
      const PauliExpBox& peb = static_cast<const PauliExpBox&>(*op);
      pauli_rots.push_back({peb.get_paulis(), peb.get_phase()});
      break;
    }
    case OpType::PauliExpPairBox: {
      const PauliExpPairBox& peb = static_cast<const PauliExpPairBox&>(*op);
      auto [paulis1, paulis2] = peb.get_paulis_pair();
      auto [phase1, phase2] = peb.get_phase_pair();
      pauli_rots.push_back({paulis1, phase1});
      pauli_rots.push_back({paulis2, phase2});
      break;
    }
    case OpType::PauliExpCommutingSetBox: {
      const PauliExpCommutingSetBox& peb =
          static_cast<const PauliExpCommutingSetBox&>(*op);
      for (const SymPauliTensor& pt : peb.get_pauli_gadgets()) {
        pauli_rots.push_back({pt.string, pt.coeff});
      }
      break;
    }
    default: {
      if (qbs.empty()) {
        // ops with no quantum dependencies
        PauliNode_ptr node = std::make_shared<ClassicalNode>(args, op);
        apply_node_at_end(node);
        return;
      }
      throw BadOpType("Cannot add gate to GPGraph", type);
    }
  }
  for (auto it = pauli_rots.begin(); it != pauli_rots.end(); ++it) {
    apply_pauli_at_end(
        it->first, it->second, qbs, conditional, cond_bits, cond_value);
  }
}

std::vector<GPVert> GPGraph::vertices_in_order() const {
  GPVIndex index = boost::get(boost::vertex_index, graph_);
  int i = 0;
  BGL_FORALL_VERTICES(v, graph_, GPDAG) { boost::put(index, v, i++); }
  std::vector<GPVert> vertices;
  boost::topological_sort(graph_, std::back_inserter(vertices));
  std::reverse(vertices.begin(), vertices.end());
  return vertices;
}

std::tuple<
    std::vector<std::vector<PauliNode_ptr>>, std::vector<PauliNode_ptr>,
    boost::bimap<unsigned, unsigned>>
GPGraph::get_sequence() {
  std::vector<GPVert> vertices = vertices_in_order();
  auto it = vertices.begin();
  std::vector<std::vector<PauliNode_ptr>> interior_nodes;
  while (it != vertices.end()) {
    const PauliNode_ptr& node = graph_[*it];
    std::vector<PauliNode_ptr> commuting_set;
    commuting_set.push_back(node);
    ++it;
    while (it != vertices.end()) {
      const PauliNode_ptr& u = graph_[*it];
      bool commutes_with_all = true;
      for (const PauliNode_ptr& v : commuting_set) {
        if (!nodes_commute(u, v)) {
          commutes_with_all = false;
          break;
        }
      }
      if (!commutes_with_all) break;
      commuting_set.push_back(u);
      ++it;
    }
    interior_nodes.push_back(commuting_set);
  }
  // add clifford
  std::vector<PauliNode_ptr> cliff_nodes =
      get_nodes_from_tableau(cliff_, n_qubits_);
  return {interior_nodes, cliff_nodes, measures_};
}

}  // namespace GreedyPauliSimp

}  // namespace Transforms

}  // namespace tket
