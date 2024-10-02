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

// convert a Pauli exponential to a PauliNode_ptr
static PauliNode_ptr get_node_from_exp(
    const std::vector<Pauli>& paulis, const Expr& theta,
    const qubit_vector_t& args, unsigned n) {
  // pad the Paulis
  std::vector<Pauli> string(n, Pauli::I);
  for (unsigned i = 0; i < args.size(); i++) {
    string[args[i].index().at(0)] = paulis[i];
  }
  return std::make_shared<PauliRotation>(string, theta);
}

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

// detect trivial pauli exps, if true then return the global phase
static std::pair<bool, Expr> is_trivial_pauliexp(
    const std::vector<Pauli>& paulis, const Expr& theta) {
  if (static_cast<std::size_t>(std::count(
          paulis.begin(), paulis.end(), Pauli::I)) == paulis.size()) {
    // If all identity term
    return {true, -theta / 2};
  }
  if (equiv_0(theta, 2)) {
    if (equiv_0(theta, 4)) {
      return {true, 0};
    } else {
      return {true, -1};
    }
  }
  return {false, 0};
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

std::tuple<
    Circuit, std::vector<std::vector<PauliNode_ptr>>,
    std::vector<PauliNode_ptr>, Circuit, unit_map_t>
gpg_from_circuit(const Circuit& circ) {
  // circuit for conversion
  Circuit circ_flat(circ);
  unsigned n_qubits = circ_flat.n_qubits();
  unsigned n_bits = circ_flat.n_bits();
  // empty circuit
  Circuit empty_circ(n_qubits, n_bits);
  std::optional<std::string> name = circ_flat.get_name();
  if (name != std::nullopt) {
    empty_circ.set_name(name.value());
  }
  empty_circ.add_phase(circ_flat.get_phase());
  // measurement circuit
  Circuit measure_circ(n_qubits, n_bits);

  // flatten registers before process
  unit_map_t unit_map = circ_flat.flatten_registers();
  unit_map_t rev_unit_map;
  for (const auto& pair : unit_map) {
    rev_unit_map.insert({pair.second, pair.first});
  }

  std::vector<std::vector<PauliNode_ptr>> rotation_sets;
  std::vector<Command> commands = circ_flat.get_commands();
  Circuit cliff(n_qubits);
  // extract the final clifford and the measurement circuits
  for (const Command& cmd : commands) {
    OpType optype = cmd.get_op_ptr()->get_type();
    switch (optype) {
      case OpType::Measure: {
        measure_circ.add_op<UnitID>(OpType::Measure, cmd.get_args());
        break;
      }
      default: {
        if (optype == OpType::PauliExpBox ||
            optype == OpType::PauliExpPairBox ||
            optype == OpType::PauliExpCommutingSetBox)
          break;
        TKET_ASSERT(is_clifford_type(optype) && is_gate_type(optype));
        cliff.add_op<UnitID>(optype, cmd.get_args());
      }
    }
  }
  UnitaryRevTableau tab = circuit_to_unitary_rev_tableau(cliff);
  // convert the tableau into a set of nodes
  std::vector<PauliNode_ptr> rows = get_nodes_from_tableau(tab, n_qubits);
  // extract the Pauli exps
  for (const Command& cmd : commands) {
    OpType optype = cmd.get_op_ptr()->get_type();
    switch (optype) {
      case OpType::PauliExpBox: {
        const PauliExpBox& pbox =
            static_cast<const PauliExpBox&>(*cmd.get_op_ptr());
        const Expr phase = pbox.get_phase();
        const std::vector<Pauli> paulis = pbox.get_paulis();
        auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
        if (trivial) {
          empty_circ.add_phase(global_phase);
        } else {
          rotation_sets.push_back(
              {get_node_from_exp(paulis, phase, cmd.get_qubits(), n_qubits)});
        }
        break;
      }
      case OpType::PauliExpPairBox: {
        const PauliExpPairBox& pbox =
            static_cast<const PauliExpPairBox&>(*cmd.get_op_ptr());
        const auto [paulis1, paulis2] = pbox.get_paulis_pair();
        const auto [phase1, phase2] = pbox.get_phase_pair();
        auto [trivial1, global_phase1] = is_trivial_pauliexp(paulis1, phase1);
        auto [trivial2, global_phase2] = is_trivial_pauliexp(paulis2, phase2);
        std::vector<PauliNode_ptr> rotation_set;
        if (trivial1) {
          empty_circ.add_phase(global_phase1);
        } else {
          rotation_set.push_back(
              get_node_from_exp(paulis1, phase1, cmd.get_qubits(), n_qubits));
        }
        if (trivial2) {
          empty_circ.add_phase(global_phase2);
        } else {
          rotation_set.push_back(
              get_node_from_exp(paulis2, phase2, cmd.get_qubits(), n_qubits));
        }
        if (!rotation_set.empty()) {
          rotation_sets.push_back(rotation_set);
        }
        break;
      }
      case OpType::PauliExpCommutingSetBox: {
        const PauliExpCommutingSetBox& pbox =
            static_cast<const PauliExpCommutingSetBox&>(*cmd.get_op_ptr());
        const std::vector<SymPauliTensor> gadgets = pbox.get_pauli_gadgets();
        std::vector<PauliNode_ptr> rotation_set;
        for (const SymPauliTensor& pt : gadgets) {
          const std::vector<Pauli> paulis = pt.string;
          const Expr phase = pt.coeff;
          auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
          if (trivial) {
            empty_circ.add_phase(global_phase);
          } else {
            rotation_set.push_back(
                get_node_from_exp(paulis, phase, cmd.get_qubits(), n_qubits));
          }
        }
        if (rotation_set.size() > 0) {
          rotation_sets.push_back(rotation_set);
        }
        break;
      }
      default:
        break;
    }
  }
  return {empty_circ, rotation_sets, rows, measure_circ, rev_unit_map};
}

// given a stabiliser Pauli string and an angle, return a dense string and an
// angle
static std::pair<std::vector<Pauli>, Expr> dense_pauli(
    const SpPauliStabiliser& pauli, const unsigned& n_qubits,
    const Expr& angle) {
  bool sign = cast_coeff<quarter_turns_t, Complex>(pauli.coeff) == 1.;
  std::vector<Pauli> string(n_qubits, Pauli::I);
  for (unsigned i = 0; i < n_qubits; i++) {
    string[i] = pauli.string.at(Qubit(i));
  }
  Expr new_angle = angle;
  if (!sign) {
    new_angle *= -1;
  }
  return {string, new_angle};
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
  if (n1->get_type() == PauliNodeType::Rotation &&
      n2->get_type() == PauliNodeType::Rotation) {
    const PauliRotation& rot1 = dynamic_cast<const PauliRotation&>(*n1);
    const PauliRotation& rot2 = dynamic_cast<const PauliRotation&>(*n2);
    return strings_commute(rot1.string(), rot2.string());
  } else {
    TKET_ASSERT(false);
  }
}

std::optional<PauliNode_ptr> merge_nodes(
    const PauliNode_ptr& n1, const PauliNode_ptr& n2) {
  if (n1->get_type() == PauliNodeType::Rotation &&
      n2->get_type() == PauliNodeType::Rotation) {
    const PauliRotation& rot1 = dynamic_cast<const PauliRotation&>(*n1);
    const PauliRotation& rot2 = dynamic_cast<const PauliRotation&>(*n2);
    if (rot1.string() == rot2.string()) {
      Expr angle1 = rot1.sign() ? rot1.theta() : -rot1.theta();
      Expr angle2 = rot2.sign() ? rot2.theta() : -rot2.theta();
      return std::make_shared<PauliRotation>(rot1.string(), angle1 + angle2);
    }
  }
  return std::nullopt;
}

GPGraph::GPGraph(const unsigned& n_qubits, const unsigned& n_bits)
    : n_qubits_(n_qubits), n_bits_(n_bits), cliff_(n_qubits) {}

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
      if (node->get_type() == PauliNodeType::Rotation &&
          compare_node->get_type() == PauliNodeType::Rotation) {
        const PauliRotation& rot1 = dynamic_cast<const PauliRotation&>(*node);
        const PauliRotation& rot2 =
            dynamic_cast<const PauliRotation&>(*compare_node);
        if (rot1.string() == rot2.string()) {
          boost::clear_vertex(new_vert, graph_);
          boost::remove_vertex(new_vert, graph_);
          Expr merged_angle = (rot1.sign() ? rot1.theta() : -rot1.theta()) +
                              (rot2.sign() ? rot2.theta() : -rot2.theta());
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

void GPGraph::apply_gate_at_end(const Command& cmd) {
  const Op_ptr op = cmd.get_op_ptr();
  unit_vector_t args = cmd.get_args();
  qubit_vector_t qbs = {args.begin(), args.end()};
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

  switch (type) {
    case OpType::Measure: {
      measures_.insert({args.at(0).index().at(0), args.at(1).index().at(0)});
      break;
    }
    case OpType::Z:
    case OpType::X:
    case OpType::Y:
    case OpType::S:
    case OpType::Sdg:
    case OpType::V:
    case OpType::Vdg:
    case OpType::H:
    case OpType::CX:
    case OpType::CY:
    case OpType::CZ:
    case OpType::SWAP:
    case OpType::noop:
    case OpType::Phase: {
      cliff_.apply_gate_at_end(type, qbs);
      break;
    }
    case OpType::Rz: {
      SpPauliStabiliser pauli = cliff_.get_zrow(qbs.at(0));
      Expr angle = op->get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        for (unsigned i = 0; i < cliff_angle.value(); i++) {
          cliff_.apply_gate_at_end(OpType::S, qbs);
        }
      } else {
        auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, angle);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::Rx: {
      SpPauliStabiliser pauli = cliff_.get_xrow(qbs.at(0));
      Expr angle = op->get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        for (unsigned i = 0; i < cliff_angle.value(); i++) {
          cliff_.apply_gate_at_end(OpType::V, qbs);
        }
      } else {
        auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, angle);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::Ry: {
      Expr angle = op->get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          cliff_.apply_gate_at_end(OpType::V, qbs);
          for (unsigned i = 0; i < cliff_angle.value(); i++) {
            cliff_.apply_gate_at_end(OpType::S, qbs);
          }
          cliff_.apply_gate_at_end(OpType::Vdg, qbs);
        }
      } else {
        SpPauliStabiliser ypauli =
            cliff_.get_row_product(SpPauliStabiliser(qbs.at(0), Pauli::Y));
        auto [pauli_dense, theta] = dense_pauli(ypauli, n_qubits_, angle);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::PhasedX: {
      Expr alpha = op->get_params().at(0);
      Expr beta = op->get_params().at(1);
      std::optional<unsigned> cliff_alpha = equiv_Clifford(alpha);
      std::optional<unsigned> cliff_beta = equiv_Clifford(beta);
      // Rz(-b)
      if (cliff_beta) {
        for (unsigned i = 0; i < cliff_beta.value(); i++) {
          cliff_.apply_gate_at_end(OpType::Sdg, qbs);
        }
      } else {
        SpPauliStabiliser zpauli = cliff_.get_zrow(qbs.at(0));
        auto [pauli_dense, theta] = dense_pauli(zpauli, n_qubits_, -beta);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      // Rx(a)
      if (cliff_alpha) {
        for (unsigned i = 0; i < cliff_alpha.value(); i++) {
          cliff_.apply_gate_at_end(OpType::V, qbs);
        }
      } else {
        SpPauliStabiliser xpauli = cliff_.get_xrow(qbs.at(0));
        auto [pauli_dense, theta] = dense_pauli(xpauli, n_qubits_, alpha);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      // Rz(b)
      if (cliff_beta) {
        for (unsigned i = 0; i < cliff_beta.value(); i++) {
          cliff_.apply_gate_at_end(OpType::S, qbs);
        }
      } else {
        SpPauliStabiliser zpauli = cliff_.get_zrow(qbs.at(0));
        auto [pauli_dense, theta] = dense_pauli(zpauli, n_qubits_, beta);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::T: {
      SpPauliStabiliser pauli = cliff_.get_zrow(qbs.at(0));
      auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, 0.25);
      PauliNode_ptr node = std::make_shared<PauliRotation>(pauli_dense, theta);
      apply_node_at_end(node);
      break;
    }
    case OpType::Tdg: {
      SpPauliStabiliser pauli = cliff_.get_zrow(qbs.at(0));
      auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, -0.25);
      PauliNode_ptr node = std::make_shared<PauliRotation>(pauli_dense, theta);
      apply_node_at_end(node);
      break;
    }
    case OpType::ZZMax: {
      cliff_.apply_gate_at_end(OpType::S, {qbs.at(0)});
      cliff_.apply_gate_at_end(OpType::Z, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::S, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::V, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::S, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::CX, qbs);
      cliff_.apply_gate_at_end(OpType::S, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::V, {qbs.at(1)});
      break;
    }
    case OpType::PhaseGadget:
    case OpType::ZZPhase: {
      Expr angle = op->get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          QubitPauliMap qpm;
          for (const Qubit& q : qbs) qpm.insert({q, Pauli::Z});
          cliff_.apply_pauli_at_end(SpPauliStabiliser(qpm), *cliff_angle);
        }
      } else {
        QubitPauliMap qpm;
        for (const Qubit& q : qbs) qpm.insert({q, Pauli::Z});
        SpPauliStabiliser pauli =
            cliff_.get_row_product(SpPauliStabiliser(qpm));
        auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, angle);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::XXPhase: {
      Expr angle = op->get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          cliff_.apply_pauli_at_end(
              SpPauliStabiliser({{qbs.at(0), Pauli::X}, {qbs.at(1), Pauli::X}}),
              *cliff_angle);
        }
      } else {
        SpPauliStabiliser pauli = cliff_.get_row_product(
            SpPauliStabiliser({{qbs.at(0), Pauli::X}, {qbs.at(1), Pauli::X}}));
        auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, angle);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::YYPhase: {
      Expr angle = op->get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          cliff_.apply_pauli_at_end(
              SpPauliStabiliser({{qbs.at(0), Pauli::Y}, {qbs.at(1), Pauli::Y}}),
              *cliff_angle);
        }
      } else {
        SpPauliStabiliser pauli = cliff_.get_row_product(
            SpPauliStabiliser({{qbs.at(0), Pauli::Y}, {qbs.at(1), Pauli::Y}}));
        auto [pauli_dense, theta] = dense_pauli(pauli, n_qubits_, angle);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
      break;
    }
    case OpType::PauliExpBox: {
      const PauliExpBox& peb = static_cast<const PauliExpBox&>(*op);
      std::vector<Pauli> paulis = peb.get_paulis();
      Expr phase = peb.get_phase();
      QubitPauliMap qpm;
      for (unsigned i = 0; i != args.size(); ++i)
        qpm.insert({Qubit(args[i]), paulis[i]});
      SpPauliStabiliser qpt = cliff_.get_row_product(SpPauliStabiliser(qpm));
      auto [pauli_dense, theta] = dense_pauli(qpt, n_qubits_, phase);
      PauliNode_ptr node = std::make_shared<PauliRotation>(pauli_dense, theta);
      apply_node_at_end(node);
    }
    case OpType::PauliExpPairBox: {
      const PauliExpPairBox& peb = static_cast<const PauliExpPairBox&>(*op);
      auto [paulis1, paulis2] = peb.get_paulis_pair();
      auto [phase1, phase2] = peb.get_phase_pair();
      QubitPauliMap qpm1, qpm2;
      for (unsigned i = 0; i != args.size(); ++i) {
        qpm1.insert({Qubit(args[i]), paulis1[i]});
        qpm2.insert({Qubit(args[i]), paulis2[i]});
      }
      SpPauliStabiliser qpt1 = cliff_.get_row_product(SpPauliStabiliser(qpm1));
      SpPauliStabiliser qpt2 = cliff_.get_row_product(SpPauliStabiliser(qpm2));
      auto [pauli_dense1, theta1] = dense_pauli(qpt1, n_qubits_, phase1);
      PauliNode_ptr node1 =
          std::make_shared<PauliRotation>(pauli_dense1, theta1);
      apply_node_at_end(node1);
      auto [pauli_dense2, theta2] = dense_pauli(qpt2, n_qubits_, phase2);
      PauliNode_ptr node2 =
          std::make_shared<PauliRotation>(pauli_dense2, theta2);
      apply_node_at_end(node2);
    }
    case OpType::PauliExpCommutingSetBox: {
      const PauliExpCommutingSetBox& peb =
          static_cast<const PauliExpCommutingSetBox&>(*op);
      for (const SymPauliTensor& pt : peb.get_pauli_gadgets()) {
        const std::vector<Pauli> paulis = pt.string;
        const Expr phase = pt.coeff;
        QubitPauliMap qpm;
        for (unsigned i = 0; i != args.size(); ++i) {
          qpm.insert({Qubit(args[i]), paulis[i]});
        }
        SpPauliStabiliser qpt = cliff_.get_row_product(SpPauliStabiliser(qpm));
        auto [pauli_dense, theta] = dense_pauli(qpt, n_qubits_, phase);
        PauliNode_ptr node =
            std::make_shared<PauliRotation>(pauli_dense, theta);
        apply_node_at_end(node);
      }
    }
    default: {
      throw BadOpType("Cannot add gate to GPGraph", type);
    }
  }
}

}  // namespace GreedyPauliSimp

}  // namespace Transforms

}  // namespace tket
