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

#include <tkassert/Assert.hpp>

#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"
#include "tket/Gate/Gate.hpp"
#include "tket/Ops/OpJsonFactory.hpp"
#include "tket/PauliGraph/ConjugatePauliFunctions.hpp"
#include "tket/PauliGraphRefactor/Converters.hpp"
#include "tket/Transformations/PauliOptimisation.hpp"

namespace tket {

using namespace pg;

// TODO:: CHECK GLOBAL PHASES FOR EACH OF THESE OPS WITH TESTS
std::pair<std::vector<PGOp_ptr>, Expr> op_to_pgops(
    const Op_ptr& op, const unit_vector_t& args, UnitaryRevTableau& tab,
    ChoiAPState& ap, bool allow_tableau) {
  if (allow_tableau && is_clifford_type(op->get_type())) {
    qubit_vector_t qs;
    for (const UnitID& a : args) qs.push_back(Qubit(a));
    tab.apply_gate_at_end(op->get_type(), qs);
    ap.apply_gate(op->get_type(), qs, ChoiAPState::TableauSegment::Output);
    return {{}, 0};
  }
  switch (op->get_type()) {
    case OpType::Z: {
      Qubit q(args.front());
      return {{std::make_shared<PGCliffordRot>(tab.get_zrow(q), 2)}, 0};
    }
    case OpType::X: {
      Qubit q(args.front());
      return {{std::make_shared<PGCliffordRot>(tab.get_xrow(q), 2)}, 0};
    }
    case OpType::Y: {
      Qubit q(args.front());
      return {
          {std::make_shared<PGCliffordRot>(
              tab.get_row_product(SpPauliStabiliser(q, Pauli::Y)), 2)},
          0};
    }
    case OpType::S: {
      Qubit q(args.front());
      return {{std::make_shared<PGCliffordRot>(tab.get_zrow(q), 1)}, 0};
    }
    case OpType::Sdg: {
      Qubit q(args.front());
      return {{std::make_shared<PGCliffordRot>(tab.get_zrow(q), 3)}, 0};
    }
    case OpType::V: {
      Qubit q(args.front());
      return {{std::make_shared<PGCliffordRot>(tab.get_xrow(q), 1)}, -0.25};
    }
    case OpType::Vdg: {
      Qubit q(args.front());
      return {{std::make_shared<PGCliffordRot>(tab.get_xrow(q), 3)}, 0.25};
    }
    case OpType::H: {
      Qubit q(args.front());
      SpPauliStabiliser zq = tab.get_zrow(q);
      return {
          {std::make_shared<PGCliffordRot>(zq, 1),
           std::make_shared<PGCliffordRot>(tab.get_xrow(q), 1),
           std::make_shared<PGCliffordRot>(zq, 1)},
          0};
    }
    case OpType::CX: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser xt = tab.get_xrow(t);
      return {
          {std::make_shared<PGCliffordRot>(zc, 3),
           std::make_shared<PGCliffordRot>(xt, 3),
           std::make_shared<PGCliffordRot>(zc * xt, 1)},
          0};
    }
    case OpType::CY: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser yt =
          tab.get_row_product(SpPauliStabiliser(t, Pauli::Y));
      return {
          {std::make_shared<PGCliffordRot>(zc, 3),
           std::make_shared<PGCliffordRot>(yt, 3),
           std::make_shared<PGCliffordRot>(zc * yt, 1)},
          0};
    }
    case OpType::CZ: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser zt = tab.get_zrow(t);
      return {
          {std::make_shared<PGCliffordRot>(zc, 3),
           std::make_shared<PGCliffordRot>(zt, 3),
           std::make_shared<PGCliffordRot>(zc * zt, 1)},
          0};
    }
    case OpType::ZZMax: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser zt = tab.get_zrow(t);
      return {{std::make_shared<PGCliffordRot>(zc * zt, 1)}, 0};
    }
    case OpType::SWAP: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser zt = tab.get_zrow(t);
      SpPauliStabiliser xc = tab.get_xrow(c);
      SpPauliStabiliser xt = tab.get_xrow(t);
      SpPauliStabiliser yy = zc * xc * zt * xt;
      yy.coeff = (yy.coeff + 2) % 4;
      return {
          {std::make_shared<PGCliffordRot>(zc * zt, 1),
           std::make_shared<PGCliffordRot>(xc * xt, 1),
           std::make_shared<PGCliffordRot>(yy, 1)},
          0};
    }
    case OpType::T: {
      Qubit q(args.front());
      return {{std::make_shared<PGRotation>(tab.get_zrow(q), 0.25)}, 0};
    }
    case OpType::Tdg: {
      Qubit q(args.front());
      return {{std::make_shared<PGRotation>(tab.get_zrow(q), -0.25)}, 0};
    }
    case OpType::Rz: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {
          {std::make_shared<PGRotation>(tab.get_zrow(q), g.get_params().at(0))},
          g.get_params().at(0)};
    }
    case OpType::Rx: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {
          {std::make_shared<PGRotation>(tab.get_xrow(q), g.get_params().at(0))},
          0};
    }
    case OpType::Ry: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {
          {std::make_shared<PGRotation>(
              tab.get_row_product(SpPauliStabiliser(q, Pauli::Y)),
              g.get_params().at(0))},
          0};
    }
    case OpType::TK1: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser zq = tab.get_zrow(q);
      return {
          {std::make_shared<PGRotation>(zq, g.get_params().at(2)),
           std::make_shared<PGRotation>(tab.get_xrow(q), g.get_params().at(1)),
           std::make_shared<PGRotation>(zq, g.get_params().at(0))},
          0};
    }
    case OpType::PhaseGadget: {
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliMap qpm;
      for (const UnitID& a : args) qpm.insert({Qubit(a), Pauli::Z});
      return {
          {std::make_shared<PGRotation>(
              tab.get_row_product(SpPauliStabiliser(qpm)),
              g.get_params().at(0))},
          0};
    }
    case OpType::ZZPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser z0 = tab.get_zrow(q0);
      SpPauliStabiliser z1 = tab.get_zrow(q1);
      return {{std::make_shared<PGRotation>(z0 * z1, g.get_params().at(0))}, 0};
    }
    case OpType::XXPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser x0 = tab.get_xrow(q0);
      SpPauliStabiliser x1 = tab.get_xrow(q1);
      return {{std::make_shared<PGRotation>(x0 * x1, g.get_params().at(0))}, 0};
    }
    case OpType::YYPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser yy = tab.get_row_product(
          SpPauliStabiliser({{q0, Pauli::Y}, {q1, Pauli::Y}}));
      return {{std::make_shared<PGRotation>(yy, g.get_params().at(0))}, 0};
    }
    case OpType::TK2: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser z0 = tab.get_zrow(q0);
      SpPauliStabiliser z1 = tab.get_zrow(q1);
      SpPauliStabiliser x0 = tab.get_xrow(q0);
      SpPauliStabiliser x1 = tab.get_xrow(q1);
      SpPauliStabiliser yy = z0 * x0 * z1 * x1;
      yy.coeff = (yy.coeff + 2) % 4;
      return {
          {std::make_shared<PGRotation>(x0 * x1, g.get_params().at(0)),
           std::make_shared<PGRotation>(yy, g.get_params().at(1)),
           std::make_shared<PGRotation>(z0 * z1, g.get_params().at(2))},
          0};
    }
    case OpType::Measure: {
      return {
          {std::make_shared<PGMeasure>(
              tab.get_zrow(Qubit(args.at(0))), Bit(args.at(1)))},
          0};
    }
    case OpType::Collapse: {
      return {
          {std::make_shared<PGDecoherence>(tab.get_zrow(Qubit(args.front())))},
          0};
    }
    case OpType::Reset: {
      Qubit q(args.front());
      return {{std::make_shared<PGReset>(tab.get_zrow(q), tab.get_xrow(q))}, 0};
    }
    case OpType::PauliExpBox: {
      const PauliExpBox& box = dynamic_cast<const PauliExpBox&>(*op);
      const std::vector<Pauli>& paulis = box.get_paulis();
      QubitPauliMap qpm;
      for (unsigned i = 0; i < args.size(); ++i)
        qpm.insert({Qubit(args.at(i)), paulis.at(i)});
      return {
          {std::make_shared<PGRotation>(
              tab.get_row_product(SpPauliStabiliser(qpm)), box.get_phase())},
          0};
    }
    case OpType::QControlBox: {
      const QControlBox& box = dynamic_cast<const QControlBox&>(*op);
      unsigned n_controls = box.get_n_controls();
      std::vector<bool> control_state = box.get_control_state();
      unit_vector_t inner_args;
      for (unsigned i = n_controls; i < args.size(); ++i)
        inner_args.push_back(args.at(i));
      std::pair<std::vector<PGOp_ptr>, Expr> inner_ops =
          op_to_pgops(box.get_op(), inner_args, tab, ap, false);
      std::vector<SpPauliStabiliser> control_paulis;
      for (unsigned i = 0; i < n_controls; ++i)
        control_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
      std::vector<PGOp_ptr> ret_ops;
      for (const PGOp_ptr& inn_op : inner_ops.first)
        ret_ops.push_back(std::make_shared<PGQControl>(
            inn_op, control_paulis, control_state));
      if (inner_ops.second == 0 || n_controls == 0)
        return {ret_ops, inner_ops.second};
      else {
        SpPauliStabiliser last_control = control_paulis.back();
        control_paulis.pop_back();
        bool last_control_state = control_state.back();
        control_state.pop_back();
        if (last_control_state) {
          // apply phase on -1 eigenstate, same as a rotation
          ret_ops.push_back(std::make_shared<PGQControl>(
              std::make_shared<PGRotation>(last_control, inner_ops.second),
              control_paulis, control_state));
          return {ret_ops, 0};
        } else {
          // apply phase on +1 eigenstate, so instead apply phase globally and
          // apply negative phase on the -1 eigenstate
          ret_ops.push_back(std::make_shared<PGQControl>(
              std::make_shared<PGRotation>(last_control, -inner_ops.second),
              control_paulis, control_state));
          return {ret_ops, inner_ops.second};
        }
      }
    }
    case OpType::MultiplexorBox: {
      const MultiplexorBox& box = dynamic_cast<const MultiplexorBox&>(*op);
      std::map<std::vector<bool>, std::vector<Op_ptr>> op_map;
      for (const std::pair<const std::vector<bool>, Op_ptr>& qcase :
           box.get_ops())
        op_map.insert({qcase.first, {qcase.second}});
      unsigned n_controls = box.get_op_map().begin()->first.size();
      std::vector<SpPauliStabiliser> control_paulis;
      for (unsigned i = 0; i < n_controls; ++i)
        control_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
      std::vector<SpPauliStabiliser> target_paulis;
      for (unsigned i = n_controls; i < args.size(); ++i) {
        target_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
        target_paulis.push_back(tab.get_xrow(Qubit(args.at(i))));
      }
      return {
          {std::make_shared<PGMultiplexedTensoredBox>(
              op_map, control_paulis, target_paulis)},
          0};
    }
    case OpType::MultiplexedRotationBox: {
      const MultiplexedRotationBox& box =
          dynamic_cast<const MultiplexedRotationBox&>(*op);
      std::map<std::vector<bool>, Expr> angle_map;
      for (const std::pair<const std::vector<bool>, Op_ptr>& qcase :
           box.get_ops()) {
        const Gate& g = dynamic_cast<const Gate&>(*qcase.second);
        angle_map.insert({qcase.first, g.get_params().at(0)});
      }
      unsigned n_controls = box.get_op_map().begin()->first.size();
      std::vector<SpPauliStabiliser> control_paulis;
      for (unsigned i = 0; i < n_controls; ++i)
        control_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
      OpType rotation_basis = box.get_ops().begin()->second->get_type();
      switch (rotation_basis) {
        case OpType::Rx: {
          return {
              {std::make_shared<PGMultiplexedRotation>(
                  angle_map, control_paulis,
                  tab.get_xrow(Qubit(args.at(n_controls))))},
              0};
        }
        case OpType::Ry: {
          return {
              {std::make_shared<PGMultiplexedRotation>(
                  angle_map, control_paulis,
                  tab.get_row_product(SpPauliStabiliser(
                      Qubit(args.at(n_controls)), Pauli::Y)))},
              0};
        }
        case OpType::Rz: {
          return {
              {std::make_shared<PGMultiplexedRotation>(
                  angle_map, control_paulis,
                  tab.get_zrow(Qubit(args.at(n_controls))))},
              0};
        }
        default: {
          throw PGError(
              "Invalid MultiplexedRotationBox encountered when converting to "
              "PauliGraph");
        }
      }
    }
    case OpType::MultiplexedU2Box: {
      const MultiplexedU2Box& box = dynamic_cast<const MultiplexedU2Box&>(*op);
      std::map<std::vector<bool>, std::vector<Op_ptr>> op_map;
      for (const std::pair<const std::vector<bool>, Op_ptr>& qcase :
           box.get_ops())
        op_map.insert({qcase.first, {qcase.second}});
      unsigned n_controls = box.get_op_map().begin()->first.size();
      std::vector<SpPauliStabiliser> control_paulis;
      for (unsigned i = 0; i < n_controls; ++i)
        control_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
      std::vector<SpPauliStabiliser> target_paulis;
      target_paulis.push_back(tab.get_zrow(Qubit(args.at(n_controls))));
      target_paulis.push_back(tab.get_xrow(Qubit(args.at(n_controls))));
      return {
          {std::make_shared<PGMultiplexedTensoredBox>(
              op_map, control_paulis, target_paulis)},
          0};
    }
    case OpType::MultiplexedTensoredU2Box: {
      const MultiplexedTensoredU2Box& box =
          dynamic_cast<const MultiplexedTensoredU2Box&>(*op);
      unsigned n_controls = box.get_op_map().begin()->first.size();
      std::vector<SpPauliStabiliser> control_paulis;
      for (unsigned i = 0; i < n_controls; ++i)
        control_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
      std::vector<SpPauliStabiliser> target_paulis;
      for (unsigned i = n_controls; i < args.size(); ++i) {
        target_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
        target_paulis.push_back(tab.get_xrow(Qubit(args.at(i))));
      }
      return {
          {std::make_shared<PGMultiplexedTensoredBox>(
              box.get_op_map(), control_paulis, target_paulis)},
          0};
    }
    case OpType::StabiliserAssertionBox: {
      const StabiliserAssertionBox& box =
          dynamic_cast<const StabiliserAssertionBox&>(*op);
      PauliStabiliserVec stabs = box.get_stabilisers();
      unsigned n_qbs = args.size() - stabs.size();
      Qubit anc(args.at(n_qbs - 1));
      SpPauliStabiliser anc_z = tab.get_zrow(anc);
      SpPauliStabiliser anc_x = tab.get_xrow(anc);
      std::vector<PGOp_ptr> ops;
      for (unsigned i = 0; i < stabs.size(); ++i) {
        const PauliStabiliser& stab = stabs.at(i);
        Bit target(args.at(n_qbs + i));
        QubitPauliMap qpm;
        for (unsigned q = 0; q < stab.string.size(); ++q) {
          qpm.insert({Qubit(args.at(q)), stab.string.at(q)});
        }
        SpPauliStabiliser qpt(qpm, stab.coeff);
        SpPauliStabiliser prod = tab.get_row_product(qpt);
        ops.push_back(
            std::make_shared<PGStabAssertion>(prod, anc_z, anc_x, target));
      }
      return {ops, 0};
    }
    case OpType::Conditional: {
      const Conditional& cond = dynamic_cast<const Conditional&>(*op);
      bit_vector_t cond_bits;
      unit_vector_t inner_args;
      for (unsigned i = 0; i < cond.get_width(); ++i)
        cond_bits.push_back(Bit(args.at(i)));
      for (unsigned i = cond.get_width(); i < args.size(); ++i)
        inner_args.push_back(args.at(i));
      // Ignore global phase between conditional branches
      std::vector<PGOp_ptr> inner_ops =
          op_to_pgops(cond.get_op(), inner_args, tab, ap, false).first;
      std::vector<PGOp_ptr> ret;
      for (const PGOp_ptr& inn_op : inner_ops)
        ret.push_back(std::make_shared<PGConditional>(
            inn_op, cond_bits, cond.get_value()));
      return {ret, 0};
    }
    default: {
      std::vector<SpPauliStabiliser> paulis;
      for (const UnitID& uid : args) {
        if (uid.type() == UnitType::Qubit) {
          Qubit q(uid);
          paulis.push_back(tab.get_zrow(q));
          paulis.push_back(tab.get_xrow(q));
        }
      }
      return {{std::make_shared<PGBox>(op, args, paulis)}, 0};
    }
  }
}

pg::PauliGraph circuit_to_pauli_graph3(
    const Circuit& circ, bool collect_cliffords) {
  pg::PauliGraph res;
  ChoiMixTableau initial(circ.all_qubits());
  for (const Qubit& q : circ.created_qubits())
    initial.post_select(q, ChoiMixTableau::TableauSegment::Input);
  ChoiAPState initial_ap = cm_tableau_to_choi_apstate(initial);
  initial_ap.ap_.phase = circ.get_phase();
  res.add_vertex_at_end(std::make_shared<PGInputTableau>(initial, initial_ap));
  UnitaryRevTableau final_u(circ.all_qubits());
  ChoiAPState final_ap(circ.all_qubits());
  for (const Command& com : circ) {
    unit_vector_t args = com.get_args();
    std::pair<std::vector<PGOp_ptr>, Expr> pgops = op_to_pgops(
        com.get_op_ptr(), args, final_u, final_ap, collect_cliffords);
    for (const PGOp_ptr& pgop : pgops.first) res.add_vertex_at_end(pgop);
  }
  std::list<ChoiMixTableau::row_tensor_t> final_rows;
  qubit_map_t implicit_perm = circ.implicit_qubit_permutation();
  for (const Qubit& q : final_u.get_qubits()) {
    Qubit out_q = implicit_perm.at(q);
    final_rows.push_back(
        {final_u.get_zrow(q), SpPauliStabiliser(out_q, Pauli::Z)});
    final_rows.push_back(
        {final_u.get_xrow(q), SpPauliStabiliser(out_q, Pauli::X)});
  }
  ChoiMixTableau final_cm(final_rows);
  for (const Qubit& q : circ.discarded_qubits()) final_cm.discard_qubit(q);
  for (const Qubit& q : circ.discarded_qubits()) final_ap.discard_qubit(q);
  final_ap.normal_form();
  res.add_vertex_at_end(std::make_shared<PGOutputTableau>(final_cm, final_ap));
  return res;
}

/**
 * A helper struct to represent the circuit used to synthesise a single PGOp.
 *
 * The general structure used for individual synthesis is C^\dagger O C, where O
 * is some op (in this case, op applied to args), and C is a Clifford circuit.
 *
 * This Clifford circuit is not unique, but must be a solution to a particular
 * diagonalisation problem, reducing the active Pauli strings of the PGOp to Z
 * or X on particular qubits, just so that the op is conjugated correctly.
 *
 * We represent this diagonalisation problem with a tableau mapping active Pauli
 * strings on the input space to args on the output space.
 *
 * In the case of recursive PGOps such as conditionals and qcontrols, the
 * tableau must, in particular, solve the constraints of the inner PGOp, so it
 * is enough to just add new rows to solve the controls in addition.
 */
struct IndividualSynth {
  ChoiMixTableau tableau;
  Op_ptr op;
  unit_vector_t args;
};

enum class MultiplexorCase { Rotation, U2, TensoredU2, General };

IndividualSynth synth_pgop(const PGOp_ptr& pgop) {
  switch (pgop->get_type()) {
    case PGOpType::Rotation: {
      PGRotation& rot = dynamic_cast<PGRotation&>(*pgop);
      return IndividualSynth{
          ChoiMixTableau(
              {{rot.port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}}),
          get_op_ptr(OpType::U1, rot.get_angle()),
          {Qubit(0)}};
    }
    case PGOpType::CliffordRot: {
      PGCliffordRot& crot = dynamic_cast<PGCliffordRot&>(*pgop);
      Op_ptr op;
      switch (crot.get_angle() % 4) {
        case 1: {
          op = get_op_ptr(OpType::S);
          break;
        }
        case 2: {
          op = get_op_ptr(OpType::Z);
          break;
        }
        case 3: {
          op = get_op_ptr(OpType::Sdg);
          break;
        }
        default: {
          op = get_op_ptr(OpType::noop);
        }
      }
      return IndividualSynth{
          ChoiMixTableau(
              {{crot.port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}}),
          op,
          {Qubit(0)}};
    }
    case PGOpType::Measure: {
      PGMeasure& meas = dynamic_cast<PGMeasure&>(*pgop);
      return IndividualSynth{
          ChoiMixTableau(
              {{meas.port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}}),
          get_op_ptr(OpType::Measure),
          {Qubit(0), meas.get_target()}};
    }
    case PGOpType::Decoherence: {
      return IndividualSynth{
          ChoiMixTableau(
              {{pgop->port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}}),
          get_op_ptr(OpType::Collapse),
          {Qubit(0)}};
    }
    case PGOpType::Reset: {
      PGReset& reset = dynamic_cast<PGReset&>(*pgop);
      return IndividualSynth{
          ChoiMixTableau(
              {{reset.get_stab(), SpPauliStabiliser(Qubit(0), Pauli::Z)},
               {reset.get_destab(), SpPauliStabiliser(Qubit(0), Pauli::X)}}),
          get_op_ptr(OpType::Reset),
          {Qubit(0)}};
    }
    case PGOpType::Conditional: {
      PGConditional& cond = dynamic_cast<PGConditional&>(*pgop);
      IndividualSynth inner_synth = synth_pgop(cond.get_inner_op());
      unit_vector_t args;
      for (const Bit& b : cond.get_args()) args.push_back(b);
      for (const UnitID& u : inner_synth.args) args.push_back(u);
      return IndividualSynth{
          inner_synth.tableau,
          std::make_shared<Conditional>(
              inner_synth.op, cond.get_args().size(), cond.get_value()),
          args};
    }
    case PGOpType::QControl: {
      PGQControl& qctrl = dynamic_cast<PGQControl&>(*pgop);
      IndividualSynth synth = synth_pgop(qctrl.get_inner_op());
      unit_vector_t ctrl_qubits;
      unsigned qi = synth.args.size();
      for (const SpPauliStabiliser& cp : qctrl.get_control_paulis()) {
        Qubit q(qi);
        synth.tableau.add_row({cp, SpPauliStabiliser(q, Pauli::Z)});
        ctrl_qubits.push_back(q);
        ++qi;
      }
      synth.args.insert(
          synth.args.begin(), ctrl_qubits.begin(), ctrl_qubits.end());
      synth.op = std::make_shared<QControlBox>(
          synth.op, qctrl.get_value().size(), qctrl.get_value());
      return synth;
    }
    case PGOpType::MultiplexedTensoredBox: {
      PGMultiplexedTensoredBox& mptb =
          dynamic_cast<PGMultiplexedTensoredBox&>(*pgop);
      IndividualSynth synth;
      const std::vector<SpPauliStabiliser>& control_paulis =
          mptb.get_control_paulis();
      for (unsigned qi = 0; qi < control_paulis.size(); ++qi) {
        Qubit q(qi);
        synth.tableau.add_row(
            {control_paulis.at(qi), SpPauliStabiliser(q, Pauli::Z)});
        synth.args.push_back(q);
      }
      const std::vector<SpPauliStabiliser>& target_paulis =
          mptb.get_target_paulis();
      for (unsigned qi = 0; 2 * qi < target_paulis.size(); ++qi) {
        Qubit q(control_paulis.size() + qi);
        synth.tableau.add_row(
            {target_paulis.at(2 * qi), SpPauliStabiliser(q, Pauli::Z)});
        synth.tableau.add_row(
            {target_paulis.at(2 * qi + 1), SpPauliStabiliser(q, Pauli::X)});
        synth.args.push_back(q);
      }
      bool all_single_qubit = true;
      for (const std::pair<const std::vector<bool>, std::vector<Op_ptr>>&
               qcase : mptb.get_op_map()) {
        for (const Op_ptr& op : qcase.second) {
          if (op->n_qubits() != 1) {
            all_single_qubit = false;
            break;
          }
        }
        if (!all_single_qubit) break;
      }
      if (all_single_qubit) {
        // Fits into either a MultiplexedU2Box or a MultiplexedTensoredU2Box
        bool all_one_op = true;
        for (const std::pair<const std::vector<bool>, std::vector<Op_ptr>>&
                 qcase : mptb.get_op_map()) {
          if (qcase.second.size() != 1) {
            all_one_op = false;
            break;
          }
        }
        if (all_one_op) {
          // Fits into a MultiplexedU2Box
          ctrl_op_map_t op_map;
          for (const std::pair<const std::vector<bool>, std::vector<Op_ptr>>&
                   qcase : mptb.get_op_map())
            op_map.insert({qcase.first, qcase.second.front()});
          synth.op = std::make_shared<MultiplexedU2Box>(op_map);
          return synth;
        } else {
          // Fits into a MultiplexedTensoredU2Box
          synth.op =
              std::make_shared<MultiplexedTensoredU2Box>(mptb.get_op_map());
          return synth;
        }
      } else {
        // Catch-all case, encode into a MultiplexorBox
        ctrl_op_map_t op_map;
        for (const std::pair<const std::vector<bool>, std::vector<Op_ptr>>&
                 qcase : mptb.get_op_map()) {
          if (qcase.second.size() == 1) {
            op_map.insert({qcase.first, qcase.second.front()});
          } else {
            // Encode a tensor of operations into a CircBox to reduce the
            // collection to a single Op
            Circuit c;
            for (const Op_ptr& op : qcase.second) {
              unsigned nq = op->n_qubits();
              qubit_vector_t args;
              unsigned c_size = c.n_qubits();
              for (unsigned qi = 0; qi < nq; ++qi) {
                c.add_qubit(Qubit(c_size + qi));
                args.push_back(Qubit(c_size + qi));
              }
              c.add_op<Qubit>(op, args);
            }
            op_map.insert({qcase.first, std::make_shared<CircBox>(c)});
          }
        }
        synth.op = std::make_shared<MultiplexorBox>(op_map);
        return synth;
      }
    }
    case PGOpType::MultiplexedRotation: {
      PGMultiplexedRotation& mpr = dynamic_cast<PGMultiplexedRotation&>(*pgop);
      IndividualSynth synth;
      const std::vector<SpPauliStabiliser>& control_paulis =
          mpr.get_control_paulis();
      for (unsigned qi = 0; qi < control_paulis.size(); ++qi) {
        Qubit q(qi);
        synth.tableau.add_row(
            {control_paulis.at(qi), SpPauliStabiliser(q, Pauli::Z)});
        synth.args.push_back(q);
      }
      synth.tableau.add_row(
          {mpr.get_target_pauli(),
           SpPauliStabiliser(Qubit(control_paulis.size()), Pauli::Z)});
      synth.args.push_back(Qubit(control_paulis.size()));
      ctrl_op_map_t op_map;
      for (const std::pair<const std::vector<bool>, Expr>& qcase :
           mpr.get_angle_map()) {
        op_map.insert({qcase.first, get_op_ptr(OpType::Rz, qcase.second)});
      }
      synth.op = std::make_shared<MultiplexedRotationBox>(op_map);
      return synth;
    }
    case PGOpType::Box: {
      PGBox& box = dynamic_cast<PGBox&>(*pgop);
      IndividualSynth res;
      for (unsigned i = 0; 2 * i < box.n_paulis(); ++i) {
        res.tableau.add_row(
            {box.port(2 * i), SpPauliStabiliser(Qubit(i), Pauli::Z)});
        res.tableau.add_row(
            {box.port(2 * i + 1), SpPauliStabiliser(Qubit(i), Pauli::X)});
      }
      res.op = box.get_op();
      unsigned qb_i = 0;
      for (const UnitID& u : box.get_args()) {
        if (u.type() == UnitType::Qubit) {
          res.args.push_back(Qubit(qb_i));
          ++qb_i;
        } else
          res.args.push_back(u);
      }
      return res;
    }
    case PGOpType::StabAssertion: {
      PGStabAssertion& stab = dynamic_cast<PGStabAssertion&>(*pgop);
      return IndividualSynth{
          ChoiMixTableau({
              {stab.get_stab(), SpPauliStabiliser(Qubit(0), Pauli::Z)},
              {stab.get_anc_z(), SpPauliStabiliser(Qubit(1), Pauli::Z)},
              {stab.get_anc_x(), SpPauliStabiliser(Qubit(1), Pauli::X)},
          }),
          std::make_shared<StabiliserAssertionBox>(
              PauliStabiliserVec{PauliStabiliser({Pauli::Z}, 0)}),
          {Qubit(0), Qubit(1), stab.get_target()}};
    }
    default: {
      throw std::logic_error(
          "Error during PauliGraph synthesis: unexpected PGOpType in "
          "synth_pgop");
    }
  }
}

Circuit pgop_to_circuit(const PGOp_ptr& pgop) {
  IndividualSynth synth_data = synth_pgop(pgop);
  auto [diag_circ, qb_map] =
      cm_tableau_to_unitary_extension_circuit(synth_data.tableau);
  Circuit diag_dag = diag_circ.dagger();
  unit_vector_t new_args;
  for (const UnitID& u : synth_data.args) {
    if (u.type() == UnitType::Qubit) {
      new_args.push_back(qb_map.at(Qubit(u)));
    } else {
      diag_circ.add_bit(Bit(u), false);
      new_args.push_back(u);
    }
  }
  diag_circ.add_op<UnitID>(synth_data.op, new_args);
  diag_circ.append(diag_dag);
  return diag_circ;
}

Circuit pauli_graph3_to_circuit_individual(
    const pg::PauliGraph& pg, CXConfigType cx_config) {
  const std::set<Qubit>& qubits = pg.get_qubits();
  const std::set<Bit>& bits = pg.get_bits();
  Circuit circ(
      qubit_vector_t{qubits.begin(), qubits.end()},
      bit_vector_t{bits.begin(), bits.end()});
  std::list<PGOp_ptr> pgop_sequence = pg.pgop_sequence();
  for (const PGOp_ptr& pgop : pgop_sequence) {
    if (pgop->get_type() == PGOpType::InputTableau) {
      PGInputTableau& tab_op = dynamic_cast<PGInputTableau&>(*pgop);
      const ChoiAPState& ap = tab_op.get_apstate();
      // Set this as the new circuit so that any initialisations are applied
      // without adding Reset operations; need to readd all bits and any qubits
      // that weren't included in the tableau
      qubit_map_t perm;
      std::tie(circ, perm) = choi_apstate_to_exact_circuit(ap, cx_config);
      qubit_map_t rev_perm;
      for (const std::pair<const Qubit, Qubit>& p : perm)
        rev_perm.insert({p.second, p.first});
      circ.permute_boundary_output(perm);
      for (const Qubit& q : qubits) circ.add_qubit(q, false);
      for (const Bit& b : bits) circ.add_bit(b, false);
    } else if (pgop->get_type() == PGOpType::OutputTableau) {
      PGOutputTableau& tab_op = dynamic_cast<PGOutputTableau&>(*pgop);
      ChoiAPState ap = tab_op.get_apstate();
      auto [tab_circ, perm] = choi_apstate_to_exact_circuit(ap, cx_config);
      qubit_map_t rev_perm;
      for (const std::pair<const Qubit, Qubit>& p : perm)
        rev_perm.insert({p.second, p.first});
      tab_circ.permute_boundary_output(rev_perm);
      circ.append(tab_circ);
    } else {
      circ.append(pgop_to_circuit(pgop));
    }
  }
  return circ;
}

Circuit pauli_graph3_to_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config) {
  const std::set<Qubit>& qubits = pg.get_qubits();
  const std::set<Bit>& bits = pg.get_bits();
  Circuit circ(
      qubit_vector_t{qubits.begin(), qubits.end()},
      bit_vector_t{bits.begin(), bits.end()});
  std::list<std::list<PGOp_ptr>> commuting_sets = pg.pgop_commuting_sets();
  // Synthesise input tableau
  std::optional<PGVert> itab_v = pg.get_input_tableau();
  if (itab_v) {
    PGOp_ptr pgop = pg.get_vertex_PGOp_ptr(*itab_v);
    PGInputTableau& tab_op = dynamic_cast<PGInputTableau&>(*pgop);
    const ChoiAPState& ap = tab_op.get_apstate();
    // Set this as the new circuit so that any initialisations are applied
    // without adding Reset operations; need to read all bits and any qubits
    // that weren't included in the tableau
    qubit_map_t perm;
    std::tie(circ, perm) = choi_apstate_to_exact_circuit(ap, cx_config);
    qubit_map_t rev_perm;
    for (const std::pair<const Qubit, Qubit>& p : perm)
      rev_perm.insert({p.second, p.first});
    circ.permute_boundary_output(rev_perm);
    for (const Qubit& q : qubits) circ.add_qubit(q, false);
    for (const Bit& b : bits) circ.add_bit(b, false);
    // Remove the input tableau from the first commuting set
    auto first_set = commuting_sets.begin();
    for (auto it = first_set->begin(); it != first_set->end(); ++it) {
      if ((*it)->get_type() == PGOpType::InputTableau) {
        first_set->erase(it);
        break;
      }
    }
    if (first_set->empty()) commuting_sets.erase(first_set);
  }
  std::optional<PGVert> otab_v = pg.get_output_tableau();
  if (otab_v) {
    // Remove the output tableau from the last commuting set
    auto last_set = commuting_sets.rbegin();
    for (auto it = last_set->begin(); it != last_set->end(); ++it) {
      if ((*it)->get_type() == PGOpType::OutputTableau) {
        last_set->erase(it);
        break;
      }
    }
    if (last_set->empty()) {
      // erase doesn't work with reverse iterators, so have to cast down to
      // regular iterators, adjusting for the off-by-one positional difference
      // between forwards and reverse iterators
      last_set = std::list<std::list<PGOp_ptr>>::reverse_iterator(
          commuting_sets.erase(std::next(last_set).base()));
    }
  }

  boost::bimap<Qubit, unsigned> qubit_indices;
  boost::bimap<Bit, unsigned> bit_indices;
  unit_vector_t args;
  for (const Qubit& q : qubits) {
    qubit_indices.insert({q, (unsigned)qubit_indices.size()});
    args.push_back(q);
  }
  for (const Bit& b : bits) {
    bit_indices.insert({b, (unsigned)bit_indices.size()});
    args.push_back(b);
  }
  // Synthesise each interior commuting set
  for (const std::list<PGOp_ptr>& cset : commuting_sets) {
    Op_ptr box = std::make_shared<PGOpCommutingSetBox>(
        std::vector<PGOp_ptr>{cset.begin(), cset.end()}, qubit_indices,
        bit_indices, cx_config);
    // Append to the circuit
    circ.add_op<UnitID>(box, args);
  }

  // Synthesise output tableau
  if (otab_v) {
    PGOp_ptr pgop = pg.get_vertex_PGOp_ptr(*otab_v);
    PGOutputTableau& tab_op = dynamic_cast<PGOutputTableau&>(*pgop);
    ChoiAPState ap = tab_op.get_apstate();
    auto [tab_circ, perm] = choi_apstate_to_exact_circuit(ap, cx_config);
    qubit_map_t rev_perm;
    for (const std::pair<const Qubit, Qubit>& p : perm)
      rev_perm.insert({p.second, p.first});
    tab_circ.permute_boundary_output(rev_perm);
    circ.append(tab_circ);
  }

  return circ;
}

/*******************************************************************************
 * LEGACY SYNTHESIS METHODS
 ******************************************************************************/

// A helper method which, given the output tableau as a UnitaryTableau,
// separates the ops into rotations and measures, commuting all measurements to
// the end and checking they can all be performed simultaneously
std::pair<std::list<PGOp_ptr>, std::map<Qubit, Bit>> rotations_and_end_measures(
    const pg::PauliGraph& pg, const std::optional<UnitaryTableau>& out_tab) {
  std::list<PGOp_ptr> all_pgops = pg.pgop_sequence_boost();
  std::list<PGOp_ptr> rotations;
  std::list<PGOp_ptr> measures;
  for (const PGOp_ptr& pgp : all_pgops) {
    switch (pgp->get_type()) {
      case PGOpType::InputTableau:
      case PGOpType::OutputTableau: {
        // Start or end of list
        break;
      }
      case PGOpType::Rotation:
      case PGOpType::CliffordRot: {
        // Check all measures encountered so far can commute through
        for (const PGOp_ptr& m : measures) {
          if (!m->commutes_with(*pgp))
            throw PGError(
                "In legacy synthesis, cannot commute " + m->get_name() +
                " through " + pgp->get_name());
        }
        // Add a rotation to the end of the list
        rotations.push_back(pgp);
        break;
      }
      case PGOpType::Measure: {
        // Add to measures list
        measures.push_back(pgp);
        break;
      }
      default: {
        // Cannot synthesise other vertex kinds using legacy methods
        throw PGError(
            "Cannot synthesise PGOp using legacy synthesis: " +
            pgp->get_name());
      }
    }
  }
  std::map<Qubit, Bit> measure_map;
  std::set<Bit> bits_used;
  for (const PGOp_ptr& m : measures) {
    //  Push through output tableau to give measurement at the end
    const PGMeasure& pgm = dynamic_cast<PGMeasure&>(*m);
    SpPauliStabiliser paulis = pgm.get_tensor();
    if (out_tab) paulis = out_tab->get_row_product(paulis);
    paulis.compress();
    Bit target = pgm.get_target();
    // Assert measurement is Z on a single qubit
    if (paulis.size() != 1 || paulis.string.begin()->second != Pauli::Z)
      throw PGError(
          "In legacy synthesis, an end-of-circuit measurement is not a simple "
          "Z measurement");
    // Assert qubits and bits in measurements are disjoint
    bool new_qubit =
        measure_map.insert({paulis.string.begin()->first, target}).second;
    bool new_bit = bits_used.insert(target).second;
    if (!new_qubit || !new_bit)
      throw PGError(
          "In legacy synthesis, measurements cannot be performed "
          "simultaneously");
  }
  return {rotations, measure_map};
}

SpSymPauliTensor gadget_from_rotation(const PGOp_ptr& pgop) {
  switch (pgop->get_type()) {
    case PGOpType::Rotation: {
      PGRotation& r = dynamic_cast<PGRotation&>(*pgop);
      const SpPauliStabiliser& pauli = r.get_tensor();
      const Expr& angle = r.get_angle();
      return SpSymPauliTensor(pauli) * SpSymPauliTensor({}, angle);
    }
    case PGOpType::CliffordRot: {
      PGCliffordRot& r = dynamic_cast<PGCliffordRot&>(*pgop);
      const SpPauliStabiliser& pauli = r.get_tensor();
      unsigned angle = r.get_angle();
      return SpSymPauliTensor(pauli) * SpSymPauliTensor({}, angle * 0.5);
    }
    default: {
      // Only Rotation and CliffordRot are identified as rotations by
      // rotations_and_end_measures
      TKET_ASSERT(false);
    }
  }
}

void append_rotations_to_circuit_individually(
    Circuit& circ, const std::list<PGOp_ptr>& rotations,
    CXConfigType cx_config) {
  for (const PGOp_ptr& pgop : rotations) {
    append_single_pauli_gadget_as_pauli_exp_box(
        circ, gadget_from_rotation(pgop), cx_config);
  }
}

void append_rotations_to_circuit_pairwise(
    Circuit& circ, const std::list<PGOp_ptr>& rotations,
    CXConfigType cx_config) {
  auto it = rotations.begin();
  while (it != rotations.end()) {
    SpSymPauliTensor gadget0 = gadget_from_rotation(*it);
    ++it;
    if (it == rotations.end()) {
      append_single_pauli_gadget_as_pauli_exp_box(circ, gadget0, cx_config);
    } else {
      SpSymPauliTensor gadget1 = gadget_from_rotation(*it);
      ++it;
      // The new PauliGraph does not automatically merge rotations on
      // construction, leaving this for an explicit rewrite. However, the
      // diagonalisation in pairwise synthesis fails if identical strings are
      // provided. If this is the case, merge them into one gadget here.
      if (SpPauliString(gadget0) == SpPauliString(gadget1)) {
        gadget0.coeff += gadget1.coeff;
        append_single_pauli_gadget_as_pauli_exp_box(circ, gadget0, cx_config);
      } else {
        append_pauli_gadget_pair_as_box(circ, gadget0, gadget1, cx_config);
      }
    }
  }
}

void append_rotations_to_circuit_setwise(
    Circuit& circ, const std::list<PGOp_ptr>& rotations,
    CXConfigType cx_config) {
  auto it = rotations.begin();
  while (it != rotations.end()) {
    std::map<SpPauliString, SpSymPauliTensor> gadget_map;
    SpSymPauliTensor gadget = gadget_from_rotation(*it);
    gadget_map.insert({gadget.string, gadget});
    ++it;
    while (it != rotations.end()) {
      SpSymPauliTensor gadget = gadget_from_rotation(*it);
      bool commutes_with_all = true;
      for (const std::pair<const SpPauliString, SpSymPauliTensor>& g :
           gadget_map) {
        if (!gadget.commutes_with(g.first)) {
          commutes_with_all = false;
          break;
        }
      }
      if (!commutes_with_all) break;
      // Merge any gadgets with identical strings which may not have been merged
      // by explicit rewrites on the PauliGraph.
      auto inserted = gadget_map.insert({gadget.string, gadget});
      if (!inserted.second) inserted.first->second.coeff += gadget.coeff;
      ++it;
    }
    if (gadget_map.size() == 1) {
      SpSymPauliTensor g = gadget_map.begin()->second;
      append_single_pauli_gadget_as_pauli_exp_box(circ, g, cx_config);
    } else if (gadget_map.size() == 2) {
      SpSymPauliTensor g0 = gadget_map.begin()->second;
      SpSymPauliTensor g1 = (++gadget_map.begin())->second;
      append_pauli_gadget_pair_as_box(circ, g0, g1, cx_config);
    } else {
      std::list<SpSymPauliTensor> gadget_list;
      for (const std::pair<const SpPauliString, SpSymPauliTensor>& g :
           gadget_map)
        gadget_list.push_back(g.second);
      append_commuting_pauli_gadget_set_as_box(circ, gadget_list, cx_config);
    }
  }
}

Circuit pauli_graph3_to_circuit_legacy(
    const pg::PauliGraph& pg, CXConfigType cx_config,
    Transforms::PauliSynthStrat strat) {
  // Assert incoming and outgoing tableau is unitary and get UnitaryTableaus
  const std::set<Qubit>& qubits = pg.get_qubits();
  const std::set<Bit>& bits = pg.get_bits();
  Circuit circ(
      qubit_vector_t{qubits.begin(), qubits.end()},
      bit_vector_t{bits.begin(), bits.end()});
  // Synthesise input tableau
  std::optional<PGVert> itab_v = pg.get_input_tableau();
  if (itab_v) {
    PGOp_ptr pgop = pg.get_vertex_PGOp_ptr(*itab_v);
    PGInputTableau& tab_op = dynamic_cast<PGInputTableau&>(*pgop);
    ChoiMixTableau cmtab = tab_op.to_cm_tableau();
    UnitaryTableau in_utab = cm_tableau_to_unitary_tableau(cmtab);
    Circuit in_cliff_circuit = unitary_tableau_to_circuit(in_utab);
    ChoiAPState target = tab_op.get_apstate();
    target.canonical_column_order();
    target.normal_form();
    ChoiAPState recreated = circuit_to_choi_apstate(in_cliff_circuit);
    recreated.canonical_column_order();
    recreated.normal_form();
    in_cliff_circuit.add_phase(target.ap_.phase - recreated.ap_.phase);
    circ.append(in_cliff_circuit);
  }
  // Assert form is legacy-compatible and push measures to the end
  std::optional<UnitaryTableau> out_utab = std::nullopt;
  std::optional<PGVert> otab_v = pg.get_output_tableau();
  if (otab_v) {
    PGOp_ptr pgop = pg.get_vertex_PGOp_ptr(*otab_v);
    PGOutputTableau& tab_op = dynamic_cast<PGOutputTableau&>(*pgop);
    ChoiMixTableau cmtab = tab_op.to_cm_tableau();
    out_utab = cm_tableau_to_unitary_tableau(cmtab);
  }
  std::list<PGOp_ptr> rotations;
  std::map<Qubit, Bit> end_measures;
  std::tie(rotations, end_measures) = rotations_and_end_measures(pg, out_utab);
  // Synthesise rotations in order
  switch (strat) {
    case Transforms::PauliSynthStrat::Individual: {
      append_rotations_to_circuit_individually(circ, rotations, cx_config);
      break;
    }
    case Transforms::PauliSynthStrat::Pairwise: {
      append_rotations_to_circuit_pairwise(circ, rotations, cx_config);
      break;
    }
    case Transforms::PauliSynthStrat::Sets: {
      append_rotations_to_circuit_setwise(circ, rotations, cx_config);
      break;
    }
    default: {
      TKET_ASSERT(false);
    }
  }
  // Synthesise output tableau
  if (out_utab) {
    PGOp_ptr pgop = pg.get_vertex_PGOp_ptr(*otab_v);
    PGOutputTableau& tab_op = dynamic_cast<PGOutputTableau&>(*pgop);
    Circuit out_cliff_circuit = unitary_tableau_to_circuit(*out_utab);
    ChoiAPState target = tab_op.get_apstate();
    target.canonical_column_order();
    target.normal_form();
    ChoiAPState recreated = circuit_to_choi_apstate(out_cliff_circuit);
    recreated.canonical_column_order();
    recreated.normal_form();
    out_cliff_circuit.add_phase(target.ap_.phase - recreated.ap_.phase);
    circ.append(out_cliff_circuit);
  }
  // Add measures
  for (const std::pair<const Qubit, Bit>& m : end_measures) {
    circ.add_measure(m.first, m.second);
  }
  return circ;
}

Circuit pauli_graph3_to_pauli_exp_box_circuit_individually(
    const pg::PauliGraph& pg, CXConfigType cx_config) {
  return pauli_graph3_to_circuit_legacy(
      pg, cx_config, Transforms::PauliSynthStrat::Individual);
}

Circuit pauli_graph3_to_pauli_exp_box_circuit_pairwise(
    const pg::PauliGraph& pg, CXConfigType cx_config) {
  return pauli_graph3_to_circuit_legacy(
      pg, cx_config, Transforms::PauliSynthStrat::Pairwise);
}

Circuit pauli_graph3_to_pauli_exp_box_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config) {
  return pauli_graph3_to_circuit_legacy(
      pg, cx_config, Transforms::PauliSynthStrat::Sets);
}

/*******************************************************************************
 * PGOpCommutingSetBox Implementation
 ******************************************************************************/

PGOpCommutingSetBox::PGOpCommutingSetBox(
    const std::vector<PGOp_ptr>& pgops,
    const boost::bimap<Qubit, unsigned>& qubit_indices,
    const boost::bimap<Bit, unsigned>& bit_indices, CXConfigType cx_config)
    : Box(OpType::PGOpCommutingSetBox),
      pgops_(),
      qubit_indices_(qubit_indices),
      bit_indices_(bit_indices),
      cx_config_(cx_config) {
  for (auto it = pgops.begin(); it != pgops.end(); ++it) {
    // Check that all PGOps commute
    auto it2 = it;
    ++it2;
    while (it2 != pgops.end()) {
      if (!(*it)->commutes_with(*(*it2)))
        throw PGError(
            "The PGOps within a PGCommutingSetBox must all commute with each "
            "other");
      ++it2;
    }
    // Deep copy PGOps since PGOp_ptr is not const
    pgops_.push_back((*it)->clone());
    // Check every qubit is in qubit_indices
    for (unsigned port = 0; port < (*it)->n_paulis(); ++port) {
      const SpPauliStabiliser& pauli = (*it)->port(port);
      for (const std::pair<const Qubit, Pauli>& qp : pauli.string) {
        if (qubit_indices_.left.find(qp.first) == qubit_indices_.left.end())
          throw PGError(
              "All qubits in a PGCommutingSetBox need to be assigned indices");
      }
    }
    for (const Bit& b : (*it)->read_bits()) {
      if (bit_indices_.left.find(b) == bit_indices_.left.end())
        throw PGError(
            "All bits in a PGCommutingSetBox need to be assigned indices");
    }
    for (const Bit& b : (*it)->write_bits()) {
      if (bit_indices_.left.find(b) == bit_indices_.left.end())
        throw PGError(
            "All bits in a PGCommutingSetBox need to be assigned indices");
    }
  }
  unsigned n_qubits = qubit_indices_.size();
  for (const auto& pair : qubit_indices_) {
    if (pair.right >= n_qubits)
      throw PGError(
          "Cannot create PGCommutingSetBox: qubit indices do not span the "
          "range [0..n-1]");
  }
  unsigned n_bits = bit_indices_.size();
  for (const auto& pair : bit_indices_) {
    if (pair.right >= n_bits)
      throw PGError(
          "Cannot create PGCommutingSetBox: bit indices do not span the range "
          "[0..n-1]");
  }
  signature_ = op_signature_t(n_qubits, EdgeType::Quantum);
  op_signature_t bit_sig(n_bits, EdgeType::Classical);
  signature_.insert(signature_.end(), bit_sig.begin(), bit_sig.end());
}

PGOpCommutingSetBox::PGOpCommutingSetBox(const PGOpCommutingSetBox& other)
    : Box(other),
      pgops_(),
      qubit_indices_(other.qubit_indices_),
      bit_indices_(other.bit_indices_),
      cx_config_(other.cx_config_) {
  // Deep copy PGOps since PGOp_ptr is not const
  for (auto it = other.pgops_.begin(); it != pgops_.end(); ++it)
    pgops_.push_back((*it)->clone());
}

PGOpCommutingSetBox::PGOpCommutingSetBox() : PGOpCommutingSetBox({}) {}

bool PGOpCommutingSetBox::is_clifford() const {
  return std::all_of(pgops_.begin(), pgops_.end(), [](const PGOp_ptr& pgop) {
    return pgop->get_type() == PGOpType::CliffordRot;
  });
}

SymSet PGOpCommutingSetBox::free_symbols() const {
  SymSet sset;
  for (const PGOp_ptr& pgop : pgops_) {
    SymSet op_sset = pgop->free_symbols();
    sset.insert(op_sset.begin(), op_sset.end());
  }
  return sset;
}

Op_ptr PGOpCommutingSetBox::dagger() const {
  throw PGError("Dagger of PGOps is not yet implemented.");
}

Op_ptr PGOpCommutingSetBox::transpose() const {
  throw PGError("Transpose of PGOps is not yet implemented.");
}

Op_ptr PGOpCommutingSetBox::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  std::vector<PGOp_ptr> sub_ops;
  for (const PGOp_ptr& pgop : pgops_) {
    PGOp_ptr sub_op = pgop->symbol_substitution(sub_map);
    if (sub_op) {
      sub_ops.push_back(sub_op);
    } else {
      // Deep copy of PGOp
      sub_ops.push_back(pgop->clone());
    }
  }
  return std::make_shared<PGOpCommutingSetBox>(
      sub_ops, qubit_indices_, bit_indices_, cx_config_);
}

/**
 * GraySynthRecord describes the arguments to a call of the GraySynth algorithm,
 * allowing us to unfold the recursion into a loop by adding GraySynthRecords to
 * a stack in place of each recursive call.
 *
 * We have to generalise from phase polynomials to generic PGOps indexed by
 * diagonal Pauli strings. We represent each goal string by the set of Qubits on
 * which it acts as Z. Having more than just rotations, we cannot guarantee that
 * all PGOps indexed by the same goal string can be merged into one, so once the
 * target string is reduced to a single qubit we synthesise all PGOps from a
 * list.
 */
struct GraySynthRecord {
  std::map<std::set<Qubit>, std::list<PGOp_ptr>> terms;
  std::set<Qubit> remaining_qubits;
  std::optional<Qubit> target;
};

void PGOpCommutingSetBox::generate_circuit() const {
  /**
   * Each PGOp contains a number of active Pauli strings, generating the
   * subspace on which it acts non-trivially. The PGOps in the box must all
   * mutually commute, so each active Pauli must commute with all others from
   * other PGOps, but need not commute with all active Paulis in the same PGOp
   * (e.g. the active subspace of a Reset is generated by Z and X on the same
   * qubit).
   *
   * PGOp::pauli_signature() splits the active Pauli strings of a PGOp into
   * pairs of strings which anti-commute with each other but commute with all
   * other active Paulis, and the remaining strings which commute with every
   * active Pauli. For each anti-commuting pair, we can reduce them to Z and X
   * on a single qubit, at which point we know (by commutation) that no other
   * PGOp in the set acts on that qubit. These qubits then just become ancillas
   * used at the point of synthesising that particular PGOp, and can otherwise
   * be ignored.
   *
   * This leaves the commuting Pauli strings, which can all be mutually
   * diagonalised.
   *
   * If a PGOp has zero such commuting strings (e.g. reset or Box PGOps), we can
   * synthesise it at any point since it acts on a completely distinct set of
   * qubits from all others in the box.
   *
   * For the PGOps that contain exactly one commuting string, we can treat this
   * string like a term in a phase polynomial and feed it into GraySynth. At the
   * point where the term is reduced to Z on a single qubit, we can synthesise
   * the PGOp, factoring in any ancillas we prepared earlier.
   *
   * If a PGOp has multiple commuting strings, synthesis requires these to be
   * reduced to Z on distinct qubits at the same time. GraySynth makes no
   * guarantees about any qubits other than the immediate target, so it is hard
   * to fit this in to GraySynth in a sensible way. For simplicity, we will also
   * handle these PGOps separately from the GraySynth algorithm, alongside those
   * with zero commuting strings. This is likely to yield suboptimal solutions,
   * but these examples are not very common - e.g. controlled-operations and
   * multiplexors, which can theoretically be reduced to a large collection of
   * simple rotations before synthesis anyway.
   */

  // inner_circ will store the circuit conjugated by the diagonalisation
  // Cliffords
  std::set<Qubit> qubits;
  Circuit inner_circ;
  for (const auto& pair : qubit_indices_) {
    qubits.insert(pair.left);
    inner_circ.add_qubit(pair.left);
  }
  for (const auto& pair : bit_indices_) {
    inner_circ.add_bit(pair.left);
  }

  /**
   * Partition PGOps into those with one commuting term (a suitable target for
   * GraySynth), versus those with zero/multiple which will be synthesised
   * individually (individual synthesis for those with zero is optimal as the
   * Clifford conjugation is guaranteed to leave it on disjoint qubits from all
   * other PGOps).
   *
   * At the same time, we form the rows for the diagonalisation problem, split
   * into `ac_pair_rows` (to reduce anti-commuting pairs down to individual
   * qubits, from which they can be treated as ancilla qubits at the point of
   * synthesis during GraySynth) and `c_rows` (to diagonalise any remaining
   * active Paulis, making them viable parity terms for GraySynth).
   */
  std::list<ChoiMixTableau::row_tensor_t> ac_pair_rows;
  std::list<ChoiMixTableau::row_tensor_t> c_rows;
  std::set<Qubit>::const_iterator qb_it = qubits.begin();
  std::list<std::pair<SpPauliStabiliser, PGOp_ptr>> single_diag_ops;
  std::list<PGOp_ptr> multi_diag_ops;
  for (const PGOp_ptr& pgop : pgops_) {
    PGOp_signature sig = pgop->pauli_signature();
    for (const std::pair<SpPauliStabiliser, SpPauliStabiliser>& ac_pair :
         sig.anti_comm_pairs) {
      /**
       * All of the anti-commuting pairs are guaranteed to appear uniquely
       * since the PGOps must commute.
       *
       * The choice of qubit name on the output space here (i.e. `qb_it`) is
       * irrelevant so long as they are independent for each pair of rows. This
       * is because `cm_tableau_to_unitary_extension_circuit` solves the
       * diagonalisation circuit up to a permutation over the outputs which is
       * returned, so the choice of qubit names used here will only affect the
       * final permutation. We just pick qubit names that already exist in the
       * circuit for convenience.
       */
      ac_pair_rows.push_back(
          {ac_pair.first, SpPauliStabiliser(*qb_it, Pauli::Z)});
      ac_pair_rows.push_back(
          {ac_pair.second, SpPauliStabiliser(*qb_it, Pauli::X)});
      ++qb_it;
    }
    for (const SpPauliStabiliser& c_pauli : sig.comm_set) {
      // However, the commuting rows may appear multiple times or be generated
      // by a subset; we will filter them by gaussian elimination and removing
      // empty rows
      c_rows.push_back({c_pauli, {}});
    }
    // Clone so we don't accidentally change the Box when we conjugate
    PGOp_ptr cl = pgop->clone();
    // Filter based on recorded commuting ports
    if (sig.comm_set.size() == 1)
      single_diag_ops.push_back({*sig.comm_set.begin(), cl});
    else
      multi_diag_ops.push_back(cl);
  }

  /**
   * cm_tableau_to_unitary_extension_circuit only guarantees the tableau rows UP
   * TO additional Zs for initialised/post-selected qubits. When we have
   * ac_pairs, these may end up not being reduced to individual qubits, but
   * strings (ZZZ..., XZZ...).
   *
   * Rather than use one diag_tab, build different ones for the ac_pairs and the
   * commuting terms. Solve the ac_pairs first, then update the commuting terms
   * based on the unitary implementation and solve them. Since the ac_pairs
   * tableau doesn't involve any initialisations or post-selections, those rows
   * will be guaranteed exactly in the unitary implementation, then by
   * commutation the commuting terms will exist on disjoint qubits so the second
   * tableau won't mess with these rows. This should not be a significant
   * structural change to the solution circuits since ChoiMixTableau synthesis
   * solves id channels first.
   */
  ChoiMixTableau ac_tab(ac_pair_rows);
  auto [diag_circ, ac_qb_map] =
      cm_tableau_to_unitary_extension_circuit(ac_tab, {}, {}, cx_config_);
  // Push the remaining c_rows through this circuit so that when we diagonalise
  // the new strings, the composite circuit diagonalises the original ones
  UnitaryTableau u_ac_tab = circuit_to_unitary_tableau(diag_circ);
  for (ChoiMixTableau::row_tensor_t& row : c_rows) {
    row.first = u_ac_tab.get_row_product(row.first);
  }
  ChoiMixTableau c_tab(c_rows);
  // The strings are not guaranteed to be linearly-independent, so pick a
  // generating set for the space to give a proper tableau
  c_tab.canonical_column_order();
  c_tab.gaussian_form();
  for (unsigned r = c_tab.get_n_rows(); r > 0; r--) {
    if (c_tab.get_row(r - 1) == ChoiMixTableau::row_tensor_t{{}, {}})
      c_tab.remove_row(r - 1);
    else
      break;
  }
  // Diagonalisation is encoded into ChoiMixTableau synthesis via the unitary
  // extension of a post-selection circuit. We need to provide available names
  // for the qubits over which the diagonal results will act
  std::vector<Qubit> remaining_qubit_names;
  remaining_qubit_names.insert(
      remaining_qubit_names.end(), qb_it, qubits.end());
  auto [c_circ, c_qb_map] = cm_tableau_to_unitary_extension_circuit(
      c_tab, {}, remaining_qubit_names, cx_config_);
  // Compose to give the full diagonalisation circuit which both reduces each
  // ac_pair to an independent qubit and diagonalises all commuting strings
  diag_circ.append(c_circ);

  // Build a UnitaryTableau for the conjugation circuit to map paulis through
  UnitaryTableau u_tab = circuit_to_unitary_tableau(diag_circ);
  // Set up initial GraySynth state
  GraySynthRecord init_rec{{}, qubits, std::nullopt};
  for (std::pair<SpPauliStabiliser, PGOp_ptr>& term : single_diag_ops) {
    SpPauliStabiliser diag_string = u_tab.get_row_product(term.first);
    std::set<Qubit> z_qubits;
    diag_string.compress();
    for (const std::pair<const Qubit, Pauli>& qp : diag_string.string)
      z_qubits.insert(qp.first);
    for (unsigned p = 0; p < term.second->n_paulis(); ++p)
      term.second->port(p) = u_tab.get_row_product(term.second->port(p));
    auto [it, inserted] = init_rec.terms.insert({z_qubits, {term.second}});
    if (!inserted) it->second.push_back(term.second);
  }
  for (PGOp_ptr& pgop : multi_diag_ops) {
    for (unsigned p = 0; p < pgop->n_paulis(); ++p)
      pgop->port(p) = u_tab.get_row_product(pgop->port(p));
  }

  // Apply GraySynth to single_diag_ops
  std::list<GraySynthRecord> Q;
  Q.push_front(init_rec);
  // Atab builds up the CX circuit. Once a PGOp is ready to be synthesised, we
  // update it once by just pushing it through the entire tableau
  UnitaryTableau Atab(qubit_vector_t{qubits.begin(), qubits.end()});
  while (!Q.empty()) {
    GraySynthRecord R = Q.front();
    Q.pop_front();

    if (R.terms.size() == 0)
      continue;
    else if (R.terms.size() == 1 && R.target) {
      // When there is only a single term, we apply CXs to reduce it down to the
      // target qubit
      Qubit target = *(R.target);
      auto& [pstr, pgops] = *(R.terms.begin());
      for (const Qubit& control : pstr) {
        if (control != target) {
          inner_circ.add_op<Qubit>(OpType::CX, {control, target});
          for (GraySynthRecord& GSR : Q) {
            std::map<std::set<Qubit>, std::list<PGOp_ptr>> old_terms =
                GSR.terms;
            GSR.terms = {};
            for (std::pair<std::set<Qubit>, std::list<PGOp_ptr>> term :
                 old_terms) {
              if (term.first.find(target) != term.first.end()) {
                auto [it, inserted] = term.first.insert(control);
                if (!inserted) term.first.erase(it);
              }
              GSR.terms.insert(term);
            }
          }
          Atab.apply_CX_at_end(control, target);
        }
      }
      // Synthesise each PGOp associated to this term
      for (const PGOp_ptr& pgop : pgops) {
        for (unsigned p = 0; p < pgop->n_paulis(); ++p)
          pgop->port(p) = Atab.get_row_product(pgop->port(p));
        IndividualSynth synth = synth_pgop(pgop);
        // synth.tableau should describe a qubit permutation, telling us which
        // Qubits to place synth.op on
        std::map<Qubit, Qubit> position_map;
        for (unsigned i = 0; i < synth.tableau.get_n_rows(); ++i) {
          ChoiMixTableau::row_tensor_t row = synth.tableau.get_row(i);
          row.first.compress();
          TKET_ASSERT(row.first.size() == 1);
          row.second.compress();
          TKET_ASSERT(row.second.size() == 1);
          position_map[row.second.string.begin()->first] =
              row.first.string.begin()->first;
        }
        for (UnitID& u : synth.args) {
          if (u.type() == UnitType::Qubit) u = position_map.at(Qubit(u));
        }
        inner_circ.add_op<UnitID>(synth.op, synth.args);
      }
    } else if (!R.remaining_qubits.empty()) {
      // Recursive case: find the best qubit to split on (with the greatest
      // difference between the strings containing it and the strings not
      // containing it)
      int max = -1;
      std::optional<Qubit> max_q = std::nullopt;
      for (const Qubit& q : R.remaining_qubits) {
        int num_ones = std::count_if(
            R.terms.begin(), R.terms.end(),
            [=](const std::pair<const std::set<Qubit>, std::list<PGOp_ptr>>&
                    term) { return term.first.find(q) != term.first.end(); });
        int num_zeros = R.terms.size() - num_ones;
        if (num_zeros > max || num_ones > max) {
          max = (num_zeros > num_ones) ? num_zeros : num_ones;
          max_q = q;
        }
      }
      R.remaining_qubits.erase(*max_q);
      // Partition the terms into the two sets and recurse on each one
      GraySynthRecord R0{{}, R.remaining_qubits, R.target};
      GraySynthRecord R1{{}, R.remaining_qubits, R.target ? R.target : max_q};
      for (const std::pair<const std::set<Qubit>, std::list<PGOp_ptr>>& term :
           R.terms) {
        if (term.first.find(*max_q) == term.first.end())
          R0.terms.insert(term);
        else
          R1.terms.insert(term);
      }
      Q.push_front(R1);
      Q.push_front(R0);
    }
  }
  // We know Atab is just a CX circuit, so it can be completely described by one
  // quadrant and solved by Gaussian elimination
  DiagMatrix m(
      Atab.tab_.zmat.block(qubits.size(), 0, qubits.size(), qubits.size()));
  CXMaker cxmaker(qubits.size(), false);
  m.gauss(cxmaker);
  Circuit cx_circ = cxmaker._circ.dagger();
  // cxmaker yields a circuit with just default qubit names; place according to
  // Atab qubits
  unit_map_t cx_umap;
  for (const auto& pair : Atab.qubits_) {
    cx_umap.insert({Qubit(pair.right), pair.left});
  }
  cx_circ.rename_units(cx_umap);
  inner_circ.append(cx_circ);

  // Synthesise multi_diag_ops individually
  for (const PGOp_ptr& pgop : multi_diag_ops)
    inner_circ.append(pgop_to_circuit(pgop));

  // Combine into one Circuit (We cannot use ConjugationBox as the inner circuit
  // may include classical data)
  Circuit circ = diag_circ;
  circ.append(inner_circ);
  circ.append(diag_circ.dagger());

  // Rename units into flat register to hook up with Box arguments correctly
  unit_map_t umap;
  for (const auto& pair : qubit_indices_) {
    umap.insert({pair.left, Qubit(pair.right)});
  }
  for (const auto& pair : bit_indices_) {
    umap.insert({pair.left, Bit(pair.right)});
  }
  circ.rename_units(umap);
  circ_ = std::make_shared<Circuit>(circ);
}

bool PGOpCommutingSetBox::is_equal(const Op& op_other) const {
  const PGOpCommutingSetBox& other =
      dynamic_cast<const PGOpCommutingSetBox&>(op_other);
  if (id_ == other.get_id()) return true;
  if (cx_config_ != other.cx_config_) return false;
  if (qubit_indices_ != other.qubit_indices_) return false;
  if (bit_indices_ != other.bit_indices_) return false;
  return std::equal(
      pgops_.begin(), pgops_.end(), other.pgops_.begin(), other.pgops_.end(),
      [](const PGOp_ptr& a, const PGOp_ptr& b) { return *a == *b; });
}

nlohmann::json PGOpCommutingSetBox::to_json(const Op_ptr&) {
  throw PGError("PGOp serialisation is not yet implemented");
}

Op_ptr PGOpCommutingSetBox::from_json(const nlohmann::json&) {
  throw PGError("PGOp deserialisation is not yet implemented");
}

REGISTER_OPFACTORY(PGOpCommutingSetBox, PGOpCommutingSetBox)

}  // namespace tket
