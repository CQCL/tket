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

#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Clifford/ChoiMixTableau.hpp"
#include "tket/Clifford/UnitaryTableau.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"
#include "tket/Gate/Gate.hpp"
#include "tket/PauliGraph/ConjugatePauliFunctions.hpp"
#include "tket/PauliGraphRefactor/Converters.hpp"
#include "tket/Transformations/PauliOptimisation.hpp"

namespace tket {

using namespace pg;

std::vector<PGOp_ptr> op_to_pgops(
    const Op_ptr& op, const unit_vector_t& args, UnitaryRevTableau& tab,
    bool allow_tableau) {
  if (allow_tableau && is_clifford_type(op->get_type())) {
    qubit_vector_t qs;
    for (const UnitID& a : args) qs.push_back(Qubit(a));
    tab.apply_gate_at_end(op->get_type(), qs);
    return {};
  }
  switch (op->get_type()) {
    case OpType::Z: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(tab.get_zrow(q), 2)};
    }
    case OpType::X: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(tab.get_xrow(q), 2)};
    }
    case OpType::Y: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(
          tab.get_row_product(SpPauliStabiliser(q, Pauli::Y)), 2)};
    }
    case OpType::S: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(tab.get_zrow(q), 1)};
    }
    case OpType::Sdg: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(tab.get_zrow(q), 3)};
    }
    case OpType::V: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(tab.get_xrow(q), 1)};
    }
    case OpType::Vdg: {
      Qubit q(args.front());
      return {std::make_shared<PGCliffordRot>(tab.get_xrow(q), 3)};
    }
    case OpType::H: {
      Qubit q(args.front());
      SpPauliStabiliser zq = tab.get_zrow(q);
      return {
          std::make_shared<PGCliffordRot>(zq, 1),
          std::make_shared<PGCliffordRot>(tab.get_xrow(q), 1),
          std::make_shared<PGCliffordRot>(zq, 1)};
    }
    case OpType::CX: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser xt = tab.get_xrow(t);
      return {
          std::make_shared<PGCliffordRot>(zc, 3),
          std::make_shared<PGCliffordRot>(xt, 3),
          std::make_shared<PGCliffordRot>(zc * xt, 1)};
    }
    case OpType::CY: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser yt =
          tab.get_row_product(SpPauliStabiliser(t, Pauli::Y));
      return {
          std::make_shared<PGCliffordRot>(zc, 3),
          std::make_shared<PGCliffordRot>(yt, 3),
          std::make_shared<PGCliffordRot>(zc * yt, 1)};
    }
    case OpType::CZ: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser zt = tab.get_zrow(t);
      return {
          std::make_shared<PGCliffordRot>(zc, 3),
          std::make_shared<PGCliffordRot>(zt, 3),
          std::make_shared<PGCliffordRot>(zc * zt, 1)};
    }
    case OpType::ZZMax: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      SpPauliStabiliser zc = tab.get_zrow(c);
      SpPauliStabiliser zt = tab.get_zrow(t);
      return {std::make_shared<PGCliffordRot>(zc * zt, 1)};
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
          std::make_shared<PGCliffordRot>(zc * zt, 1),
          std::make_shared<PGCliffordRot>(xc * xt, 1),
          std::make_shared<PGCliffordRot>(yy, 1)};
    }
    case OpType::T: {
      Qubit q(args.front());
      return {std::make_shared<PGRotation>(tab.get_zrow(q), 0.25)};
    }
    case OpType::Tdg: {
      Qubit q(args.front());
      return {std::make_shared<PGRotation>(tab.get_zrow(q), -0.25)};
    }
    case OpType::Rz: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {
          std::make_shared<PGRotation>(tab.get_zrow(q), g.get_params().at(0))};
    }
    case OpType::Rx: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {
          std::make_shared<PGRotation>(tab.get_xrow(q), g.get_params().at(0))};
    }
    case OpType::Ry: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {std::make_shared<PGRotation>(
          tab.get_row_product(SpPauliStabiliser(q, Pauli::Y)),
          g.get_params().at(0))};
    }
    case OpType::TK1: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser zq = tab.get_zrow(q);
      return {
          std::make_shared<PGRotation>(zq, g.get_params().at(2)),
          std::make_shared<PGRotation>(tab.get_xrow(q), g.get_params().at(1)),
          std::make_shared<PGRotation>(zq, g.get_params().at(0))};
    }
    case OpType::PhaseGadget: {
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliMap qpm;
      for (const UnitID& a : args) qpm.insert({Qubit(a), Pauli::Z});
      return {std::make_shared<PGRotation>(
          tab.get_row_product(SpPauliStabiliser(qpm)), g.get_params().at(0))};
    }
    case OpType::ZZPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser z0 = tab.get_zrow(q0);
      SpPauliStabiliser z1 = tab.get_zrow(q1);
      return {std::make_shared<PGRotation>(z0 * z1, g.get_params().at(0))};
    }
    case OpType::XXPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser x0 = tab.get_xrow(q0);
      SpPauliStabiliser x1 = tab.get_xrow(q1);
      return {std::make_shared<PGRotation>(x0 * x1, g.get_params().at(0))};
    }
    case OpType::YYPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      SpPauliStabiliser yy = tab.get_row_product(
          SpPauliStabiliser({{q0, Pauli::Y}, {q1, Pauli::Y}}));
      return {std::make_shared<PGRotation>(yy, g.get_params().at(0))};
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
          std::make_shared<PGRotation>(x0 * x1, g.get_params().at(0)),
          std::make_shared<PGRotation>(yy, g.get_params().at(1)),
          std::make_shared<PGRotation>(z0 * z1, g.get_params().at(2))};
    }
    case OpType::Measure: {
      return {std::make_shared<PGMeasure>(
          tab.get_zrow(Qubit(args.at(0))), Bit(args.at(1)))};
    }
    case OpType::Collapse: {
      return {
          std::make_shared<PGDecoherence>(tab.get_zrow(Qubit(args.front())))};
    }
    case OpType::Reset: {
      Qubit q(args.front());
      return {std::make_shared<PGReset>(tab.get_zrow(q), tab.get_xrow(q))};
    }
    case OpType::PauliExpBox: {
      const PauliExpBox& box = dynamic_cast<const PauliExpBox&>(*op);
      const std::vector<Pauli>& paulis = box.get_paulis();
      QubitPauliMap qpm;
      for (unsigned i = 0; i < args.size(); ++i)
        qpm.insert({Qubit(args.at(i)), paulis.at(i)});
      return {std::make_shared<PGRotation>(
          tab.get_row_product(SpPauliStabiliser(qpm)), box.get_phase())};
    }
    case OpType::QControlBox: {
      const QControlBox& box = dynamic_cast<const QControlBox&>(*op);
      unsigned n_controls = box.get_n_controls();
      std::vector<bool> control_state = box.get_control_state();
      unit_vector_t inner_args;
      for (unsigned i = n_controls; i < args.size(); ++i)
        inner_args.push_back(args.at(i));
      std::vector<PGOp_ptr> inner_ops =
          op_to_pgops(box.get_op(), inner_args, tab, false);
      std::vector<SpPauliStabiliser> control_paulis;
      for (unsigned i = 0; i < n_controls; ++i)
        control_paulis.push_back(tab.get_zrow(Qubit(args.at(i))));
      std::vector<PGOp_ptr> ret;
      for (const PGOp_ptr& inn_op : inner_ops)
        ret.push_back(std::make_shared<PGQControl>(
            inn_op, control_paulis, control_state));
      return ret;
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
      return {std::make_shared<PGMultiplexedTensoredBox>(
          op_map, control_paulis, target_paulis)};
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
          return {std::make_shared<PGMultiplexedRotation>(
              angle_map, control_paulis,
              tab.get_xrow(Qubit(args.at(n_controls))))};
        }
        case OpType::Ry: {
          return {std::make_shared<PGMultiplexedRotation>(
              angle_map, control_paulis,
              tab.get_row_product(
                  SpPauliStabiliser(Qubit(args.at(n_controls)), Pauli::Y)))};
        }
        case OpType::Rz: {
          return {std::make_shared<PGMultiplexedRotation>(
              angle_map, control_paulis,
              tab.get_zrow(Qubit(args.at(n_controls))))};
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
      return {std::make_shared<PGMultiplexedTensoredBox>(
          op_map, control_paulis, target_paulis)};
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
      return {std::make_shared<PGMultiplexedTensoredBox>(
          box.get_op_map(), control_paulis, target_paulis)};
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
      return ops;
    }
    case OpType::Conditional: {
      const Conditional& cond = dynamic_cast<const Conditional&>(*op);
      bit_vector_t cond_bits;
      unit_vector_t inner_args;
      for (unsigned i = 0; i < cond.get_width(); ++i)
        cond_bits.push_back(Bit(args.at(i)));
      for (unsigned i = cond.get_width(); i < args.size(); ++i)
        inner_args.push_back(args.at(i));
      std::vector<PGOp_ptr> inner_ops =
          op_to_pgops(cond.get_op(), inner_args, tab, false);
      std::vector<PGOp_ptr> ret;
      for (const PGOp_ptr& inn_op : inner_ops)
        ret.push_back(std::make_shared<PGConditional>(
            inn_op, cond_bits, cond.get_value()));
      return ret;
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
      return {std::make_shared<PGBox>(op, args, paulis)};
    }
  }
}

pg::PauliGraph circuit_to_pauli_graph3(
    const Circuit& circ, bool collect_cliffords) {
  pg::PauliGraph res;
  ChoiMixTableau initial(circ.all_qubits());
  for (const Qubit& q : circ.created_qubits())
    initial.post_select(q, ChoiMixTableau::TableauSegment::Input);
  res.add_vertex_at_end(std::make_shared<PGInputTableau>(initial));
  UnitaryRevTableau final_u(circ.all_qubits());
  for (const Command& com : circ) {
    unit_vector_t args = com.get_args();
    for (const PGOp_ptr& pgop :
         op_to_pgops(com.get_op_ptr(), args, final_u, collect_cliffords))
      res.add_vertex_at_end(pgop);
  }
  std::list<ChoiMixTableau::row_tensor_t> final_rows;
  for (const Qubit& q : final_u.get_qubits()) {
    final_rows.push_back({final_u.get_zrow(q), SpPauliStabiliser(q, Pauli::Z)});
    final_rows.push_back({final_u.get_xrow(q), SpPauliStabiliser(q, Pauli::X)});
  }
  ChoiMixTableau final_cm(final_rows);
  for (const Qubit& q : circ.discarded_qubits()) final_cm.discard_qubit(q);
  res.add_vertex_at_end(std::make_shared<PGOutputTableau>(final_cm));
  return res;
}

struct IndividualSynth {
  std::list<ChoiMixTableau::row_tensor_t> rows;
  Op_ptr op;
  unit_vector_t args;
};

enum class MultiplexorCase { Rotation, U2, TensoredU2, General };

IndividualSynth synth_pgop(const PGOp_ptr& pgop) {
  switch (pgop->get_type()) {
    case PGOpType::Rotation: {
      PGRotation& rot = dynamic_cast<PGRotation&>(*pgop);
      return IndividualSynth{
          {{rot.port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}},
          get_op_ptr(OpType::Rz, rot.get_angle()),
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
          {{crot.port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}},
          op,
          {Qubit(0)}};
    }
    case PGOpType::Measure: {
      PGMeasure& meas = dynamic_cast<PGMeasure&>(*pgop);
      return IndividualSynth{
          {{meas.port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}},
          get_op_ptr(OpType::Measure),
          {Qubit(0), meas.get_target()}};
    }
    case PGOpType::Decoherence: {
      return IndividualSynth{
          {{pgop->port(0), SpPauliStabiliser(Qubit(0), Pauli::Z)}},
          get_op_ptr(OpType::Collapse),
          {Qubit(0)}};
    }
    case PGOpType::Reset: {
      PGReset& reset = dynamic_cast<PGReset&>(*pgop);
      return IndividualSynth{
          {{reset.get_stab(), SpPauliStabiliser(Qubit(0), Pauli::Z)},
           {reset.get_destab(), SpPauliStabiliser(Qubit(0), Pauli::X)}},
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
          inner_synth.rows,
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
        synth.rows.push_back({cp, SpPauliStabiliser(q, Pauli::Z)});
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
        synth.rows.push_back(
            {control_paulis.at(qi), SpPauliStabiliser(q, Pauli::Z)});
        synth.args.push_back(q);
      }
      const std::vector<SpPauliStabiliser>& target_paulis =
          mptb.get_target_paulis();
      for (unsigned qi = 0; 2 * qi < target_paulis.size(); ++qi) {
        Qubit q(control_paulis.size() + qi);
        synth.rows.push_back(
            {target_paulis.at(2 * qi), SpPauliStabiliser(q, Pauli::Z)});
        synth.rows.push_back(
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
      PGMultiplexedRotation& mpr =
          dynamic_cast<PGMultiplexedRotation&>(*pgop);
      IndividualSynth synth;
      const std::vector<SpPauliStabiliser>& control_paulis =
          mpr.get_control_paulis();
      for (unsigned qi = 0; qi < control_paulis.size(); ++qi) {
        Qubit q(qi);
        synth.rows.push_back(
            {control_paulis.at(qi), SpPauliStabiliser(q, Pauli::Z)});
        synth.args.push_back(q);
      }
      synth.rows.push_back(
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
        res.rows.push_back(
            {box.port(2 * i), SpPauliStabiliser(Qubit(i), Pauli::Z)});
        res.rows.push_back(
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
          {
              {stab.get_stab(), SpPauliStabiliser(Qubit(0), Pauli::Z)},
              {stab.get_anc_z(), SpPauliStabiliser(Qubit(1), Pauli::Z)},
              {stab.get_anc_x(), SpPauliStabiliser(Qubit(1), Pauli::X)},
          },
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
  ChoiMixTableau diag_tab(synth_data.rows);
  auto [diag_circ, qb_map] = cm_tableau_to_unitary_extension_circuit(diag_tab);
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
    if (pgop->get_type() == PGOpType::InputTableau ||
        pgop->get_type() == PGOpType::OutputTableau) {
      ChoiMixTableau tab(0);
      if (pgop->get_type() == PGOpType::InputTableau) {
        PGInputTableau& tab_op = dynamic_cast<PGInputTableau&>(*pgop);
        tab = tab_op.to_cm_tableau();
      } else {
        PGOutputTableau& tab_op = dynamic_cast<PGOutputTableau&>(*pgop);
        tab = tab_op.to_cm_tableau();
      }
      std::pair<Circuit, qubit_map_t> tab_circ =
          cm_tableau_to_exact_circuit(tab, cx_config);
      qubit_map_t perm;
      for (const std::pair<const Qubit, Qubit>& p : tab_circ.second)
        perm.insert({p.second, p.first});
      tab_circ.first.permute_boundary_output(perm);
      circ.append(tab_circ.first);
    } else {
      circ.append(pgop_to_circuit(pgop));
    }
  }
  return circ;
}

Circuit pauli_graph3_to_circuit_sets(const pg::PauliGraph& pg, CXConfigType) {
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
    ChoiMixTableau cmtab = tab_op.to_cm_tableau();
    UnitaryTableau in_utab = cm_tableau_to_unitary_tableau(cmtab);
    Circuit in_cliff_circuit = unitary_tableau_to_circuit(in_utab);
    circ.append(in_cliff_circuit);
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

  // Synthesise each interior commuting set
  for (const std::list<PGOp_ptr>& cset : commuting_sets) {
    // TODO::Implement
    // TODO::I can somewhat see what to do when the active paulis of an op can
    // be partitioned into a set of pairs of anti-commuting paulis for exact
    // reduction to specific qubits, and up to one pauli for diagonalisation.
    // Will we ever have a kind of op that will have multiple strings that
    // require simultaneous diagonalisation and reduction to single qubits at
    // the same time e.g. a quantum-controlled unitary op? Graysynth targets one
    // string at a time, so how do we generalise it to target multiple strings
    // in parallel?

    // For each PGOp, split the active paulis into pairs of anti-commuting
    // paulis which will be reduced to a specific qubit, and others which will
    // just be diagonalised and then handled by graysynth
    std::list<std::pair<SpPauliStabiliser, SpPauliStabiliser>> ac_pairs;
    std::list<SpPauliStabiliser> to_diag;
    for (const PGOp_ptr& pgop : cset) {
      std::list<SpPauliStabiliser> commuting;
      for (const SpPauliStabiliser& p : pgop->active_paulis()) {
        bool anti_found = false;
        for (std::list<SpPauliStabiliser>::iterator it = commuting.begin();
             it != commuting.end(); ++it) {
          if (!p.commutes_with(*it)) {
            anti_found = true;
            ac_pairs.push_back({p, *it});
            commuting.erase(it);
            break;
          }
        }
        if (!anti_found) {
          commuting.push_back(p);
        }
      }
      to_diag.insert(to_diag.end(), commuting.begin(), commuting.end());
    }

    // Pick an independent generating set of the commuting paulis to be
    // diagonalised
    std::list<ChoiMixTableau::row_tensor_t> diag_tab_rows;
    std::set<Qubit>::const_iterator qb_it = qubits.begin();
    for (const std::pair<SpPauliStabiliser, SpPauliStabiliser>& pair :
         ac_pairs) {
      diag_tab_rows.push_back(
          {pair.first, SpPauliStabiliser(*qb_it, Pauli::Z)});
      diag_tab_rows.push_back(
          {pair.second, SpPauliStabiliser(*qb_it, Pauli::Z)});
      ++qb_it;
    }
    for (const SpPauliStabiliser& p : to_diag) {
      diag_tab_rows.push_back({p, {}});
    }
    ChoiMixTableau diag_tab(
        diag_tab_rows);  // TODO:: GAUSSIAN FORM AND REMOVE EMPTY ROWS
    // TODO:: RETURN TO OTHER SYNTHESIS METHODS TO INCORPORATE QControl and
    // Multiplexor

    // Obtain the conjugation circuit as the unitary extension of a
    // ChoiMixTableau

    // Build a UnitaryTableau for the conjugation circuit to map paulis through

    // Sort PGOps into groups whose diagonal components are identical (i.e. they
    // will be synthesised at the same point during GraySynth)

    // Conjugate

    // Append to the circuit
  }

  // Synthesise output tableau
  if (otab_v) {
    PGOp_ptr pgop = pg.get_vertex_PGOp_ptr(*otab_v);
    PGOutputTableau& tab_op = dynamic_cast<PGOutputTableau&>(*pgop);
    ChoiMixTableau cmtab = tab_op.to_cm_tableau();
    UnitaryTableau out_utab = cm_tableau_to_unitary_tableau(cmtab);
    Circuit out_cliff_circuit = unitary_tableau_to_circuit(out_utab);
    circ.append(out_cliff_circuit);
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
  std::list<PGOp_ptr> all_pgops = pg.pgop_sequence();
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
    Circuit out_cliff_circuit = unitary_tableau_to_circuit(*out_utab);
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

}  // namespace tket
