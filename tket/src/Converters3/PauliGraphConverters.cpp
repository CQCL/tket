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

#include "Circuit/Boxes.hpp"
#include "Clifford/ChoiMixTableau.hpp"
#include "Clifford/UnitaryTableau.hpp"
#include "Converters.hpp"
#include "Gate/Gate.hpp"

namespace tket {

namespace pg {

Op_ptr PGBox::get_op() const { return op_; }

const unit_vector_t& PGBox::get_args() const { return args_; }

PGBox::PGBox(
    const Op_ptr& op, const unit_vector_t& args,
    const std::vector<QubitPauliTensor>& paulis)
    : PGOp(PGOpType::Box), op_(op), args_(args), paulis_(paulis) {
  op_signature_t sig = op_->get_signature();
  unsigned nqs = 0;
  for (const EdgeType& et : sig) {
    if (et == EdgeType::Quantum) ++nqs;
  }
  if (paulis.size() != 2 * nqs)
    throw PGError(
        "Cannot create PGBox; number of QubitPauliTensors must match twice the "
        "number of qubits in the op");
  if (args.size() != sig.size())
    throw PGError(
        "Cannot create PGBox; number of arguments must match the signature of "
        "the op");
  // Could consider checking commutation properties of the paulis to ensure they
  // are in anticommuting pairs for each qubit
}

SymSet PGBox::free_symbols() const { return op_->free_symbols(); }

PGOp_ptr PGBox::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  Op_ptr new_inner = op_->symbol_substitution(sub_map);
  if (new_inner)
    return std::make_shared<PGBox>(new_inner, args_, paulis_);
  else
    return PGOp_ptr();
}

std::string PGBox::get_name(bool latex) const {
  std::stringstream str;
  str << op_->get_name(latex) << "(";
  for (const UnitID& u : args_) {
    str << u.repr() << ", ";
  }
  str << "\b\b)";
  return str.str();
}

bool PGBox::is_equal(const PGOp& op_other) const {
  const PGBox& other = dynamic_cast<const PGBox&>(op_other);
  return (args_ == other.args_) && (op_ == other.op_);
}

unsigned PGBox::n_paulis() const { return paulis_.size(); }

std::vector<QubitPauliTensor> PGBox::active_paulis() const {
  std::vector<QubitPauliTensor> active;
  for (const UnitID& u : args_) {
    if (u.type() == UnitType::Qubit) {
      Qubit q(u);
      active.push_back(QubitPauliTensor(q, Pauli::Z));
      active.push_back(QubitPauliTensor(q, Pauli::X));
    }
  }
  return active;
}

QubitPauliTensor& PGBox::port(unsigned p) {
  if (p >= paulis_.size())
    throw PGError(
        "Cannot dereference port " + std::to_string(p) +
        " on PGBox: " + this->get_name());
  return paulis_.at(p);
}

bit_vector_t PGBox::read_bits() const {
  op_signature_t sig = op_->get_signature();
  bit_vector_t read;
  for (unsigned i = 0; i < sig.size(); ++i) {
    if (sig.at(i) == EdgeType::Boolean) read.push_back(Bit(args_.at(i)));
  }
  return read;
}

bit_vector_t PGBox::write_bits() const {
  op_signature_t sig = op_->get_signature();
  bit_vector_t writes;
  for (unsigned i = 0; i < sig.size(); ++i) {
    if (sig.at(i) == EdgeType::Classical) writes.push_back(Bit(args_.at(i)));
  }
  return writes;
}

}  // namespace pg

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
          tab.get_row_product(QubitPauliTensor(q, Pauli::Y)), 2)};
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
      QubitPauliTensor zq = tab.get_zrow(q);
      return {
          std::make_shared<PGCliffordRot>(zq, 1),
          std::make_shared<PGCliffordRot>(tab.get_xrow(q), 1),
          std::make_shared<PGCliffordRot>(zq, 1)};
    }
    case OpType::CX: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      QubitPauliTensor zc = tab.get_zrow(c);
      QubitPauliTensor xt = tab.get_xrow(t);
      return {
          std::make_shared<PGCliffordRot>(zc, 3),
          std::make_shared<PGCliffordRot>(xt, 3),
          std::make_shared<PGCliffordRot>(zc * xt, 1)};
    }
    case OpType::CY: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      QubitPauliTensor zc = tab.get_zrow(c);
      QubitPauliTensor yt = tab.get_row_product(QubitPauliTensor(t, Pauli::Y));
      return {
          std::make_shared<PGCliffordRot>(zc, 3),
          std::make_shared<PGCliffordRot>(yt, 3),
          std::make_shared<PGCliffordRot>(zc * yt, 1)};
    }
    case OpType::CZ: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      QubitPauliTensor zc = tab.get_zrow(c);
      QubitPauliTensor zt = tab.get_zrow(t);
      return {
          std::make_shared<PGCliffordRot>(zc, 3),
          std::make_shared<PGCliffordRot>(zt, 3),
          std::make_shared<PGCliffordRot>(zc * zt, 1)};
    }
    case OpType::ZZMax: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      QubitPauliTensor zc = tab.get_zrow(c);
      QubitPauliTensor zt = tab.get_zrow(t);
      return {std::make_shared<PGCliffordRot>(zc * zt, 1)};
    }
    case OpType::SWAP: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      QubitPauliTensor zc = tab.get_zrow(c);
      QubitPauliTensor zt = tab.get_zrow(t);
      QubitPauliTensor xc = tab.get_xrow(c);
      QubitPauliTensor xt = tab.get_xrow(t);
      return {
          std::make_shared<PGCliffordRot>(zc * zt, 1),
          std::make_shared<PGCliffordRot>(xc * xt, 1),
          std::make_shared<PGCliffordRot>(-1. * zc * xc * zt * xt, 1)};
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
          tab.get_row_product(QubitPauliTensor(q, Pauli::Y)),
          g.get_params().at(0))};
    }
    case OpType::TK1: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor zq = tab.get_zrow(q);
      return {
          std::make_shared<PGRotation>(zq, g.get_params().at(0)),
          std::make_shared<PGRotation>(tab.get_xrow(q), g.get_params().at(1)),
          std::make_shared<PGRotation>(zq, g.get_params().at(2))};
    }
    case OpType::PhaseGadget: {
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliMap qpm;
      for (const UnitID& a : args) qpm.insert({Qubit(a), Pauli::Z});
      return {std::make_shared<PGRotation>(
          tab.get_row_product(QubitPauliTensor(qpm)), g.get_params().at(0))};
    }
    case OpType::ZZPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor z0 = tab.get_zrow(q0);
      QubitPauliTensor z1 = tab.get_zrow(q1);
      return {std::make_shared<PGRotation>(z0 * z1, g.get_params().at(0))};
    }
    case OpType::XXPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor x0 = tab.get_xrow(q0);
      QubitPauliTensor x1 = tab.get_xrow(q1);
      return {std::make_shared<PGRotation>(x0 * x1, g.get_params().at(0))};
    }
    case OpType::YYPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor yy = tab.get_row_product(
          QubitPauliTensor({{q0, Pauli::Y}, {q1, Pauli::Y}}));
      return {std::make_shared<PGRotation>(yy, g.get_params().at(0))};
    }
    case OpType::TK2: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor z0 = tab.get_zrow(q0);
      QubitPauliTensor z1 = tab.get_zrow(q1);
      QubitPauliTensor x0 = tab.get_xrow(q0);
      QubitPauliTensor x1 = tab.get_xrow(q1);
      return {
          std::make_shared<PGRotation>(x0 * x1, g.get_params().at(0)),
          std::make_shared<PGRotation>(
              -1. * z0 * x0 * z1 * x1, g.get_params().at(1)),
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
          tab.get_row_product(QubitPauliTensor(qpm)), box.get_phase())};
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
      std::vector<QubitPauliTensor> paulis;
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

PauliGraph circuit_to_pauli_graph3(const Circuit& circ) {
  PauliGraph res;
  ChoiMixTableau initial(circ.all_qubits());
  for (const Qubit& q : circ.created_qubits())
    initial.post_select(q, ChoiMixTableau::TableauSegment::Input);
  res.add_vertex_at_end(std::make_shared<PGInputTableau>(initial));
  UnitaryRevTableau final_u(circ.all_qubits());
  for (const Command& com : circ) {
    unit_vector_t args = com.get_args();
    for (const PGOp_ptr& pgop :
         op_to_pgops(com.get_op_ptr(), args, final_u, true))
      res.add_vertex_at_end(pgop);
  }
  std::list<ChoiMixTableau::row_tensor_t> final_rows;
  for (const Qubit& q : final_u.get_qubits()) {
    final_rows.push_back({final_u.get_zrow(q), QubitPauliTensor(q, Pauli::Z)});
    final_rows.push_back({final_u.get_xrow(q), QubitPauliTensor(q, Pauli::X)});
  }
  ChoiMixTableau final_cm(final_rows);
  for (const Qubit& q : circ.discarded_qubits()) final_cm.discard_qubit(q);
  res.add_vertex_at_end(std::make_shared<PGOutputTableau>(final_cm));
  return res;
}

}  // namespace tket
