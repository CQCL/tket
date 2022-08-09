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

#include "Converters.hpp"
#include "Circuit/Boxes.hpp"
#include "Gate/Gate.hpp"

namespace tket {

namespace pg {

Op_ptr PGBox::get_op() const { return op_; }

const unit_vector_t& PGBox::get_args() const { return args_; }

PGBox::PGBox(const Op_ptr& op, const unit_vector_t& args) : PGOp(PGOpType::Box), op_(op), args_(args) {}

SymSet PGBox::free_symbols() const { return op_->free_symbols(); }

PGOp_ptr PGBox::symbol_substitution(const SymEngine::map_basic_basic& sub_map) const {
  Op_ptr new_inner = op_->symbol_substitution(sub_map);
  if (new_inner) return std::make_shared<const PGBox>(new_inner, args_);
  else return PGOp_ptr();
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

std::vector<PGOp_ptr> solve_qubit_with_gadgets(const QubitPauliTensor& zrow, const QubitPauliTensor& xrow, const Qubit& qb);

QubitPauliTensor stab_to_tensor(const PauliStabiliser& stab, const tableau_col_index_t& qbs) {
  QubitPauliMap qpm;
  for (unsigned i = 0; i < stab.string.size(); ++i) {
    Pauli p = stab.string.at(i);
    if (p != Pauli::I) qpm.insert({qbs.right.at(i), p});
  }
  return QubitPauliTensor(qpm, stab.coeff ? 1. : -1.);
}

PauliStabiliser tensor_to_stab(const QubitPauliTensor& ten, const tableau_col_index_t& qbs) {
  std::vector<Pauli> ps;
  for (unsigned i = 0; i < qbs.size(); ++i) {
    Qubit qb = qbs.right.at(i);
    ps.push_back(ten.string.get(qb));
  }
  return PauliStabiliser(ps, (ten.coeff == 1.));
}

std::vector<PGOp_ptr> op_to_pgops(const Op_ptr& op, const unit_vector_t& args, PauliGraph& pg, bool allow_tableau) {
  switch (op->get_type()) {
    case OpType::Z: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.phase_(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})) ^= true;
        return {};
      }
      else return {std::make_shared<const PGCliffordRot>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_), 2)};
    }
    case OpType::X: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.phase_(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})) ^= true;
        return {};
      }
      else return {std::make_shared<const PGCliffordRot>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_), 2)};
    }
    case OpType::Y: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.phase_(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})) ^= true;
        pg.final_.phase_(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})) ^= true;
        return {};
      }
      else {
        QubitPauliTensor zten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_);
        QubitPauliTensor xten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_);
        return {std::make_shared<const PGCliffordRot>(i_ * xten * zten, 2)};
      }
    }
    case OpType::S: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow}), i_);
        return {};
      }
      else return {std::make_shared<const PGCliffordRot>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_), 1)};
    }
    case OpType::Sdg: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow}), -i_);
        return {};
      }
      else return {std::make_shared<const PGCliffordRot>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_), 3)};
    }
    case OpType::V: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow}), i_);
        return {};
      }
      else return {std::make_shared<const PGCliffordRot>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_), 1)};
    }
    case OpType::Vdg: {
      Qubit q(args.front());
      if (allow_tableau) {
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow}), -i_);
        return {};
      }
      else return {std::make_shared<const PGCliffordRot>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_), 3)};
    }
    case OpType::H: {
      Qubit q(args.front());
      if (allow_tableau) {
        auto found = pg.final_rows_.left.find(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow});
        unsigned x = found->second;
        pg.final_rows_.left.erase(found);
        found = pg.final_rows_.left.find(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow});
        unsigned z = found->second;
        pg.final_rows_.left.erase(found);
        pg.final_rows_.insert({{q, TableauRowType::ZRow}, x});
        pg.final_rows_.insert({{q, TableauRowType::XRow}, z});
        return {};
      }
      else {
        PGOp_ptr s = op_to_pgops(get_op_ptr(OpType::S), args, pg, false).front();
        PGOp_ptr v = op_to_pgops(get_op_ptr(OpType::V), args, pg, false).front();
        return {s, v, s};
      }
    }
    case OpType::CX: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      if (allow_tableau) {
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::XRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::XRow}));
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::ZRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::ZRow}));
        return {};
      }
      else {
        QubitPauliTensor cten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::ZRow})), pg.final_cols_);
        QubitPauliTensor tten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::XRow})), pg.final_cols_);
        return {std::make_shared<const PGCliffordRot>(cten, 3), std::make_shared<const PGCliffordRot>(tten, 3), std::make_shared<const PGCliffordRot>(cten * tten, 1)};
      }
    }
    case OpType::CY: {
      Qubit c(args.at(0));
      unsigned cz = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::ZRow});
      unsigned cx = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::XRow});
      Qubit t(args.at(1));
      unsigned tz = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::ZRow});
      unsigned tx = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::XRow});
      if (allow_tableau) {
        pg.final_.row_mult(tz, tx, i_);
        pg.final_.row_mult(tx, cx);
        pg.final_.row_mult(cz, tz);
        pg.final_.row_mult(tz, tx, -i_);
        return {};
      }
      else {
        QubitPauliTensor cten = stab_to_tensor(pg.final_.get_pauli(cz), pg.final_cols_);
        QubitPauliTensor txten = stab_to_tensor(pg.final_.get_pauli(tx), pg.final_cols_);
        QubitPauliTensor tzten = stab_to_tensor(pg.final_.get_pauli(tz), pg.final_cols_);
        QubitPauliTensor tyten = i_ * txten * tzten;
        return {std::make_shared<const PGCliffordRot>(cten, 3), std::make_shared<const PGCliffordRot>(tyten, 3), std::make_shared<const PGCliffordRot>(cten * tyten, 1)};
      }
    }
    case OpType::CZ: {
      Qubit c(args.at(0));
      Qubit t(args.at(1));
      if (allow_tableau) {
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::ZRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::XRow}));
        pg.final_.row_mult(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::ZRow}), pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::XRow}));
        return {};
      }
      else {
        QubitPauliTensor cten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::ZRow})), pg.final_cols_);
        QubitPauliTensor tten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::ZRow})), pg.final_cols_);
        return {std::make_shared<const PGCliffordRot>(cten, 3), std::make_shared<const PGCliffordRot>(tten, 3), std::make_shared<const PGCliffordRot>(cten * tten, 1)};
      }
    }
    case OpType::ZZMax: {
      Qubit c(args.at(0));
      unsigned cz = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::ZRow});
      unsigned cx = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{c, TableauRowType::XRow});
      Qubit t(args.at(1));
      unsigned tz = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::ZRow});
      unsigned tx = pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{t, TableauRowType::XRow});
      if (allow_tableau) {
        pg.final_.row_mult(tz, cx);
        pg.final_.row_mult(cz, tx);
        pg.final_.row_mult(cz, cx, i_);
        pg.final_.row_mult(tz, tx, i_);
        return {};
      }
      else {
        QubitPauliTensor cten = stab_to_tensor(pg.final_.get_pauli(cz), pg.final_cols_);
        QubitPauliTensor tten = stab_to_tensor(pg.final_.get_pauli(tz), pg.final_cols_);
        return {std::make_shared<const PGCliffordRot>(cten * tten, 1)};
      }
    }
    case OpType::SWAP: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      if (allow_tableau) {
        auto found = pg.final_rows_.left.find(std::pair<Qubit, TableauRowType>{q0, TableauRowType::XRow});
        unsigned x0 = found->second;
        pg.final_rows_.left.erase(found);
        found = pg.final_rows_.left.find(std::pair<Qubit, TableauRowType>{q0, TableauRowType::ZRow});
        unsigned z0 = found->second;
        pg.final_rows_.left.erase(found);
        found = pg.final_rows_.left.find(std::pair<Qubit, TableauRowType>{q1, TableauRowType::XRow});
        unsigned x1 = found->second;
        pg.final_rows_.left.erase(found);
        found = pg.final_rows_.left.find(std::pair<Qubit, TableauRowType>{q1, TableauRowType::ZRow});
        unsigned z1 = found->second;
        pg.final_rows_.left.erase(found);
        pg.final_rows_.insert({{q0, TableauRowType::ZRow}, z1});
        pg.final_rows_.insert({{q0, TableauRowType::XRow}, x1});
        pg.final_rows_.insert({{q1, TableauRowType::ZRow}, z0});
        pg.final_rows_.insert({{q1, TableauRowType::XRow}, x0});
        return {};
      }
      else {
        QubitPauliTensor z0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::ZRow})), pg.final_cols_);
        QubitPauliTensor x0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::XRow})), pg.final_cols_);
        QubitPauliTensor z1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::ZRow})), pg.final_cols_);
        QubitPauliTensor x1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::XRow})), pg.final_cols_);
        return {std::make_shared<const PGCliffordRot>(z0 * z1, 1), std::make_shared<const PGCliffordRot>(x0 * x1, 1), std::make_shared<const PGCliffordRot>(-1. * z0 * x0 * z1 * x1, 1)};
      }
    }
    case OpType::Rz: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {std::make_shared<const PGRotation>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_), g.get_params().at(0))};
    }
    case OpType::Rx: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      return {std::make_shared<const PGRotation>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_), g.get_params().at(0))};
    }
    case OpType::Ry: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor zten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor xten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_);
      return {std::make_shared<const PGRotation>(i_ * xten * zten, g.get_params().at(0))};
    }
    case OpType::TK1: {
      Qubit q(args.front());
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor zten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor xten = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_);
      return {std::make_shared<const PGRotation>(zten, g.get_params().at(0)), std::make_shared<const PGRotation>(xten, g.get_params().at(1)), std::make_shared<const PGRotation>(zten, g.get_params().at(2))};
    }
    case OpType::PhaseGadget: {
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor comb{};
      for (const UnitID& a : args) {
        Qubit q(a);
        comb = comb * stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_);
      }
      return {std::make_shared<const PGRotation>(comb, g.get_params().at(0))};
    }
    case OpType::ZZPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor z0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor z1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::ZRow})), pg.final_cols_);
      return {std::make_shared<const PGRotation>(z0 * z1, g.get_params().at(0))};
    }
    case OpType::XXPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor x0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::XRow})), pg.final_cols_);
      QubitPauliTensor x1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::XRow})), pg.final_cols_);
      return {std::make_shared<const PGRotation>(x0 * x1, g.get_params().at(0))};
    }
    case OpType::YYPhase: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor z0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor z1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor x0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::XRow})), pg.final_cols_);
      QubitPauliTensor x1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::XRow})), pg.final_cols_);
      return {std::make_shared<const PGRotation>(-1. * z0 * x0 * z1 * x1, g.get_params().at(0))};
    }
    case OpType::TK2: {
      Qubit q0(args.at(0));
      Qubit q1(args.at(1));
      const Gate& g = dynamic_cast<const Gate&>(*op);
      QubitPauliTensor z0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor z1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::ZRow})), pg.final_cols_);
      QubitPauliTensor x0 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q0, TableauRowType::XRow})), pg.final_cols_);
      QubitPauliTensor x1 = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q1, TableauRowType::XRow})), pg.final_cols_);
      return {std::make_shared<const PGRotation>(x0 * x1, g.get_params().at(0)), std::make_shared<const PGRotation>(-1. * z0 * x0 * z1 * x1, g.get_params().at(1)), std::make_shared<const PGRotation>(z0 * z1, g.get_params().at(2))};
    }
    case OpType::Measure: {
      return {std::make_shared<const PGMeasure>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{args.at(0), TableauRowType::ZRow})), pg.final_cols_), Bit(args.at(1)))};
    }
    case OpType::Collapse: {
      return {std::make_shared<const PGDecoherence>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{args.front(), TableauRowType::ZRow})), pg.final_cols_))};
    }
    case OpType::Reset: {
      return {std::make_shared<const PGReset>(stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{args.front(), TableauRowType::ZRow})), pg.final_cols_), stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{args.front(), TableauRowType::XRow})), pg.final_cols_))};
    }
    case OpType::PauliExpBox: {
      const PauliExpBox& box = dynamic_cast<const PauliExpBox&>(*op);
      QubitPauliTensor comb{};
      const std::vector<Pauli>& paulis = box.get_paulis();
      for (unsigned i = 0; i < args.size(); ++i) {
        Qubit q(args.at(i));
        switch (paulis.at(i)) {
          case Pauli::I: {
            break;
          }
          case Pauli::X: {
            comb = comb * stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_);
            break;
          }
          case Pauli::Y: {
            comb = i_ *comb * stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::XRow})), pg.final_cols_) * stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_);
            break;
          }
          case Pauli::Z: {
            comb = comb * stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{q, TableauRowType::ZRow})), pg.final_cols_);
            break;
          }
        }
      }
      return {std::make_shared<const PGRotation>(comb, box.get_phase())};
    }
    case OpType::Conditional: {
      const Conditional& cond = dynamic_cast<const Conditional&>(*op);
      bit_vector_t cond_bits;
      unit_vector_t inner_args;
      for (unsigned i = 0; i < cond.get_width(); ++i) cond_bits.push_back(Bit(args.at(i)));
      for (unsigned i = cond.get_width(); i < args.size(); ++i) inner_args.push_back(args.at(i));
      std::vector<PGOp_ptr> inner_ops = op_to_pgops(cond.get_op(), inner_args, pg, false);
      std::vector<PGOp_ptr> ret;
      for (const PGOp_ptr& inn_op : inner_ops) ret.push_back(std::make_shared<const PGConditional>(inn_op, cond_bits, cond.get_value()));
      return ret;
    }
    default: {
      std::vector<PGOp_ptr> ret;
      for (const UnitID& uid : args) {
        if (uid.type() == UnitType::Qubit) {
          Qubit qb(uid);
          QubitPauliTensor zrow = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{qb, TableauRowType::ZRow})), pg.final_cols_);
          QubitPauliTensor xrow = stab_to_tensor(pg.final_.get_pauli(pg.final_rows_.left.at(std::pair<Qubit, TableauRowType>{qb, TableauRowType::XRow})), pg.final_cols_);
          std::vector<PGOp_ptr> solve_qb = solve_qubit_with_gadgets(zrow, xrow, qb);
          for (const PGOp_ptr& pgop : solve_qb) {
            const PGCliffordRot& cr = dynamic_cast<const PGCliffordRot&>(*pgop);
            pg.final_.apply_pauli_gadget(tensor_to_stab(cr.get_tensor(), pg.final_cols_), cr.get_angle());
            ret.push_back(pgop);
          }
        }
      }
      ret.push_back(std::make_shared<PGBox>(op, args));
      return ret;
    }
  }
}

/**
 * CLifford Unitary Synthesis by Pauli Gadgets:
 * 
 * For each qubit, it has Z row P and X row Q. The only relation we know is that P and Q must anti-commute and both commute with every other row for other qubits.
 * We can apply a small number of Clifford rotations to map P to Z and Q to X. The commutation relations mean that this will also make all other rows act as I on this qubit, completely solving it and removing it from the tableau.
 * We repeat this to solve as many of the qubits as we need.
 * The choice of gadgets to apply may be dependent on how P and Q act on our chosen qubit.
 * This will not give a minimal number to synthesise the full unitary, e.g. as given in https://arxiv.org/abs/2102.11380, but does give a simple and convenient way to solve only a few qubits as needed
 * 
 * P-I, Q-I: iPQY, QX (phase to solve for Z), PZ (phase to solve for X)
 * P-I, Q-X: iQZ, PQY (phase to solve for Z after next), Y (phase to solve for X)
 * P-I, Q-Y: iQX (phase to solve for X), iPQY (phase to solve for Z)
 * P-I, Q-Z: iQX (phase to solve for X), PQY (phase to solve for Z)
 * P-X, Q-I: iPZ (phase to solve for Z), PQY (phase to solve for X)
 * P-X, Q-X: iPZ (phase to solve for Z), iQY, Z (phase to solve for X)
 * P-X, Q-Y: iPZ (phase to solve for Z), iQX (phase to solve for X)
 * P-X, Q-Z: iQY, iPZ (phase to solve for Z), Z (phase to solve for X)
 * P-Y, Q-I: iPZ (phase to solve for Z), iPQY (phase to solve for X)
 * P-Y, Q-X: iPZ (phase to solve for Z), iQY, Z (phase to solve for X)
 * P-Y, Q-Y: iPZ (phase to solve for Z), iQX (phase to solve for X)
 * P-Y, Q-Z: iPZ (phase to solve for Z), PQY (phase to solve for X)
 * P-Z, Q-I: iPX, PQY (phase to solve for X after next), Y (phase to solve for Z)
 * P-Z, Q-X: PQZ, iQY, Y (phase to solve for Z), Z (phase to solve for X)
 * P-Z, Q-Y: iQX (phase to solve for X), iPY, X (phase to solve for Z)
 * P-Z, Q-Z: iQX (phase to solve for X), iPY, X (phase to solve for Z)
 */

std::pair<Complex, Complex> verify_solution_to_phase(QubitPauliTensor zrow, QubitPauliTensor xrow, const Qubit& qb, const std::vector<QubitPauliTensor>& tens) {
  for (const QubitPauliTensor& t : tens) {
    if (!zrow.commutes_with(t)) {
      zrow = i_ * zrow * t;
    }
    if (!xrow.commutes_with(t)) {
      xrow = i_ * xrow * t;
    }
  }
  if (zrow.string != QubitPauliString(qb, Pauli::Z)) throw std::logic_error("Tableau solving via gadgets did not reduce row to Z");
  if (xrow.string != QubitPauliString(qb, Pauli::X)) throw std::logic_error("Tableau solving via gadgets did not reduce row to X");
  return {zrow.coeff, xrow.coeff};
}

constexpr unsigned switch_pauli_pair(Pauli p0, Pauli p1) {
  return ((unsigned) p0 << 8) + (unsigned) p1;
}

std::vector<PGOp_ptr> solve_qubit_with_gadgets(const QubitPauliTensor& zrow, const QubitPauliTensor& xrow, const Qubit& qb) {
  Pauli qz = zrow.string.get(qb);
  Pauli qx = xrow.string.get(qb);
  switch (switch_pauli_pair(qz, qx)) {
    case switch_pauli_pair(Pauli::I, Pauli::I): {
      QubitPauliTensor t0 = zrow * xrow * QubitPauliTensor(qb, Pauli::Y, i_);
      QubitPauliTensor t1 = xrow * QubitPauliTensor(qb, Pauli::X);
      QubitPauliTensor t2 = zrow * QubitPauliTensor(qb, Pauli::Z);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1, t2});
      t1.coeff *= phases.first;
      t2.coeff *= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1), std::make_shared<const PGCliffordRot>(t2, 1)};
    }
    case switch_pauli_pair(Pauli::I, Pauli::X): {
      QubitPauliTensor t0 = xrow * QubitPauliTensor(qb, Pauli::Z, i_);
      QubitPauliTensor t1 = zrow * xrow * QubitPauliTensor(qb, Pauli::Y);
      QubitPauliTensor t2 = QubitPauliTensor(qb, Pauli::Y);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1, t2});
      t1.coeff *= phases.first * phases.second;
      t2.coeff *= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1), std::make_shared<const PGCliffordRot>(t2, 1)};
    }
    case switch_pauli_pair(Pauli::I, Pauli::Y):
    case switch_pauli_pair(Pauli::I, Pauli::Z): {
      QubitPauliTensor t0 = xrow * QubitPauliTensor(qb, Pauli::X, i_);
      QubitPauliTensor t1 = zrow * xrow * QubitPauliTensor(qb, Pauli::Y);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1});
      t0.coeff *= phases.second;
      t1.coeff /= phases.first;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1)};
    }
    case switch_pauli_pair(Pauli::X, Pauli::I):
    case switch_pauli_pair(Pauli::Y, Pauli::I):
    case switch_pauli_pair(Pauli::Y, Pauli::Z): {
      QubitPauliTensor t0 = zrow * QubitPauliTensor(qb, Pauli::Z, i_);
      QubitPauliTensor t1 = zrow * xrow * QubitPauliTensor(qb, Pauli::Y);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1});
      t0.coeff *= phases.first;
      t1.coeff /= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1)};
    }
    case switch_pauli_pair(Pauli::X, Pauli::X):
    case switch_pauli_pair(Pauli::Y, Pauli::X): {
      QubitPauliTensor t0 = zrow * QubitPauliTensor(qb, Pauli::Z, i_);
      QubitPauliTensor t1 = xrow * QubitPauliTensor(qb, Pauli::Y, i_);
      QubitPauliTensor t2 = QubitPauliTensor(qb, Pauli::Z);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1, t2});
      t0.coeff *= phases.first;
      t2.coeff *= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1), std::make_shared<const PGCliffordRot>(t2, 1)};
    }
    case switch_pauli_pair(Pauli::X, Pauli::Y):
    case switch_pauli_pair(Pauli::Y, Pauli::Y): {
      QubitPauliTensor t0 = zrow * QubitPauliTensor(qb, Pauli::Z, i_);
      QubitPauliTensor t1 = xrow * QubitPauliTensor(qb, Pauli::X, i_);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1});
      t0.coeff *= phases.first;
      t1.coeff *= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1)};
    }
    case switch_pauli_pair(Pauli::X, Pauli::Z): {
      QubitPauliTensor t0 = xrow * QubitPauliTensor(qb, Pauli::Y, i_);
      QubitPauliTensor t1 = zrow * QubitPauliTensor(qb, Pauli::Z, i_);
      QubitPauliTensor t2 = QubitPauliTensor(qb, Pauli::Z);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1, t2});
      t1.coeff *= phases.first;
      t2.coeff *= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1), std::make_shared<const PGCliffordRot>(t2, 1)};
    }
    case switch_pauli_pair(Pauli::Z, Pauli::I): {
      QubitPauliTensor t0 = zrow * QubitPauliTensor(qb, Pauli::X, i_);
      QubitPauliTensor t1 = zrow * xrow * QubitPauliTensor(qb, Pauli::Y);
      QubitPauliTensor t2 = QubitPauliTensor(qb, Pauli::Y);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1, t2});
      t1.coeff *= phases.first * phases.second;
      t2.coeff *= phases.second;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1), std::make_shared<const PGCliffordRot>(t2, 1)};
    }
    case switch_pauli_pair(Pauli::Z, Pauli::X): {
      std::vector<PGOp_ptr> pgops;
      bool zphase = false;
      if (zrow.string != QubitPauliString(qb, Pauli::Z)) {
        QubitPauliTensor t0 = zrow * QubitPauliTensor(qb, Pauli::Y, i_);
        QubitPauliTensor t1 = QubitPauliTensor(qb, Pauli::X);
        std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, QubitPauliTensor(qb, Pauli::X), qb, {t0, t1});
        t1.coeff *= phases.first;
        pgops.push_back(std::make_shared<const PGCliffordRot>(t0, 1));
        pgops.push_back(std::make_shared<const PGCliffordRot>(t1, 1));
      }
      else if (zrow.coeff == -1.) {
        zphase = true;
      }
      if (xrow.string != QubitPauliString(qb, Pauli::X)) {
        QubitPauliTensor t0 = xrow * QubitPauliTensor(qb, Pauli::Y, i_);
        QubitPauliTensor t1 = QubitPauliTensor(qb, Pauli::Z);
        std::pair<Complex, Complex> phases = verify_solution_to_phase(QubitPauliTensor(qb, Pauli::Z), xrow, qb, {t0, t1});
        t1.coeff *= phases.second;
        pgops.push_back(std::make_shared<const PGCliffordRot>(t0, 1));
        pgops.push_back(std::make_shared<const PGCliffordRot>(t1, 1));
      }
      else if (xrow.coeff == -1.) {
        if (zphase) pgops.push_back(std::make_shared<const PGCliffordRot>(QubitPauliTensor(qb, Y), 2));
        else pgops.push_back(std::make_shared<const PGCliffordRot>(QubitPauliTensor(qb, Z), 2));
      }
      else if (zphase) {
        pgops.push_back(std::make_shared<const PGCliffordRot>(QubitPauliTensor(qb, X), 2));
      }
      return pgops;
    }
    case switch_pauli_pair(Pauli::Z, Pauli::Y):
    case switch_pauli_pair(Pauli::Z, Pauli::Z): {
      QubitPauliTensor t0 = xrow * QubitPauliTensor(qb, Pauli::X, i_);
      QubitPauliTensor t1 = zrow * QubitPauliTensor(qb, Pauli::Y, i_);
      QubitPauliTensor t2 = QubitPauliTensor(qb, Pauli::X);
      std::pair<Complex, Complex> phases = verify_solution_to_phase(zrow, xrow, qb, {t0, t1, t2});
      t0.coeff *= phases.second;
      t2.coeff *= phases.first;
      return {std::make_shared<const PGCliffordRot>(t0, 1), std::make_shared<const PGCliffordRot>(t1, 1), std::make_shared<const PGCliffordRot>(t2, 1)};
    }
    default: throw std::logic_error("Unexpected Paulis when solving tableau with gadgets");
  }
}

PauliGraph circuit_to_pauli_graph(const Circuit& circ) {
  PauliGraph res(circ.all_qubits(), circ.all_bits());
  for (const Command& com : circ) {
    unit_vector_t args = com.get_args();
    for (const PGOp_ptr& pgop : op_to_pgops(com.get_op_ptr(), args, res, true)) {
      res.add_vertex_at_end(pgop);
    }
  }
  PauliStabiliserList initial_strings;
  tableau_row_index_t initial_rows;
  unsigned initial_i = 0;
  PauliStabiliserList final_strings;
  tableau_row_index_t final_rows;
  unsigned final_i = 0;
  qubit_map_t perm = circ.implicit_qubit_permutation();
  for (const Qubit& qb : circ.all_qubits()) {
    initial_strings.push_back(res.initial_.get_pauli(res.initial_rows_.left.at(std::pair<Qubit, TableauRowType>{qb, TableauRowType::ZRow})));
    initial_rows.insert({{qb, TableauRowType::ZRow}, initial_i});
    ++initial_i;
    if (!circ.is_created(qb)) {
      initial_strings.push_back(res.initial_.get_pauli(res.initial_rows_.left.at(std::pair<Qubit, TableauRowType>{qb, TableauRowType::XRow})));
      initial_rows.insert({{qb, TableauRowType::XRow}, initial_i});
      ++initial_i;
    }
    if (!circ.is_discarded(qb)) {
      final_strings.push_back(res.final_.get_pauli(res.final_rows_.left.at(std::pair<Qubit, TableauRowType>{qb, TableauRowType::ZRow})));
      final_rows.insert({{perm.at(qb), TableauRowType::ZRow}, final_i});
      ++final_i;
      final_strings.push_back(res.final_.get_pauli(res.final_rows_.left.at(std::pair<Qubit, TableauRowType>{qb, TableauRowType::XRow})));
      final_rows.insert({{perm.at(qb), TableauRowType::XRow}, final_i});
      ++final_i;
    }
  }
  res.initial_ = SymplecticTableau(initial_strings);
  res.initial_rows_ = initial_rows;
  res.final_ = SymplecticTableau(final_strings);
  res.final_rows_ = final_rows;
  return res;
}

// void synthesise_pgop(const PGOp_ptr& op, Circuit& circ, CXConfigType cx_config) {
//   switch (op->get_type()) {
//     case PGOpType::Rotation:
//     case PGOpType::CliffordRot:
//     case PGOpType::Measure:
//     case PGOpType::Decoherence:
//     case PGOpType::Reset:
//     case PGOpType::Conditional:
//     case PGOpType::Box:
//     default:
//   }
// }

// Circuit pauli_graph_to_circuit_individual(const PauliGraph& pg, CXConfigType cx_config) {
//   Circuit circ;
//   for (tableau_col_index_t::left_const_iterator it = pg.initial_cols_.left.begin(); it != pg.initial_cols_.left.end(); ++it) {
//     circ.add_qubit(it->first);
//   }
//   for (const Bit &b : pg.bits_) circ.add_bit(b);
//   for (PauliGraph::TopSortIterator it = pg.begin(); it != pg.end(); ++it) {
//     PGOp_ptr pgop = pg.graph_[*it].op_;
//     synthesise_pgop(pgop, circ, cx_config);
//   }
// }

}  // namespace tket
