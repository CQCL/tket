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

#include "tket/PauliGraphRefactor/PGOp.hpp"

namespace tket {
namespace pg {

/**
 * PGOp Abstract Implementation
 */

PGOpType PGOp::get_type() const { return type_; }

bool PGOp::operator==(const PGOp& other) const {
  return (type_ == other.type_) && is_equal(other);
}

PGOp::~PGOp() {}

PGOp::PGOp(PGOpType type) : type_(type) {}

bool PGOp::commutes_with(const PGOp& other) const {
  for (const SpPauliStabiliser& t : active_paulis()) {
    for (const SpPauliStabiliser& ot : other.active_paulis()) {
      if (!t.commutes_with(ot)) return false;
    }
  }
  for (const Bit& b : write_bits()) {
    for (const Bit& ob : other.write_bits()) {
      if (b == ob) return false;
    }
    for (const Bit& ob : other.read_bits()) {
      if (b == ob) return false;
    }
  }
  for (const Bit& b : read_bits()) {
    for (const Bit& ob : other.write_bits()) {
      if (b == ob) return false;
    }
  }
  return true;
}

unsigned PGOp::n_paulis() const { return 1; }

bit_vector_t PGOp::read_bits() const { return {}; }

bit_vector_t PGOp::write_bits() const { return {}; }

/**
 * PGRotation Implementation
 */

const SpPauliStabiliser& PGRotation::get_tensor() const { return tensor_; }

const Expr& PGRotation::get_angle() const { return angle_; }

PGRotation::PGRotation(const SpPauliStabiliser& tensor, const Expr& angle)
    : PGOp(PGOpType::Rotation), tensor_(tensor), angle_(angle) {
  if (tensor.is_real_negative()) {
    angle_ *= -1;
    tensor_.coeff = 0;
  }
}

SymSet PGRotation::free_symbols() const { return expr_free_symbols(angle_); }

PGOp_ptr PGRotation::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  return std::make_shared<PGRotation>(tensor_, angle_.subs(sub_map));
}

std::string PGRotation::get_name(bool) const {
  std::stringstream str;
  str << "Rot(" << tensor_.to_str() << "; " << angle_ << ")";
  return str.str();
}

bool PGRotation::is_equal(const PGOp& op_other) const {
  const PGRotation& other = dynamic_cast<const PGRotation&>(op_other);
  return (tensor_ == other.tensor_) && equiv_expr(angle_, other.angle_, 2.);
}

std::vector<SpPauliStabiliser> PGRotation::active_paulis() const {
  return {tensor_};
}

SpPauliStabiliser& PGRotation::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGRotation: " + std::to_string(p));
  return tensor_;
}

/**
 * PGCliffordRot Implementation
 */

const SpPauliStabiliser& PGCliffordRot::get_tensor() const { return tensor_; }

unsigned PGCliffordRot::get_angle() const { return angle_; }

PGCliffordRot::PGCliffordRot(const SpPauliStabiliser& tensor, unsigned angle)
    : PGOp(PGOpType::CliffordRot), tensor_(tensor), angle_(angle) {
  if (tensor.is_real_negative()) {
    angle_ = (4 - angle) % 4;
    tensor_.coeff = 0;
  }
}

SymSet PGCliffordRot::free_symbols() const { return {}; }

PGOp_ptr PGCliffordRot::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGCliffordRot::get_name(bool) const {
  std::stringstream str;
  str << "ClfRot(" << tensor_.to_str() << "; " << (angle_ * 0.5) << ")";
  return str.str();
}

bool PGCliffordRot::is_equal(const PGOp& op_other) const {
  const PGCliffordRot& other = dynamic_cast<const PGCliffordRot&>(op_other);
  return (tensor_ == other.tensor_) && (angle_ % 4 == other.angle_ % 4);
}

std::vector<SpPauliStabiliser> PGCliffordRot::active_paulis() const {
  return {tensor_};
}

SpPauliStabiliser& PGCliffordRot::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGCliffordRot: " + std::to_string(p));
  return tensor_;
}

/**
 * PGMeasure Implementation
 */

const SpPauliStabiliser& PGMeasure::get_tensor() const { return tensor_; }

const Bit& PGMeasure::get_target() const { return target_; }

PGMeasure::PGMeasure(const SpPauliStabiliser& tensor, const Bit& target)
    : PGOp(PGOpType::Measure), tensor_(tensor), target_(target) {
  // Assert that tensor has a real coefficient
  tensor.is_real_negative();
}

SymSet PGMeasure::free_symbols() const { return {}; }

PGOp_ptr PGMeasure::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGMeasure::get_name(bool) const {
  std::stringstream str;
  str << "Meas(" << tensor_.to_str() << " -> " << target_.repr() << ")";
  return str.str();
}

bool PGMeasure::is_equal(const PGOp& op_other) const {
  const PGMeasure& other = dynamic_cast<const PGMeasure&>(op_other);
  return (tensor_ == other.tensor_) && (target_ == other.target_);
}

std::vector<SpPauliStabiliser> PGMeasure::active_paulis() const {
  return {tensor_};
}

SpPauliStabiliser& PGMeasure::port(unsigned p) {
  if (p != 0)
    throw PGError("Cannot dereference port on PGMeasure: " + std::to_string(p));
  return tensor_;
}

bit_vector_t PGMeasure::write_bits() const { return {target_}; }

/**
 * PGDecoherence Implementation
 */

const SpPauliStabiliser& PGDecoherence::get_tensor() const { return tensor_; }

PGDecoherence::PGDecoherence(const SpPauliStabiliser& tensor)
    : PGOp(PGOpType::Decoherence), tensor_(tensor) {
  if (tensor.is_real_negative()) {
    tensor_.coeff = 0;
  }
}

SymSet PGDecoherence::free_symbols() const { return {}; }

PGOp_ptr PGDecoherence::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGDecoherence::get_name(bool) const {
  std::stringstream str;
  str << "Deco(" << tensor_.to_str() << ")";
  return str.str();
}

bool PGDecoherence::is_equal(const PGOp& op_other) const {
  const PGDecoherence& other = dynamic_cast<const PGDecoherence&>(op_other);
  return tensor_ == other.tensor_;
}

std::vector<SpPauliStabiliser> PGDecoherence::active_paulis() const {
  return {tensor_};
}

SpPauliStabiliser& PGDecoherence::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGDecoherence: " + std::to_string(p));
  return tensor_;
}

/**
 * PGReset Implementation
 */

const SpPauliStabiliser& PGReset::get_stab() const { return stab_; }

const SpPauliStabiliser& PGReset::get_destab() const { return destab_; }

PGReset::PGReset(const SpPauliStabiliser& stab, const SpPauliStabiliser& destab)
    : PGOp(PGOpType::Reset), stab_(stab), destab_(destab) {
  if (destab.is_real_negative()) destab_.coeff = 0;
  // Assert that stab has a real coefficient
  stab.is_real_negative();
}

SymSet PGReset::free_symbols() const { return {}; }

PGOp_ptr PGReset::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGReset::get_name(bool) const {
  std::stringstream str;
  str << "Reset(" << stab_.to_str() << "; " << destab_.to_str() << ")";
  return str.str();
}

bool PGReset::is_equal(const PGOp& op_other) const {
  const PGReset& other = dynamic_cast<const PGReset&>(op_other);
  return (stab_ == other.stab_) && (destab_ == other.destab_);
}

unsigned PGReset::n_paulis() const { return 2; }

std::vector<SpPauliStabiliser> PGReset::active_paulis() const {
  return {stab_, destab_};
}

SpPauliStabiliser& PGReset::port(unsigned p) {
  if (p == 0)
    return stab_;
  else if (p == 1)
    return destab_;
  else
    throw PGError("Cannot dereference port on PGReset: " + std::to_string(p));
}

/**
 * PGConditional Implementation
 */

PGOp_ptr PGConditional::get_inner_op() const { return inner_; }

const bit_vector_t& PGConditional::get_args() const { return args_; }

unsigned PGConditional::get_value() const { return value_; }

PGConditional::PGConditional(
    PGOp_ptr inner, const bit_vector_t& args, unsigned value)
    : PGOp(PGOpType::Conditional), inner_(inner), args_(args), value_(value) {}

SymSet PGConditional::free_symbols() const { return inner_->free_symbols(); }

PGOp_ptr PGConditional::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  PGOp_ptr inner_sub = inner_->symbol_substitution(sub_map);
  if (inner_sub)
    return std::make_shared<PGConditional>(inner_sub, args_, value_);
  else
    return PGOp_ptr();
}

std::string PGConditional::get_name(bool latex) const {
  std::stringstream str;
  str << "[";
  bool first = true;
  for (const Bit& b : args_) {
    if (first)
      first = false;
    else
      str << ", ";
    str << b.repr();
  }
  str << "] == " << value_ << " ? " << inner_->get_name(latex);
  return str.str();
}

bool PGConditional::is_equal(const PGOp& op_other) const {
  const PGConditional& other = dynamic_cast<const PGConditional&>(op_other);
  return (value_ == other.value_) && (args_ == other.args_) &&
         (*inner_ == *other.inner_);
}

unsigned PGConditional::n_paulis() const { return inner_->n_paulis(); }

std::vector<SpPauliStabiliser> PGConditional::active_paulis() const {
  return inner_->active_paulis();
}

SpPauliStabiliser& PGConditional::port(unsigned p) { return inner_->port(p); }

bit_vector_t PGConditional::read_bits() const {
  bit_vector_t bits = inner_->read_bits();
  bits.insert(bits.end(), args_.begin(), args_.end());
  return bits;
}

bit_vector_t PGConditional::write_bits() const { return inner_->write_bits(); }

/**
 * PGQControl Implementation
 */

PGOp_ptr PGQControl::get_inner_op() const { return inner_; }

const std::vector<SpPauliStabiliser>& PGQControl::get_control_paulis() const {
  return control_paulis_;
}

std::vector<bool> PGQControl::get_value() const { return value_; }

PGQControl::PGQControl(
    PGOp_ptr inner, const std::vector<SpPauliStabiliser>& control_paulis,
    std::vector<bool> value)
    : PGOp(PGOpType::QControl),
      inner_(inner),
      control_paulis_(control_paulis),
      value_(value) {
  if (control_paulis_.size() != value_.size())
    throw PGError(
        "PGQControl: Size mismatch between number of controls and length of "
        "value");
}

SymSet PGQControl::free_symbols() const { return inner_->free_symbols(); }

PGOp_ptr PGQControl::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  PGOp_ptr inner_sub = inner_->symbol_substitution(sub_map);
  if (inner_sub)
    return std::make_shared<PGQControl>(inner_sub, control_paulis_, value_);
  else
    return PGOp_ptr();
}

std::string PGQControl::get_name(bool latex) const {
  std::stringstream str;
  str << "qif (";
  if (!control_paulis_.empty()) {
    str << (value_.at(0) ? "-" : "") << control_paulis_.at(0).to_str();
    for (unsigned i = 1; i < control_paulis_.size(); ++i) {
      str << ", " << (value_.at(i) ? "-" : "")
          << control_paulis_.at(i).to_str();
    }
  }
  str << ") " << inner_->get_name(latex);
  return str.str();
}

bool PGQControl::is_equal(const PGOp& op_other) const {
  const PGQControl& other = dynamic_cast<const PGQControl&>(op_other);
  return (value_ == other.value_) &&
         (control_paulis_ == other.control_paulis_) &&
         (*inner_ == *other.inner_);
}

unsigned PGQControl::n_paulis() const {
  return control_paulis_.size() + inner_->n_paulis();
}

std::vector<SpPauliStabiliser> PGQControl::active_paulis() const {
  std::vector<SpPauliStabiliser> aps = control_paulis_;
  std::vector<SpPauliStabiliser> inner_aps = inner_->active_paulis();
  aps.insert(aps.end(), inner_aps.begin(), inner_aps.end());
  return aps;
}

SpPauliStabiliser& PGQControl::port(unsigned p) {
  if (p < control_paulis_.size())
    return control_paulis_.at(p);
  else
    return inner_->port(p - control_paulis_.size());
}

/**
 * PGMultiplexor Implementation
 */

const std::map<std::vector<bool>, PGOp_ptr>& PGMultiplexor::get_inner_op_map()
    const {
  return op_map_;
}

const std::vector<SpPauliStabiliser>& PGMultiplexor::get_control_paulis()
    const {
  return control_paulis_;
}

PGMultiplexor::PGMultiplexor(
    const std::map<std::vector<bool>, PGOp_ptr>& op_map,
    const std::vector<SpPauliStabiliser>& control_paulis)
    : PGOp(PGOpType::Multiplexor),
      op_map_(op_map),
      control_paulis_(control_paulis) {
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    if (it->first.size() != control_paulis_.size())
      throw PGError(
          "PGMultiplexor: Size mismatch between number of controls and length "
          "of values in op_map");
  }
}

SymSet PGMultiplexor::free_symbols() const {
  SymSet sset;
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    SymSet it_sset = it->second->free_symbols();
    sset.insert(it_sset.begin(), it_sset.end());
  }
  return sset;
}

PGOp_ptr PGMultiplexor::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  std::map<std::vector<bool>, PGOp_ptr> new_op_map = op_map_;
  bool any_change = false;
  for (auto it = new_op_map.begin(); it != new_op_map.end(); ++it) {
    PGOp_ptr new_op = it->second->symbol_substitution(sub_map);
    if (new_op) {
      any_change = true;
      it->second = new_op;
    }
  }
  if (any_change)
    return std::make_shared<PGMultiplexor>(new_op_map, control_paulis_);
  else
    return PGOp_ptr();
}

std::string PGMultiplexor::get_name(bool latex) const {
  std::stringstream str;
  str << "qswitch [";
  if (!control_paulis_.empty()) {
    str << control_paulis_.at(0).to_str();
    for (unsigned i = 1; i < control_paulis_.size(); ++i) {
      str << ", " << control_paulis_.at(i).to_str();
    }
  }
  str << "]";
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    str << ", ";
    for (bool b : it->first) {
      str << (b ? "1" : "0");
    }
    str << "->[" << it->second->get_name(latex) << "]";
  }
  return str.str();
}

bool PGMultiplexor::is_equal(const PGOp& op_other) const {
  const PGMultiplexor& other = dynamic_cast<const PGMultiplexor&>(op_other);
  if (control_paulis_ != other.control_paulis_) return false;
  auto this_it = op_map_.begin();
  auto other_it = other.op_map_.end();
  while (this_it != op_map_.end()) {
    if (other_it == other.op_map_.end()) return false;
    if (this_it->first != other_it->first) return false;
    if (*(this_it->second) != *(other_it->second)) return false;
    ++this_it;
    ++other_it;
  }
  return (other_it == other.op_map_.end());
}

unsigned PGMultiplexor::n_paulis() const {
  unsigned n_paulis = control_paulis_.size();
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    n_paulis += it->second->n_paulis();
  }
  return n_paulis;
}

std::vector<SpPauliStabiliser> PGMultiplexor::active_paulis() const {
  std::vector<SpPauliStabiliser> aps = control_paulis_;
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    std::vector<SpPauliStabiliser> it_aps = it->second->active_paulis();
    aps.insert(aps.end(), it_aps.begin(), it_aps.end());
  }
  return aps;
}

SpPauliStabiliser& PGMultiplexor::port(unsigned p) {
  unsigned original_p = p;
  if (p < control_paulis_.size()) return control_paulis_.at(p);
  p -= control_paulis_.size();
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    unsigned it_paulis = it->second->n_paulis();
    if (p < it_paulis) return it->second->port(p);
    p -= it_paulis;
  }
  throw PGError(
      "Cannot dereference port of PGMultiplexor: " +
      std::to_string(original_p));
}

/**
 * PGStabAssertion Implementation
 */

const SpPauliStabiliser& PGStabAssertion::get_stab() const { return stab_; }

const SpPauliStabiliser& PGStabAssertion::get_anc_z() const { return anc_z_; }

const SpPauliStabiliser& PGStabAssertion::get_anc_x() const { return anc_x_; }

const Bit& PGStabAssertion::get_target() const { return target_; }

PGStabAssertion::PGStabAssertion(
    const SpPauliStabiliser& stab, const SpPauliStabiliser& anc_z,
    const SpPauliStabiliser& anc_x, const Bit& target)
    : PGOp(PGOpType::StabAssertion),
      stab_(stab),
      anc_z_(anc_z),
      anc_x_(anc_x),
      target_(target) {
  if (anc_x.is_real_negative()) anc_x_.coeff = 0;
  // Assert that stab and anc_z have real coefficients
  stab.is_real_negative();
  anc_z.is_real_negative();
}

SymSet PGStabAssertion::free_symbols() const { return {}; }

PGOp_ptr PGStabAssertion::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGStabAssertion::get_name(bool) const {
  std::stringstream str;
  str << "Stab(" << stab_.to_str() << " -> " << target_.repr() << "; "
      << anc_z_.to_str() << ", " << anc_x_.to_str() << ")";
  return str.str();
}

bool PGStabAssertion::is_equal(const PGOp& op_other) const {
  const PGStabAssertion& other = dynamic_cast<const PGStabAssertion&>(op_other);
  return (stab_ == other.stab_) && (anc_z_ == other.anc_z_) &&
         (anc_x_ == other.anc_x_) && (target_ == other.target_);
}

unsigned PGStabAssertion::n_paulis() const { return 3; }

std::vector<SpPauliStabiliser> PGStabAssertion::active_paulis() const {
  return {stab_, anc_z_, anc_x_};
}

SpPauliStabiliser& PGStabAssertion::port(unsigned p) {
  switch (p) {
    case 0:
      return stab_;
    case 1:
      return anc_z_;
    case 2:
      return anc_x_;
    default:
      throw PGError(
          "Cannot dereference port of PGStabAssertion: " + std::to_string(p));
  }
}

/**
 * PGInputTableau Implementation
 */

const ChoiMixTableau::row_tensor_t& PGInputTableau::get_full_row(
    unsigned p) const {
  if (p >= rows_.size())
    throw PGError(
        "Cannot dereference row on PGInputTableau: " + std::to_string(p));
  return rows_.at(p);
}

ChoiMixTableau PGInputTableau::to_cm_tableau() const {
  std::list<ChoiMixTableau::row_tensor_t> rows;
  for (unsigned i = 0; i < n_paulis(); ++i) {
    rows.push_back(get_full_row(i));
  }
  return ChoiMixTableau(rows);
}

PGInputTableau::PGInputTableau(const ChoiMixTableau& tableau)
    : PGOp(PGOpType::InputTableau), rows_() {
  for (unsigned i = 0; i < tableau.get_n_rows(); ++i)
    rows_.push_back(tableau.get_row(i));
}

SymSet PGInputTableau::free_symbols() const { return {}; }

PGOp_ptr PGInputTableau::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGInputTableau::get_name(bool) const {
  std::stringstream str;
  str << "Input(\n";
  for (const ChoiMixTableau::row_tensor_t& row : rows_)
    str << "\t" << row.first.to_str() << "\t->\t" << row.second.to_str()
        << "\n";
  str << "\n)";
  return str.str();
}

bool PGInputTableau::is_equal(const PGOp& op_other) const {
  const PGInputTableau& other = dynamic_cast<const PGInputTableau&>(op_other);
  return (rows_ == other.rows_);
}

unsigned PGInputTableau::n_paulis() const { return rows_.size(); }

std::vector<SpPauliStabiliser> PGInputTableau::active_paulis() const {
  std::vector<SpPauliStabiliser> paulis;
  for (const ChoiMixTableau::row_tensor_t& row : rows_)
    paulis.push_back(row.second);
  return paulis;
}

SpPauliStabiliser& PGInputTableau::port(unsigned p) {
  if (p >= rows_.size())
    throw PGError(
        "Cannot dereference port on PGInputTableau: " + std::to_string(p));
  return rows_.at(p).second;
}

/**
 * PGOutputTableau Implementation
 */

const ChoiMixTableau::row_tensor_t& PGOutputTableau::get_full_row(
    unsigned p) const {
  if (p >= rows_.size())
    throw PGError(
        "Cannot dereference row on PGOutputTableau: " + std::to_string(p));
  return rows_.at(p);
}

ChoiMixTableau PGOutputTableau::to_cm_tableau() const {
  std::list<ChoiMixTableau::row_tensor_t> rows;
  for (unsigned i = 0; i < n_paulis(); ++i) {
    rows.push_back(get_full_row(i));
  }
  return ChoiMixTableau(rows);
}

PGOutputTableau::PGOutputTableau(const ChoiMixTableau& tableau)
    : PGOp(PGOpType::OutputTableau), rows_() {
  for (unsigned i = 0; i < tableau.get_n_rows(); ++i)
    rows_.push_back(tableau.get_row(i));
}

SymSet PGOutputTableau::free_symbols() const { return {}; }

PGOp_ptr PGOutputTableau::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGOutputTableau::get_name(bool) const {
  std::stringstream str;
  str << "Output(\n";
  for (const ChoiMixTableau::row_tensor_t& row : rows_)
    str << "\t" << row.first.to_str() << "\t->\t" << row.second.to_str()
        << "\n";
  str << "\n)";
  return str.str();
}

bool PGOutputTableau::is_equal(const PGOp& op_other) const {
  const PGOutputTableau& other = dynamic_cast<const PGOutputTableau&>(op_other);
  return (rows_ == other.rows_);
}

unsigned PGOutputTableau::n_paulis() const { return rows_.size(); }

std::vector<SpPauliStabiliser> PGOutputTableau::active_paulis() const {
  std::vector<SpPauliStabiliser> paulis;
  for (const ChoiMixTableau::row_tensor_t& row : rows_)
    paulis.push_back(row.first);
  return paulis;
}

SpPauliStabiliser& PGOutputTableau::port(unsigned p) {
  if (p >= rows_.size())
    throw PGError(
        "Cannot dereference port on PGOutputTableau: " + std::to_string(p));
  return rows_.at(p).first;
}

}  // namespace pg
}  // namespace tket
