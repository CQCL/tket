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

#include "tket/PauliGraphRefactor/PGBox.hpp"

#include <tkassert/Assert.hpp>

namespace tket {

namespace pg {

/**
 * PGBox Implementation
 */

Op_ptr PGBox::get_op() const { return op_; }

const unit_vector_t& PGBox::get_args() const { return args_; }

PGBox::PGBox(
    const Op_ptr& op, const unit_vector_t& args,
    const std::vector<SpPauliStabiliser>& paulis)
    : PGOp(PGOpType::Box), op_(op), args_(args), paulis_(paulis) {
  op_signature_t sig = op_->get_signature();
  unsigned nqs = 0;
  for (const EdgeType& et : sig) {
    if (et == EdgeType::Quantum) ++nqs;
  }
  if (paulis.size() != 2 * nqs)
    throw PGError(
        "Cannot create PGBox; number of SpPauliStabilisers must match twice "
        "the "
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

PGOp_ptr PGBox::clone() const {
  return std::make_shared<PGBox>(op_, args_, paulis_);
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

std::vector<SpPauliStabiliser> PGBox::active_paulis() const { return paulis_; }

PGOp_signature PGBox::pauli_signature() const {
  PGOp_signature sig;
  for (auto it = paulis_.begin(); it != paulis_.end();) {
    SpPauliStabiliser zp = *it;
    ++it;
    TKET_ASSERT(it != paulis_.end());
    sig.ac_pairs.push_back({zp, *it});
    ++it;
  }
  return sig;
}

SpPauliStabiliser& PGBox::port(unsigned p) {
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

/**
 * PGMultiplexedTensoredBox Implementation
 */

const std::map<std::vector<bool>, std::vector<Op_ptr>>&
PGMultiplexedTensoredBox::get_op_map() const {
  return op_map_;
}

const std::vector<SpPauliStabiliser>&
PGMultiplexedTensoredBox::get_control_paulis() const {
  return control_paulis_;
}

const std::vector<SpPauliStabiliser>&
PGMultiplexedTensoredBox::get_target_paulis() const {
  return target_paulis_;
}

PGMultiplexedTensoredBox::PGMultiplexedTensoredBox(
    const std::map<std::vector<bool>, std::vector<Op_ptr>>& op_map,
    const std::vector<SpPauliStabiliser>& control_paulis,
    const std::vector<SpPauliStabiliser>& target_paulis)
    : PGOp(PGOpType::MultiplexedTensoredBox),
      op_map_(op_map),
      control_paulis_(control_paulis),
      target_paulis_(target_paulis) {
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    if (it->first.size() != control_paulis_.size())
      throw PGError(
          "PGMultiplexedTensoredBox: Size mismatch between number of controls "
          "and length of values in op_map");
    unsigned total_qubits = 0;
    for (const Op_ptr& op : it->second) total_qubits += op->n_qubits();
    if (total_qubits * 2 != target_paulis_.size())
      throw PGError(
          "PGMultiplexedTensoredBox: Size mismatch between the number of "
          "qubits in tensored op and number of target qubits expected.");
  }
}

SymSet PGMultiplexedTensoredBox::free_symbols() const {
  SymSet sset;
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    for (const Op_ptr& op : it->second) {
      SymSet it_sset = op->free_symbols();
      sset.insert(it_sset.begin(), it_sset.end());
    }
  }
  return sset;
}

PGOp_ptr PGMultiplexedTensoredBox::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  std::map<std::vector<bool>, std::vector<Op_ptr>> new_op_map = op_map_;
  bool any_changed = false;
  for (auto it = new_op_map.begin(); it != new_op_map.end(); ++it) {
    for (auto op_it = it->second.begin(); op_it != it->second.end(); ++op_it) {
      Op_ptr new_op = (*op_it)->symbol_substitution(sub_map);
      if (new_op) {
        *op_it = new_op;
        any_changed = true;
      }
    }
  }
  if (any_changed)
    return std::make_shared<PGMultiplexedTensoredBox>(
        new_op_map, control_paulis_, target_paulis_);
  else
    return PGOp_ptr();
}

PGOp_ptr PGMultiplexedTensoredBox::clone() const {
  return std::make_shared<PGMultiplexedTensoredBox>(
      op_map_, control_paulis_, target_paulis_);
}

std::string PGMultiplexedTensoredBox::get_name(bool latex) const {
  std::stringstream str;
  str << "qswitch [";
  if (!control_paulis_.empty()) {
    str << control_paulis_.at(0).to_str();
    for (unsigned i = 1; i < control_paulis_.size(); ++i) {
      str << ", " << control_paulis_.at(i).to_str();
    }
  }
  str << "]";
  bool first = true;
  for (auto it = op_map_.begin(); it != op_map_.end(); ++it) {
    if (!first) {
      str << ", ";
    }
    for (bool b : it->first) {
      str << (b ? "1" : "0");
    }
    str << "->[";
    bool first_op = true;
    for (auto op_it = it->second.begin(); op_it != it->second.end(); ++op_it) {
      if (!first_op) str << " x ";
      str << (*op_it)->get_name(latex);
      first_op = false;
    }
    str << "]";
    first = false;
  }
  return str.str();
}

bool PGMultiplexedTensoredBox::is_equal(const PGOp& op_other) const {
  const PGMultiplexedTensoredBox& other =
      dynamic_cast<const PGMultiplexedTensoredBox&>(op_other);
  return (control_paulis_ == other.control_paulis_) &&
         (target_paulis_ == other.target_paulis_) && (op_map_ == other.op_map_);
}

unsigned PGMultiplexedTensoredBox::n_paulis() const {
  return control_paulis_.size() + target_paulis_.size();
}

std::vector<SpPauliStabiliser> PGMultiplexedTensoredBox::active_paulis() const {
  std::vector<SpPauliStabiliser> aps = control_paulis_;
  aps.insert(aps.end(), target_paulis_.begin(), target_paulis_.end());
  return aps;
}

PGOp_signature PGMultiplexedTensoredBox::pauli_signature() const {
  PGOp_signature sig;
  for (auto it = control_paulis_.begin(); it != control_paulis_.end(); ++it) {
    sig.c.push_back(*it);
  }
  for (auto it = target_paulis_.begin(); it != target_paulis_.end();) {
    SpPauliStabiliser zp = *it;
    ++it;
    TKET_ASSERT(it != target_paulis_.end());
    sig.ac_pairs.push_back({zp, *it});
    ++it;
  }
  return sig;
}

SpPauliStabiliser& PGMultiplexedTensoredBox::port(unsigned p) {
  if (p < control_paulis_.size()) return control_paulis_.at(p);
  if (p < control_paulis_.size() + target_paulis_.size())
    return target_paulis_.at(p - control_paulis_.size());
  else
    throw PGError(
        "Cannot dereference port of PGMultiplexedTensoredBox: " +
        std::to_string(p));
}

}  // namespace pg
}  // namespace tket
