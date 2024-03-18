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

namespace tket {

namespace pg {

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

}  // namespace pg
}  // namespace tket
