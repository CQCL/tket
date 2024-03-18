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

#pragma once

#include "tket/Ops/Op.hpp"
#include "tket/PauliGraphRefactor/PGOp.hpp"

namespace tket {

namespace pg {

/**
 * PGOp for PGOpType::Box, representing an arbitrary Op conjugated by some
 * Clifford circuit. For each qubit the Op acts on, we maintain two
 * active_paulis corresponding to the Pauli operators mapped into +Z and +X by
 * the conjugating circuit. This allows it to be resynthesised as a unitary
 * extension circuit for the corresponding ChoiMixTableau.
 */
class PGBox : public PGOp {
 public:
  /**
   * Get the Op captured within the box.
   */
  Op_ptr get_op() const;

  /**
   * Get the original arguments (both Qubits and Bits) of the Op as used in the
   * original circuit. Any Qubits in this list are treated as placeholders and
   * will be replaced with new Qubits at the point of synthesis depending on
   * which are easiest for synthesising the conjugating Clifford circuit. We
   * assume any Bits could be both read and written to.
   */
  const unit_vector_t& get_args() const;

  /**
   * Constructs a black box abstraction of the given Op \p op applied to \p
   * args. Any Qubits in \p args are treated as placeholders for the Qubits used
   * at the point of synthesis, and we assume any Bits could be both read and
   * written to by \p op.
   *
   * \p paulis specifies the Clifford conjugation around \p op. Specifically,
   * for each i in [0..op->n_qubits()-1]:
   *   - paulis_[2*i] is the Pauli operator mapped into Z_q[i]
   *   - paulis_[2*i+1] is the Pauli operator mapped into X_q[i]
   */
  PGBox(
      const Op_ptr& op, const unit_vector_t& args,
      const std::vector<SpPauliStabiliser>& paulis);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;
  virtual bit_vector_t read_bits() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  Op_ptr op_;
  unit_vector_t args_;
  /**
   * Semantics of paulis_:
   * For each i in [0..op->n_qubits()-1]:
   *   - paulis_[2*i] is the Pauli string Z_q[i] is mapped to by Clifford conj
   *   - paulis_[2*i+1] is the Pauli string X_q[i] is mapped to by Clifford conj
   */
  std::vector<SpPauliStabiliser> paulis_;
};

}  // namespace pg

}  // namespace tket
