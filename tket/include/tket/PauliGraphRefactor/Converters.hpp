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

#pragma once

#include "tket/Circuit/Circuit.hpp"
#include "tket/PauliGraphRefactor/PauliGraph.hpp"

namespace tket {

// Copied from Converters/Converters.hpp to allow synthesis and testing of reset
// and boxes
ChoiMixTableau circuit_to_cm_tableau(const Circuit& circ);
std::pair<Circuit, qubit_map_t> cm_tableau_to_exact_circuit(
    const ChoiMixTableau& tab, CXConfigType cx_config = CXConfigType::Snake);
std::pair<Circuit, qubit_map_t> cm_tableau_to_unitary_extension_circuit(
    const ChoiMixTableau& tab, const std::vector<Qubit>& init_names = {},
    const std::vector<Qubit>& post_names = {},
    CXConfigType cx_config = CXConfigType::Snake);

/**
 * Converts the Circuit to a PauliGraph representing the same circuit by
 * iterating through the circuit, accumulating a UnitaryRevTableau of Clifford
 * operations and yielding PGOps for any other Op. Qubit initialisations are
 * recorded as part of the PGInputTableau, and discards are incorporated into
 * the accumulated tableau to give the PGOutputTableau.
 */
pg::PauliGraph circuit_to_pauli_graph3(
    const Circuit& circ, bool collect_cliffords = true);

/**
 * THIS IS A LEGACY METHOD DESIGNED TO REPLICATE THE SYNTHESIS METHODS FOR THE
 * EXISTING PAULIGRAPH
 *
 * Synthesises a circuit equivalent to the PauliGraph by adding each
 * pauli gadget to the circuit as a PauliExpBox individually
 * in the order given by TopSortIterator.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph3_to_pauli_exp_box_circuit_individually(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * THIS IS A LEGACY METHOD DESIGNED TO REPLICATE THE SYNTHESIS METHODS FOR THE
 * EXISTING PAULIGRAPH
 *
 * Synthesises a circuit equivalent to the PauliGraph by inserting pairs of
 * pauli gadgets as PauliExpPairBoxes into the circuit
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph3_to_pauli_exp_box_circuit_pairwise(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * THIS IS A LEGACY METHOD DESIGNED TO REPLICATE THE SYNTHESIS METHODS FOR THE
 * EXISTING PAULIGRAPH
 *
 * Synthesises a circuit equivalent to the PauliGraph by building
 * sets of mutually commuting pauli gadgets and
 * inserting them into the circuit as PauliExpCommutingSetBoxes
 * The tableau is then synthesised at the end.
 */

Circuit pauli_graph3_to_pauli_exp_box_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by synthesising each
 * vertex individually in the order given by TopSortIterator. The tableaux are
 * synthesised at each end using the default synthesis method for
 * ChoiMixTableau.
 */
Circuit pauli_graph3_to_circuit_individual(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

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
