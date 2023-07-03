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
#include "tket/PauliGraph3/PauliGraph.hpp"

namespace tket {

// Copied from Converters/Converters.hpp to allow synthesis and testing of reset
// and boxes
ChoiMixTableau circuit_to_cm_tableau(const Circuit& circ);
std::pair<Circuit, unit_map_t> cm_tableau_to_circuit(
    const ChoiMixTableau& circ, ChoiMixSynthType synth_type);

pg::PauliGraph circuit_to_pauli_graph3(const Circuit& circ);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building each
 * pauli gadget individually in the order given by TopSortIterator.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph3_to_circuit_individual(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building pairs of
 * pauli gadgets simultaneously using the method detailed in Cowtan et al.
 * Phase Gadget Synthesis for Shallow Circuits, Lemma 4.9.
 * The tableau is then synthesised at the end.
 */

Circuit pauli_graph3_to_circuit_pairwise(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building
 * sets of mutually commuting pauli gadgets and simultaneously
 * diagonalizing each gadget in a set.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph3_to_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

namespace pg {

class PGBox : public PGOp {
 public:
  Op_ptr get_op() const;
  const unit_vector_t& get_args() const;

  PGBox(
      const Op_ptr& op, const unit_vector_t& args,
      const std::vector<QubitPauliTensor>& paulis);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;
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
  std::vector<QubitPauliTensor> paulis_;
};

}  // namespace pg

}  // namespace tket
