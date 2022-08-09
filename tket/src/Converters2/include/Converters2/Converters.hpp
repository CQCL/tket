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

#include "Circuit/Circuit.hpp"
#include "PauliGraph2/PauliGraph.hpp"

namespace tket {

pg::PauliGraph circuit_to_pauli_graph2(const Circuit& circ);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building each
 * pauli gadget individually in the order given by TopSortIterator.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph2_to_circuit_individual(const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building pairs of
 * pauli gadgets simultaneously using the method detailed in Cowtan et al.
 * Phase Gadget Synthesis for Shallow Circuits, Lemma 4.9.
 * The tableau is then synthesised at the end.
 */

Circuit pauli_graph2_to_circuit_pairwise(
    const pg::PauliGraph &pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building
 * sets of mutually commuting pauli gadgets and simultaneously
 * diagonalizing each gadget in a set.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph2_to_circuit_sets(
    const pg::PauliGraph &pg, CXConfigType cx_config = CXConfigType::Snake);

namespace pg {

class PGBox : public PGOp {
 public:
  Op_ptr get_op() const;
  const unit_vector_t& get_args() const;
  // const bitt_vector_t& get_bits() const;

  PGBox(const Op_ptr& op, const unit_vector_t& args);
  // PGBox(const Op_ptr& op, const bit_vector_t& bits, const std::vector<QubitPauliTensor>& ztens, std::vector<QubitPauliTensor>& xtens);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual bit_vector_t read_bits() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  Op_ptr op_;
  unit_vector_t args_;
  // bit_vector_t bits_;
  // SymplecticTableau tab_;
  // tableau_col_index_t tab_cols_;
};

}  // namespace pg

}  // namespace tket
