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

#include "PGBox.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Clifford/UnitaryTableau.hpp"
#include "tket/PauliGraphRefactor/PauliGraph.hpp"

namespace tket {

// Copied from Converters/Converters.hpp and Converters/Gauss.hpp to allow
// synthesis and testing of reset and boxes
ChoiMixTableau circuit_to_cm_tableau(const Circuit& circ);
std::pair<Circuit, qubit_map_t> cm_tableau_to_exact_circuit(
    const ChoiMixTableau& tab, CXConfigType cx_config = CXConfigType::Snake);
std::pair<Circuit, qubit_map_t> cm_tableau_to_unitary_extension_circuit(
    const ChoiMixTableau& tab, const std::vector<Qubit>& init_names = {},
    const std::vector<Qubit>& post_names = {},
    CXConfigType cx_config = CXConfigType::Snake);
UnitaryTableau circuit_to_unitary_tableau(const Circuit& circ);
Circuit unitary_tableau_to_circuit(const UnitaryTableau& tab);
ChoiAPState cm_tableau_to_choi_apstate(const ChoiMixTableau& tab);
std::pair<Circuit, qubit_map_t> choi_apstate_to_exact_circuit(
    ChoiAPState ap, CXConfigType cx_config = CXConfigType::Snake);
ChoiAPState circuit_to_choi_apstate(const Circuit& circ);
ChoiMixTableau choi_apstate_to_cm_tableau(const ChoiAPState& ap);
class CXMaker {
 public:
  explicit CXMaker(unsigned qubits, bool reverse_cx_dirs = false)
      : _circ(qubits), _reverse_cx_dirs(reverse_cx_dirs) {}
  void row_add(unsigned r0, unsigned r1);
  Circuit _circ;
  bool _reverse_cx_dirs;
};

class DiagMatrix {
 public:
  DiagMatrix() {}
  explicit DiagMatrix(const MatrixXb& matrix) : _matrix(matrix) {}
  void row_add(unsigned r0, unsigned r1);
  void col_add(unsigned c0, unsigned c1);
  void gauss(CXMaker& cxmaker, unsigned blocksize = 6);
  friend std::ostream& operator<<(std::ostream& out, const DiagMatrix& diam);
  bool is_id() const;
  bool is_id_until_columns(unsigned limit) const;
  unsigned n_rows() const;
  unsigned n_cols() const;

  MatrixXb _matrix;
};

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

Circuit pauli_graph3_to_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Based on PauliExpCommutingSetBox, this Box represents a collection of PGOps
 * from a PauliGraph which are mutually commuting and can hence be performed
 * simultaneously.
 *
 * generate_circuit() synthesises this by mutual diagonalisation and GraySynth
 * (though PGOps with more than one active Pauli which commutes with all others
 * as solved ad-hoc, as this would otherwise require GraySynth targets to be
 * sets of Pauli strings which need to be simultaneously available)
 */
class PGOpCommutingSetBox : public Box {
 public:
  PGOpCommutingSetBox(
      const std::vector<pg::PGOp_ptr>& pgops,
      const boost::bimap<Qubit, unsigned>& qubit_indices,
      const boost::bimap<Bit, unsigned>& bit_indices,
      CXConfigType cx_config = CXConfigType::Tree);

  /**
   * Construct from the empty vector
   */
  PGOpCommutingSetBox();

  /**
   * Copy constructor
   */
  PGOpCommutingSetBox(const PGOpCommutingSetBox& other);

  ~PGOpCommutingSetBox() override {}

  bool is_clifford() const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two instances
   */
  bool is_equal(const Op& op_other) const override;

  /** Get the PGOps */
  const std::vector<pg::PGOp_ptr>& get_pgops() const { return pgops_; }

  /** Get the cx_config parameter (affects box decomposition) */
  CXConfigType get_cx_config() const { return cx_config_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;

  static Op_ptr from_json(const nlohmann::json& j);

  static nlohmann::json to_json(const Op_ptr& op);

 protected:
  void generate_circuit() const override;

 private:
  std::vector<pg::PGOp_ptr> pgops_;
  boost::bimap<Qubit, unsigned> qubit_indices_;
  boost::bimap<Bit, unsigned> bit_indices_;
  CXConfigType cx_config_;
};

}  // namespace tket
