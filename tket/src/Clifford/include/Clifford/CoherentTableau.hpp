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

#include "SymplecticTableau.hpp"

namespace tket {

// Forward declare Circuit for friend converter
class Circuit;

class CoherentTableau {
  /**
   * Represents the stabiliser group for a Clifford process with qubit
   * initialisations/post-selections and mixed initialisations/discards.
   *
   * Full presentation of the data structure and its theory is available via the
   * internal docs.
   *
   * Rows correspond to the coherent stabilisers of the mixed process.
   * Generators are not maintained for the decoherent stabilisers, since they
   * can be derived through commutation with the coherent subgroup.
   *
   * Each row is divided into its input segment and output segment. Recall that
   * a row RxS means SCR^T = C. When mapped to a sparse readable representation,
   * independent QubitPauliTensor objects are used for each segment, so we no
   * longer expect their individual phases to be +-1, instead only requiring
   * this on their product.
   *
   * Columns of the tableau are indexed by pair of Qubit id and a tag to
   * distinguish input vs output. Rows are not maintained in any particular
   * order.
   */
 public:
  enum class TableauSegment { Input, Output };
  typedef std::pair<Qubit, TableauSegment> col_key_t;
  typedef boost::bimap<col_key_t, unsigned> tableau_col_index_t;
  typedef std::pair<QubitPauliTensor, QubitPauliTensor> row_tensor_t;
  /**
   * Construct the tableau for the identity unitary over n qubits (given default
   * qubit names).
   */
  explicit CoherentTableau(unsigned n);
  /**
   * Construct the tableau for the identity unitary over specific qubits.
   */
  explicit CoherentTableau(const qubit_vector_t& qbs);
  /**
   * Construct a tableau from the underling binary matrices.
   * Qubits are given default names and mapped such that the first columns are
   * inputs and the last columns are outputs.
   * @param xmat The X component of the rows in the symplectic representation.
   * @param zmat The Z component of the rows in the symplectic representation.
   * @param phase The phases of the rows.
   */
  explicit CoherentTableau(
      const MatrixXb& xmat, const MatrixXb& zmat, const VectorXb& phase,
      unsigned n_ins = 0);
  /**
   * Construct a tableau directly from its rows.
   * Each row is represented as a product of QubitPauliTensors where the first
   * is over the input qubits and the second is over the outputs..
   */
  explicit CoherentTableau(const std::list<row_tensor_t>& rows);
  /**
   * Other required constructors
   */
  CoherentTableau(const CoherentTableau& other) = default;
  CoherentTableau(CoherentTableau&& other) = default;
  CoherentTableau& operator=(const CoherentTableau& other) = default;
  CoherentTableau& operator=(CoherentTableau&& other) = default;

  /**
   * Get the number of rows in the tableau
   */
  unsigned get_n_rows() const;
  /**
   * Get the number of total qubits/boundaries of the process.
   * The number of boundaries is the sum of inputs and outputs.
   * The number of binary columns in the full symplectic representation is
   * 2*n_boundaries+1.
   */
  unsigned get_n_boundaries() const;
  /**
   * Get the number of boundaries representing inputs to the process.
   */
  unsigned get_n_inputs() const;
  /**
   * Get the number of boundaries representing outputs to the process.
   */
  unsigned get_n_outputs() const;

  /**
   * Get the map between columns and labels, identified as a qubit ID and an
   * indicator for whether it represents an input or output of the process.
   */
  // const tableau_col_index_t& get_col_index() const;

  /**
   * Read off a row as a Pauli string
   */
  row_tensor_t get_row(unsigned i) const;
  /**
   * Combine rows into a single row
   */
  row_tensor_t get_row_product(const std::vector<unsigned>& rows) const;

  /**
   * Transform the tableau according to consuming a Clifford gate at either end
   * of the circuit. Applying a gate to the inputs really consumes its transpose
   * via the CJ isomorphism. The basic gates considered here are all
   * self-transpose, and the method for other Clifford gates handles the
   * transposition internally by reversing the order of application. For
   * multi-qubit gates, the qubits applied must either be all inputs or all
   * outputs.
   */
  void apply_S(const Qubit& qb, TableauSegment seg = TableauSegment::Output);
  void apply_V(const Qubit& qb, TableauSegment seg = TableauSegment::Output);
  void apply_CX(
      const Qubit& control, const Qubit& target,
      TableauSegment seg = TableauSegment::Output);
  void apply_gate(
      OpType type, const qubit_vector_t& qbs,
      TableauSegment seg = TableauSegment::Output);

  /**
   * Transform the tableau according to consuming a Clifford-phase Pauli gadget
   * at either end of the circuit. Transposition for inputs is handled
   * internally. The qubits applied over must either be all inputs or all
   * outputs.
   *
   * @param pauli The string of the Pauli gadget
   * @param half_pis The Clifford angle: {0, 1, 2, 3} represents {0, pi/2, pi,
   * -pi/2}
   */
  void apply_pauli(
      const QubitPauliTensor& pauli, unsigned half_pis,
      TableauSegment seg = TableauSegment::Output);

  /**
   * Post-select a qubit of the Choi state in |0>.
   * Over the input system, this corresponds to initialising in the |0> state.
   * Over the output system, this corresponds to a true post-selection.
   * Throws an exception if post-selection has a zero chance of success.
   */
  void post_select(
      const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Discards a qubit of the Choi state.
   * Over the input system, this corresponds to initialising in the
   * maximally-mixed state. Over the output system, this corresponds to a true
   * discard/partial trace.
   */
  void discard_qubit(
      const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Removes a row from the tableau.
   * The final row is shifted into its place.
   */
  void remove_row(unsigned row);

  /**
   * Combine two tableaux in sequence/parallel.
   * Any matching output qubits of first and input qubits of second will be
   * contracted, and others will be added in parallel. Throws an exception if
   * parallel composition would introduce name clashes.
   *
   * @param first First circuit in the sequential order
   * @param second Second circuit in the sequential order
   */
  static CoherentTableau compose(
      const CoherentTableau& first, const CoherentTableau& second);

  friend CoherentTableau circuit_to_partial_tableau(const Circuit& circ);
  friend Circuit partial_tableau_to_circuit(const CoherentTableau& circ);

  friend void to_json(nlohmann::json& j, const CoherentTableau& tab);
  friend void from_json(const nlohmann::json& j, CoherentTableau& tab);

  friend std::ostream& operator<<(std::ostream& os, const CoherentTableau& tab);

 private:
  /**
   * The actual binary tableau.
   */
  SymplecticTableau tab_;

  /**
   * Map between column indices and the corresponding qubit ID and type.
   */
  tableau_col_index_t col_index_;

  row_tensor_t stab_to_row_tensor(const PauliStabiliser& stab) const;
  PauliStabiliser row_tensor_to_stab(const row_tensor_t& ten) const;

  /**
   * Removes a column from the tableau.
   * The final column is shifted into its place and col_index_ updated.
   * Hidden as private since it just strips this qubit out of the table without
   * attempting to preserve correctness of the stabilizers.
   */
  void remove_col(unsigned col);
};

JSON_DECL(CoherentTableau::TableauSegment)
JSON_DECL(CoherentTableau)

std::ostream& operator<<(std::ostream& os, const CoherentTableau& tab);

}  // namespace tket