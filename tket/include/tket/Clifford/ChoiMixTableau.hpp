// Copyright 2019-2023 Cambridge Quantum Computing
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

class ChoiMixTableau {
  /**
   * Represents the stabiliser group for a Clifford process with qubit
   * initialisations/post-selections and mixed initialisations/discards.
   *
   * Based on the mixed stabiliser tableau representation of Audenaert & Pleino
   * 2005 (doi:10.1088/1367-2630/7/1/170), commonly used for representing a
   * stabilizer code in QEC. We make use of the Choi-Jamiolkovski isomorphism to
   * represent processes with inputs via the mixed tableau of their Choi state.
   *
   * Rows correspond to generators for the stabilizers of the mixed process,
   * i.e. Pauli measurements with expectation value 1. We will refer to these as
   * "coherent stabilizers" or "coherent subgroup". With mixed states, there may
   * exist additional Pauli operators which leave the state unchanged but do
   * have expectation value 0, which we will refer to as "decoherent
   * stabilizers" (i.e. the logical operators of the corresponding stabiliser
   * code). The extra generators for these are not stored since they can be
   * derived through commutation with the coherent subgroup. Similarly,
   * destabilisers (detectable errors) can be derived. A future implementation
   * may wish to always store the decoherent stabilizers and destabilisers for
   * faster updates.
   *
   * Each row is divided into its input segment and output segment. Under the CJ
   * isomorphism, a row RxS means (in matrix multiplication order) SCR^T = C.
   * When mapped to a sparse readable representation, independent
   * SpPauliStabiliser objects are used for each segment, so we no longer expect
   * their individual phases to be +-1, instead only requiring this on their
   * product. get_row() will automatically transpose the input segment term so
   * it is presented as RxS s.t. SCR = C.
   *
   * Columns of the tableau are indexed by pair of Qubit id and a tag to
   * distinguish input vs output. Rows are not maintained in any particular
   * order, though the `gaussian_form()` method will bring it into row-reduced
   * echelon form.
   */
 public:
  enum class TableauSegment { Input, Output };
  typedef std::pair<Qubit, TableauSegment> col_key_t;
  typedef boost::bimap<col_key_t, unsigned> tableau_col_index_t;
  typedef std::pair<SpPauliStabiliser, SpPauliStabiliser> row_tensor_t;

  /**
   * The actual binary tableau.
   */
  SymplecticTableau tab_;

  /**
   * Map between column indices and the corresponding qubit ID and type.
   */
  tableau_col_index_t col_index_;

  /**
   * Construct the tableau for the identity unitary over n qubits (given default
   * qubit names).
   */
  explicit ChoiMixTableau(unsigned n);
  /**
   * Construct the tableau for the identity unitary over specific qubits.
   */
  explicit ChoiMixTableau(const qubit_vector_t& qbs);
  /**
   * Construct a tableau from the underlying binary matrices.
   * Qubits are given default names and mapped such that the first columns are
   * inputs and the last columns are outputs.
   * @param xmat The X component of the rows in the symplectic representation.
   * @param zmat The Z component of the rows in the symplectic representation.
   * @param phase The phases of the rows.
   * @param n_ins The number of qubits identified as inputs.
   */
  explicit ChoiMixTableau(
      const MatrixXb& xmat, const MatrixXb& zmat, const VectorXb& phase,
      unsigned n_ins = 0);
  /**
   * Construct a tableau directly from its rows.
   * Each row is represented as a product of SpPauliStabilisers where the first
   * is over the input qubits and the second is over the outputs.
   * A row RxS is a pair s.t. SCR = C
   */
  explicit ChoiMixTableau(const std::list<row_tensor_t>& rows);
  /**
   * Other required constructors
   */
  ChoiMixTableau(const ChoiMixTableau& other) = default;
  ChoiMixTableau(ChoiMixTableau&& other) = default;
  ChoiMixTableau& operator=(const ChoiMixTableau& other) = default;
  ChoiMixTableau& operator=(ChoiMixTableau&& other) = default;

  /**
   * Get the number of rows in the tableau
   */
  unsigned get_n_rows() const;
  /**
   * Get the total number of qubits/boundaries of the process.
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
   * Get the number of boundaries representing outputs from the process.
   */
  unsigned get_n_outputs() const;
  /**
   * Get all qubit names present in the input segment.
   */
  qubit_vector_t input_qubits() const;
  /**
   * Get all qubit names present in the output segment.
   */
  qubit_vector_t output_qubits() const;

  /**
   * Read off a row as a Pauli string.
   * Returns a pair of Pauli strings RxS such that SCR = C
   */
  row_tensor_t get_row(unsigned i) const;
  /**
   * Combine rows into a single row.
   * Returns a pair of Pauli strings RxS such that SCR = C
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
  void apply_Z(const Qubit& qb, TableauSegment seg = TableauSegment::Output);
  void apply_V(const Qubit& qb, TableauSegment seg = TableauSegment::Output);
  void apply_X(const Qubit& qb, TableauSegment seg = TableauSegment::Output);
  void apply_H(const Qubit& qb, TableauSegment seg = TableauSegment::Output);
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
   * @param seg Whether to apply the Pauli gadget over the inputs or outputs
   */
  void apply_pauli(
      const SpPauliStabiliser& pauli, unsigned half_pis,
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
   * Performs the effect of an OpType::Collapse gate on the qubit.
   * I.e. non-destructively measure the qubit in the Z basis and ignore the
   * result.
   */
  void collapse_qubit(
      const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  void add_qubit(const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Removes a row from the tableau.
   * The final row is shifted into its place.
   */
  void remove_row(unsigned row);

  /**
   * Permutes columns into a canonical order: chosen segment first, subordered
   * in ILO.
   */
  void canonical_column_order(TableauSegment first = TableauSegment::Input);

  /**
   * Reduces the underlying SymplecticTableau to its Gaussian/row-echelon form.
   */
  void gaussian_form();

  /**
   * Renames qubits.
   */
  void rename_qubits(const qubit_map_t& qmap, TableauSegment seg);

  /**
   * Check whether the process described by the tableau is a unitary
   */
  bool is_unitary() const;

  /**
   * Combine two tableaux in sequence/parallel.
   * Any matching output qubits of first and input qubits of second will be
   * contracted, and others will be added in parallel. Throws an exception if
   * parallel composition would introduce name clashes.
   *
   * @param first First circuit in the sequential order
   * @param second Second circuit in the sequential order
   */
  static ChoiMixTableau compose(
      const ChoiMixTableau& first, const ChoiMixTableau& second);

  friend std::ostream& operator<<(std::ostream& os, const ChoiMixTableau& tab);
  bool operator==(const ChoiMixTableau& other) const;

 private:
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

JSON_DECL(ChoiMixTableau::TableauSegment)
JSON_DECL(ChoiMixTableau)

std::ostream& operator<<(std::ostream& os, const ChoiMixTableau& tab);

}  // namespace tket