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

#include "tket/OpType/OpType.hpp"
#include "tket/Utils/MatrixAnalysis.hpp"

namespace tket {

/**
 * APState class, giving a unique form for a (possibly mixed) Clifford state,
 * through which we can track global phase when appropriate.
 *
 * The "affine with phases" form of a Clifford state from ZX calculus (see
 * Kissinger & van de Wetering, "Picturing Quantum Software") represents n-qubit
 * Clifford states uniquely with the following data:
 * - A binary (n,n) matrix A.
 * - A binary (n) vector B.
 * - A symmetric, zero-diagonal, binary (n,n) matrix E.
 * - A n vector P of integers mod 4 describing S gates.
 *
 * This gives a canonical statevector (up to a normalisation scalar):
 * \f[\sum_{x, Ax=B} i^{\Phi(x)} |x>\f]
 *
 * This canonical statevector fixes a reference phase from which we can track
 * global phase with an additional parameter.
 *
 * We can generalise this to mixed qubit states by adding a binary (n,n) matrix
 * C, with which we now have the canonical density matrix (up to a normalisation
 * scalar):
 * \f[\sum_{x_1, x_2; Ax_1 = B = Ax_2, Cx_1 = Cx_2}
 *         i^{x_1^T E x_1 + P^T x_1} |x_1><x_2| (-i)^{x_2^T E x_2 + P^T x_2}\f]
 * where the inner products are calculated in Z4.
 *
 * We can encode this via the ZX-calculus as:
 * - A green spider for each qubit q, with phase given by P(q)pi/2.
 * - The green spiders are connected via Hadamard edges according to E.
 * - For each i, a red spider with phase B(i)pi, connected to the green spiders
 *   according to the row A(i, -).
 * - For each j, a discard connected to a red spider, connected to the green
 *   spiders according to the row C(j, -).
 *
 * Several operations on this data leave the state unchanged:
 * - We can freely add rows simultaneously in A and B whilst preserving the
 *   state, and similarly combine rows in C or add rows from A to C.
 * - If a green spider q is connected to precisely one red spider (a single 1
 *   appears in the column A(-, q), and nothing in column C(-, q)),
 *   transformations exist that allow us to set P(q), E(q, -), and E(-, q) to
 *   zero, by some sequence of bialgebra and local complementation rules.
 * - We can perform a local complementation around any discard j - concerning
 *   its neighbours (those q with C(j, q) = 1), we may increment each P(q) by 1
 *   and flip each E(q, q').
 * However, whilst a normal form exists (generalising "reduced AP-form"),
 * attempting to maintain that form massively increases the number of cases we
 * need to consider for each gate application. We instead just offer a method
 * that reduces a given APState to the normal form which can be used for
 * comparison checks.
 *
 * To mimic all behaviour that ChoiMixTableau supports, we consider the gate set
 * {CZ, S, V, Init, PostSelect, Collapse}, where V is the only difficult one to
 * derive (always comes down to a sequence of local complementations in the ZX
 * picture).
 */
class APState {
 public:
  /**
   * A binary (n,n) matrix used in describing the subspace of computational
   * basis states in the support of the state.
   */
  MatrixXb A;

  /**
   * A binary n vector used in describing the subspace of computaitonal basis
   * states in the support of the state.
   */
  VectorXb B;

  /**
   * A binary (n,n) matrix describing the incoherent constraints, relating the
   * ket and bra sides of the density matrix representation of the state.
   */
  MatrixXb C;

  /**
   * A symmetric, zero diagonal matrix whose entries indicate CZs between
   * qubits.
   */
  MatrixXb E;

  /**
   * A vector indicating S^{P(i)} on qubit i.
   */
  Eigen::VectorXi P;

  /**
   * The global phase term (in half-turns).
   */
  Expr phase;

  /**
   * Constructs a state in AP form from given data.
   */
  APState(
      const MatrixXb& A_, const VectorXb& B_, const MatrixXb& C_,
      const MatrixXb& E_, const Eigen::VectorXi& P_, const Expr& phase_);

  /**
   * Constructs the state |0>^{x n_qubits} in AP form.
   */
  APState(unsigned n_qubits);

  /**
   * Constructs the state in AP form from a given statevector.
   */
  APState(const Eigen::VectorXcd& sv);

  /**
   * Constructs the state in AP form from a given density matrix.
   */
  APState(const Eigen::MatrixXcd& dm);

  /**
   * Other required constructors
   */
  APState(const APState& other) = default;
  APState(APState&& other) = default;
  APState& operator=(const APState& other) = default;
  APState& operator=(APState& other) = default;

  /**
   * Verifies the internal correctness of the data structure. Throws an
   * exception if the data structure is invalid.
   */
  void verify() const;

  /**
   * Calculates the statevector of the state.
   *
   * Throws an exception if C is non-zero (i.e. it represents a mixed state).
   */
  Eigen::VectorXcd to_statevector() const;

  /**
   * Calculates the density matrix of the state.
   */
  Eigen::MatrixXcd to_density_matrix() const;

  /**
   * Applies a CZ gate to the state. O(1).
   */
  void apply_CZ(unsigned ctrl, unsigned trgt);

  /**
   * Applies an S gate to the state. O(1).
   */
  void apply_S(unsigned q);

  /**
   * Applies a V gate to the state. O(n^2) wrt number of qubits in the state.
   */
  void apply_V(unsigned q);

  /**
   * Applies an X gate to the state, faster than applying V twice. O(n) wrt
   * number of qubits in the state.
   */
  void apply_X(unsigned q);

  /**
   * Applies an unparameterised Clifford gate to the state on the chosen qubits.
   * O(n^2) wrt number of qubits in the state.
   */
  void apply_gate(OpType type, const std::vector<unsigned>& qbs);

  /**
   * Initialises a new qubit in the |0> state. Returns the index of the new
   * qubit. O(n) wrt number of qubits in the state.
   */
  unsigned init_qubit();

  /**
   * Post-selects the chosen qubit to <0|, removing it from the state. Moves the
   * final qubit into the place of q, returning its old index. O(n^3) wrt
   * number of qubits in the state.
   */
  unsigned post_select(unsigned q);

  /**
   * Collapses the given qubit in the Z basis. O(n^3) wrt number of qubits in
   * the state.
   */
  void collapse_qubit(unsigned q);

  /**
   * Discards the given qubit, removing it from the state. Moves the final qubit
   * into the place of q, returning its old index. O(n^3) wrt number of qubits
   * in the state.
   */
  unsigned discard_qubit(unsigned q);

  /**
   * Manipulates the state into the following normal form:
   * - A is in reduced row-echelon form. We refer to the leading columns of each
   * row of A as "leading qubits".
   * - C is in reduced row-echelon form, and furthermore is zero in any column
   * of a leading qubit. We refer to the leading columns of each row of C as
   * "mixed qubits".
   * - Each entry of E is either between a mixed qubit and a free qubit, or
   * between two free qubits.
   * - For each leading or mixed qubit, the index of P is zero.
   *
   * Takes time O(n^4) wrt number of qubits in the state.
   */
  void normal_form();

  bool operator==(const APState& other) const;
};

/**
 * A wrapper for APState to provide the same interface as ChoiMixTableau wrt
 * qubit indexing and distinction between input vs output.
 *
 * When applying gates, the methods automatically transpose anything occurring
 * over the input subspace according to the CJ-isomorphism.
 */
class ChoiAPState {
 public:
  enum class TableauSegment { Input, Output };
  typedef std::pair<Qubit, TableauSegment> col_key_t;
  typedef boost::bimap<col_key_t, unsigned> tableau_col_index_t;

  /**
   * The internal APState.
   */
  APState ap_;

  /**
   * Map between column indices and the corresponding qubit ID and type.
   */
  tableau_col_index_t col_index_;

  /**
   * Construct the identity unitary over n qubits (given default qubit names).
   */
  explicit ChoiAPState(unsigned n);
  /**
   * Construct the identity unitary over specific qubits.
   */
  explicit ChoiAPState(const qubit_vector_t& qbs);
  /**
   * Construct the APState from the underlying binary matrices.
   * Qubits are given default names and mapped such that the first columns are
   * inputs and the last columns are outputs.
   */
  explicit ChoiAPState(
      const MatrixXb& A, const VectorXb& b, const MatrixXb& C,
      const MatrixXb& E, const Eigen::VectorXi& P, const Expr& phase = 0,
      unsigned n_ins = 0);
  /**
   * Other required constructors
   */
  ChoiAPState(const ChoiAPState& other) = default;
  ChoiAPState(ChoiAPState&& other) = default;
  ChoiAPState& operator=(const ChoiAPState& other) = default;
  ChoiAPState& operator=(ChoiAPState& other) = default;

  /**
   * Get the total number of qubits/boundaries (sum of inputs and outputs) of
   * the process.
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
   * Get all qubit names present in the input segment.
   */
  qubit_vector_t input_qubits() const;

  /**
   * Get all qubit names present in the output segment.
   */
  qubit_vector_t output_qubits() const;

  /**
   * Transforms the state according to consuming a Clifford gate at either end
   * of the circuit. Applying a gate to the inputs really consumes its transpose
   * via the CJ isomorphism. This method handles the transposition internally by
   * reversing the order of basic self-transpose gates. For multi-qubit gates,
   * the qubits applied must either be all inputs or all outputs.
   */
  void apply_gate(
      OpType type, const qubit_vector_t& qbs,
      TableauSegment seg = TableauSegment::Output);

  /**
   * Initialises a new qubit in the |0> state (or <0| for inputs). O(n) wrt
   * number of qubits in the state.
   */
  void init_qubit(const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Post-selects the chosen qubit to <0| (or applies |0> on inputs), removing
   * it from the state. This DOES NOT CHECK whether or not the post-selection
   * would be successful (e.g. we can post-select <0| on a |1> state). O(n^3)
   * wrt number of qubits in the state.
   */
  void post_select(
      const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Discards the given qubit, removing it from the state. O(n^3) wrt number of
   * qubits in the state.
   */
  void discard_qubit(
      const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Collapses the given qubit in the Z basis. O(n^3) wrt number of qubits in
   * the state.
   */
  void collapse_qubit(
      const Qubit& qb, TableauSegment seg = TableauSegment::Output);

  /**
   * Permutes columns into a canonical order: chosen segment first, subordered
   * in ILO.
   */
  void canonical_column_order(TableauSegment first = TableauSegment::Input);

  /**
   * Reduces the underlying APState to its normal form wrt the current ordering
   * of qubits in its columns.
   */
  void normal_form();

  /**
   * Renames qubits.
   */
  void rename_qubits(const qubit_map_t& qmap, TableauSegment seg);

  bool operator==(const ChoiAPState& other) const;
};

JSON_DECL(APState)
JSON_DECL(ChoiAPState::TableauSegment)
JSON_DECL(ChoiAPState)

}  // namespace tket