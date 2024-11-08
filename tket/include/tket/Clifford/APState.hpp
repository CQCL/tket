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
 * APState class for Clifford states with global phase tracking
 *
 * The "affine with phases" form of a Clifford state from ZX calculus (see
 * Kissinger & van de Wetering, "Picturing Quantum Software") represents n-qubit
 * Clifford states uniquely with the following data:
 * - A binary (k,n) matrix A in reduced row echelon form. The leading columns of
 * each row are referred to as "leading qubits", and the others as "free
 * qubits". We apply H gates to all free qubits, then for each row of A we apply
 * a CX targeted on the leading qubit and control by each other entry in the
 * row.
 * - A binary (k) vector B describing X gates applied to the leading qubits.
 * - A Clifford phase polynomial Phi over the free qubits. Wlog this can be in
 * the form of a symmetric, zero-diagonal, binary (n-k,n-k) matrix E describing
 * CZ gates and a (n-k) vector P of integers mod 4 describing S gates.
 *
 * This gives a canonical statevector:
 * \sum_{x, Ax=B} i^{\Phi(x)} \ket{x}
 *
 * This canonical statevector fixes a reference phase from which we can track
 * global phase with an additional parameter.
 *
 * When applying gates, qubits may switch between being leading and free
 * (including those that aren't involved in the gate due to the need to reduce
 * to the canonical form, e.g. reduced row echelon form for A). The updates are
 * easiest to prove in the ZX calculus for the gateset (CZ, S i.e. pi/2 green
 * phase, SX i.e. pi/2 red phase).
 *
 * To save on resizing the matrices and vectors, we will keep them at full size
 * and just assert correctness that, e.g. every entry in A is either on the
 * diagonal (to indicate the qubit is leading) or in the row of a leading qubit
 * and a column of a later following qubit.
 */
class APState {
 public:
  /**
   * An upper triangular matrix whose diagonal entries indicate leading qubits
   * and off-diagonal entries represent a CX gate from the column qubit
   * (necessarily a free qubit) to the row qubit.
   *
   * If a diagonal entry is 0 (it is a free qubit), then every entry in the row
   * is 0. If a diagonal entry is 1 (it is a leading qubit), all other entries
   * in the column are 0.
   */
  MatrixXb A;

  /**
   * A vector indicating X gates on leading qubits.
   *
   * Must be 0 on every free qubit.
   */
  VectorXb B;

  /**
   * A symmetric, zero diagonal matrix whose entries indicate CZs between free
   * qubits.
   *
   * Every diagonal entry must be 0.
   * The column and row for each leading qubit must be 0.
   */
  MatrixXb E;

  /**
   * A vector indicating S^{P(i)} on qubit i which must be free.
   *
   * Must be 0 on every leading qubit.
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
      const MatrixXb& A_, const VectorXb& B_, const MatrixXb& E_,
      const Eigen::VectorXi& P_, const Expr& phase_);

  /**
   * Constructs the state \ket{0}^{\otimes n_qubits} in AP form.
   */
  APState(unsigned n_qubits);

  /**
   * Constructs the state in AP form from a given statevector.
   */
  APState(const Eigen::VectorXcd& sv);

  /**
   * Verifies the internal correctness of the data structure. Throws an
   * exception if the data structure is invalid.
   */
  void verify();

  /**
   * Calculates the statevector of the state.
   */
  Eigen::VectorXcd to_statevector();

  /**
   * Applies a CZ gate to the state. O(n^2) wrt number of qubits in the state.
   */
  void apply_CZ(unsigned ctrl, unsigned trgt);

  /**
   * Applies an S gate to the state. O(n^2) wrt number of qubits in the state.
   */
  void apply_S(unsigned q);

  /**
   * Applies a V gate to the state. O(n^2) wrt number of qubits in the state.
   */
  void apply_V(unsigned q);

  /**
   * Applies an unparameterised Clifford gate to the state on the chosen qubits.
   * O(n^2) wrt number of qubits in the state.
   */
  void apply_gate(OpType type, const std::vector<unsigned>& qbs);
};

}  // namespace tket