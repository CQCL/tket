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

#include "AssertionSynthesis.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>

#include "CircUtils.hpp"
#include "Circuit.hpp"
#include "OpType/OpType.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Constants.hpp"
#include "Utils/CosSinDecomposition.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

/**
 * Diagonalise a projector matrix P such that P = U*D*U.dag, where
 * U is a unitary and D is a diagonal binary matrix.
 * The resulting diagonal matrix always have 1s on the left.
 * The benefit of permuting 1s to the left is that the
 * diagonal matrix can always be factorised into |0><0|, |1><1| or I if the
 * projector has a valid rank. However, this is not optimised.
 *
 * @param P projector matrix
 *
 * @return D as a boolean vector, U, rank
 */
static std::tuple<VectorXb, Eigen::MatrixXcd, int> projector_diagonalisation(
    const Eigen::MatrixXcd &P) {
  // Solve the eigen values
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigen_solver(P);
  // Make a permutation matrix to reorder D and U_dag
  Eigen::PermutationMatrix<Eigen::Dynamic> perm(P.rows());
  Eigen::MatrixXcd P1 = eigen_solver.eigenvectors() *
                        eigen_solver.eigenvalues().asDiagonal() *
                        eigen_solver.eigenvectors().adjoint();
  TKET_ASSERT(P1.isApprox(P));
  // Cast eigenvalues to booleans
  VectorXb eigenvalues(eigen_solver.eigenvalues().rows());
  for (unsigned i = 0; i < eigenvalues.rows(); i++) {
    if (std::abs(eigen_solver.eigenvalues()[i]) < EPS) {
      eigenvalues[i] = 0;
    } else {
      eigenvalues[i] = 1;
    }
  }
  // Move 1 to the front
  unsigned j = 0;
  for (unsigned i = 0; i < P.rows(); i++) {
    if (eigenvalues[i] == 1) {
      perm.indices()[j] = i;
      j++;
    }
  }
  // Assign rank
  unsigned rank = j;
  // Pad 0 to the end
  for (unsigned i = 0; i < P.rows(); i++) {
    if (eigenvalues[i] == 0) {
      perm.indices()[j] = i;
      j++;
    }
  }
  VectorXb D = perm.transpose() * eigenvalues;
  Eigen::MatrixXcd U = eigen_solver.eigenvectors() * perm;
  TKET_ASSERT((U * D.cast<double>().asDiagonal() * U.adjoint()).isApprox(P));
  return {D, U, rank};
}

// Get the n_th bit of an integer
static bool get_bit(unsigned value, unsigned index) {
  return (value & (1 << index)) >> index;
}

/**
 * Tensor-factorize a binary diagonal matrix with dimension 2^n.
 * Assume D has a rank <= 2^(n-1), and the rank is power of 2
 *
 * @param D a vector representation of a diagonal binary matrix
 *
 * @return A list of 2x2 diagonal binary matrices D_0,...,D_m (as vectors)
 *  such that tensor(D_0,...,D_m) == D
 */
static std::vector<Vector2b> tensor_factorization(const VectorXb &D) {
  unsigned dim = D.rows();
  unsigned log_dim = (unsigned)log2(dim);
  std::vector<Vector2b> factorisations;
  for (unsigned i = 0; i < log_dim; i++) {
    factorisations.push_back(Vector2b::Zero());
  }
  // The diagonal matrix D can be factorised into a tensorproduct of
  // log_dim 2x2 diagonal boolean matrices in the set {\0><0|, |1><1|, I}
  // The factorisation can be determined by the diagonals in D
  // For example, if D[0,0] == 1, then all the 2x2 matices should have 1 in
  // their top left entry.
  for (unsigned i = 0; i < dim; i++) {
    if (D(i) == 1) {
      for (unsigned j = 0; j < log_dim; j++) {
        if (get_bit(i, j)) {
          factorisations[j](1) = 1;
        } else {
          factorisations[j](0) = 1;
        }
      }
    }
  }
  return factorisations;
}

static void apply_unitary(Circuit *circ, const Eigen::MatrixXcd &U) {
  if (U.rows() == 2) {
    // Add a Unitary1qBox
    Unitary1qBox ubox(U);
    circ->add_op<unsigned>(std::make_shared<Unitary1qBox>(ubox), {0});
  } else if (U.rows() == 4) {
    // Add a Unitary2qBox
    Unitary2qBox ubox(U);
    circ->add_op<unsigned>(std::make_shared<Unitary2qBox>(ubox), {0, 1});
  } else if (U.rows() == 8) {
    // Add a Unitary3qBox
    Unitary3qBox ubox(U);
    circ->add_op<unsigned>(std::make_shared<Unitary3qBox>(ubox), {0, 1, 2});

  } else {
    throw CircuitInvalidity("Only 2x2, 4x4, and 8x8 projectors are supported");
  }
}

/**
 * Apply measurements and update the expectations
 *
 * @param circ the circuit to apply measurements to
 * @param z_projectors a list of 2x2 diagonal binary matrices
 * @param debug_bit_index the next available index of the default classical
 * register
 * @param expected_readouts a list of expectations. e.g.
 * expected_readouts[i]=True indicates that the ith bit in the default register
 * expects the readout to be 1.
 *
 * @return the next available index of the default classical register
 */
static unsigned apply_z_measurements(
    Circuit &circ, const std::vector<Vector2b> z_projectors,
    unsigned debug_bit_index, std::vector<bool> &expected_readouts) {
  for (unsigned i = 0; i < z_projectors.size(); i++) {
    if (z_projectors[i](0) == 1 && z_projectors[i](1) == 1) {
      continue;
    }
    Qubit q(i);
    Bit b(debug_bit_index++);
    circ.add_bit(b);
    circ.add_op<UnitID>(OpType::Measure, {q, b});
    const bool is_1_projector =
        z_projectors[i](0) == 0 && z_projectors[i](1) == 1;
    expected_readouts.push_back(is_1_projector);
  }
  return debug_bit_index;
}

/**
 * Split a projector into two 2^(n-1)-ranked projectors that project the same
 * subspace when combined. Assume D has 1s on the left. Suppose D.rows() = 16,
 * rank = 5, then the diagonal of D is 11111000|00000000. D can't be factorised
 * into Z basis projectors, so we create two rank=2^(n-1) projectors that
 * project the same subspace when combined (intersection of two subspaces).
 * e.g. D1: 11111111|00000000, D2: 11111000|11100000.
 * After permuting D2, both diagonal matrices should have diagonal entires:
 * 11111111|00000000.
 * Similar to the projector_diagonalisation function, the permutation of the
 * eigen decomposition is not optimised.
 *
 * @param D a vector representation of a diagonal binary matrix
 * @param U a unitary matrix
 * @param rank the rank of D
 *
 * @return D1, U1, D2, U2 such that the intersection of the subspaces projected
 * by P1=U1*D1*U1.dag and P2=U2*D2*U2.dag is the subspace projected by
 * P=U*D*U.dag.
 */
static std::tuple<VectorXb, Eigen::MatrixXcd, VectorXb, Eigen::MatrixXcd>
projector_split(
    const VectorXb &D, const Eigen::MatrixXcd &U, const unsigned &rank) {
  VectorXb D_new(D);
  for (unsigned i = rank; i < D_new.rows() / 2; i++) {
    D_new(i) = 1;
  }
  // U1
  Eigen::MatrixXcd U1(U);
  // U2
  Eigen::PermutationMatrix<Eigen::Dynamic> perm(D_new.rows());
  for (unsigned i = 0; i < D_new.rows(); i++) {
    if (i >= rank && i < D_new.rows() / 2) {
      perm.indices()[i] = i + D_new.rows() / 2 - rank;
    } else if (i >= D_new.rows() / 2 && i < D_new.rows() - rank) {
      perm.indices()[i] = i - (D_new.rows() / 2 - rank);
    } else {
      perm.indices()[i] = i;
    }
  }
  Eigen::MatrixXcd U2 = U * perm;
  return {D_new, U1, D_new, U2};
}

std::tuple<Circuit, std::vector<bool>> projector_assertion_synthesis(
    const Eigen::MatrixXcd &P) {
  // Diagonoalise the projector P
  Eigen::MatrixXcd projector(P);
  unsigned n_qubits = (unsigned)log2(projector.rows());
  VectorXb D;
  Eigen::MatrixXcd U;
  unsigned rank;
  std::tie(D, U, rank) = projector_diagonalisation(projector);
  if (rank == 0) {
    throw CircuitInvalidity("The projector must have none-zero rank");
  }
  std::vector<bool> expected_readouts;
  Circuit circ(n_qubits);
  if (rank > projector.rows() / 2) {
    if (n_qubits >= 3) {
      throw CircuitInvalidity(
          "8x8 projector that requires an ancilla is not supported");
    }
    // Add an auxilary qubit to the end
    Qubit ancilla(n_qubits);
    circ.add_qubit(ancilla, false);
    circ.add_op<Qubit>(OpType::Reset, {ancilla});
    Eigen::MatrixXcd zero_projector = Eigen::MatrixXcd::Zero(2, 2);
    zero_projector(0, 0) = 1;
    Eigen::MatrixXcd new_projector =
        Eigen::KroneckerProduct(projector, zero_projector);
    std::tie(D, U, rank) = projector_diagonalisation(new_projector);
  }

  // Initialise debug bit index
  unsigned debug_bit_index = 0;
  // Check if the rank is power of 2
  // If a number is power of 2 then the count of binary 1s will be 1.
  if (rank != 0 && (rank & (rank - 1)) == 0) {
    // Implement the projection
    // Apply the transformation before the measurements
    std::vector<Vector2b> tensors = tensor_factorization(D);
    apply_unitary(&circ, U.adjoint());
    apply_z_measurements(circ, tensors, debug_bit_index, expected_readouts);
    apply_unitary(&circ, U);
  }
  // Transform this into two projectors
  else {
    auto [D1, U1, D2, U2] = projector_split(D, U, rank);
    std::vector<Vector2b> tensors1 = tensor_factorization(D1);
    apply_unitary(&circ, U1.adjoint());
    debug_bit_index = apply_z_measurements(
        circ, tensors1, debug_bit_index, expected_readouts);
    apply_unitary(&circ, U1);

    std::vector<Vector2b> tensors2 = tensor_factorization(D2);
    apply_unitary(&circ, U2.adjoint());
    apply_z_measurements(circ, tensors2, debug_bit_index, expected_readouts);
    apply_unitary(&circ, U2);
  }
  return {circ, expected_readouts};
}

static unsigned get_n_qubits_from_stabilisers(
    const PauliStabiliserList &paulis) {
  if (paulis.size() == 0) {
    throw CircuitInvalidity("Stabilisers cannot be empty");
  }
  unsigned stabiliser_len = paulis[0].string.size();
  for (unsigned i = 1; i < paulis.size(); i++) {
    if (paulis[i].string.size() != stabiliser_len) {
      throw CircuitInvalidity("Stabilisers have unequal lengths");
    }
  }
  return stabiliser_len;
}

/**
 * Apply assertion operators and update the expectations
 *
 * @param circ the circuit to apply assertion operators to
 * @param pauli a Pauli stabiliser
 * @param debug_bit_index the next available index of the default classical
 * register
 * @param expected_readouts a list of expectations. e.g.
 * expected_readouts[i]=True indicates that the ith bit in the default register
 * expects the readout to be 1.
 *
 * @return the next available index of the default classical register
 */
static unsigned add_assertion_operator(
    Circuit &circ, const PauliStabiliser &pauli, unsigned debug_bit_index,
    const Qubit &ancilla, std::vector<bool> &expected_readouts) {
  circ.add_op<Qubit>(OpType::Reset, {ancilla});
  circ.add_op<Qubit>(OpType::H, {ancilla});
  for (unsigned i = 0; i < pauli.string.size(); i++) {
    switch (pauli.string[i]) {
      case Pauli::I:
        break;
      case Pauli::X:
        circ.add_op<Qubit>(OpType::CX, {ancilla, Qubit(i)});
        break;
      case Pauli::Y:
        circ.add_op<Qubit>(OpType::CY, {ancilla, Qubit(i)});
        break;
      case Pauli::Z:
        circ.add_op<Qubit>(OpType::CZ, {ancilla, Qubit(i)});
        break;
    }
  }
  circ.add_op<Qubit>(OpType::H, {ancilla});
  Bit b(debug_bit_index++);
  circ.add_bit(b);
  circ.add_op<UnitID>(OpType::Measure, {ancilla, b});
  expected_readouts.push_back(!pauli.coeff);
  return debug_bit_index;
}

// Assume all Pauli stabilisers have equal lengths and no identity
std::tuple<Circuit, std::vector<bool>> stabiliser_assertion_synthesis(
    const PauliStabiliserList &paulis) {
  unsigned n_qubits = get_n_qubits_from_stabilisers(paulis);
  std::vector<bool> expected_readouts;
  Circuit circ(n_qubits);

  // Initialise debug bit index
  unsigned debug_bit_index = 0;
  // Add an ancilla
  Qubit ancilla(n_qubits);
  circ.add_qubit(ancilla, false);

  for (auto &pauli : paulis) {
    debug_bit_index = add_assertion_operator(
        circ, pauli, debug_bit_index, ancilla, expected_readouts);
  }

  return {circ, expected_readouts};
}

}  // namespace tket
