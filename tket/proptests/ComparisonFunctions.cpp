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

#include "ComparisonFunctions.hpp"

#include <sstream>
#include <stdexcept>

// NOTE: this is an identical copy of a file in tket-tests, in the Simulation
// folder. This is deliberate! Please keep both in sync!
// See the tket-tests file for a mathematical discussion.

namespace tket {
namespace tket_sim {

// Checks that they're both column vectors, or both square, of the same size.
// Throws if not.
static void throw_if_invalid_sizes(
    const Eigen::MatrixXcd& m1, const Eigen::MatrixXcd& m2) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    throw std::runtime_error("Different sized matrices");
  }
  // Check that it has 2^n rows for some n. Will throw if not.
  (void)get_number_of_qubits(m1.rows());
  if (m1.rows() == m1.cols() || m1.cols() == 1) {
    // Square, or a column vector.
    return;
  }
  throw std::runtime_error("Not square, and also not column vectors");
}

// Note: these should NOT be asserts, because it is conceivable
// that a really deep circuit could be tested, with so many gates
// that the numerical errors build up and make matrices
// which are not almost unitary.
static void throw_if_not_unitary_or_unnormalised_statevector(
    const Eigen::MatrixXcd& m, double tolerance) {
  const Eigen::MatrixXcd product = m.adjoint() * m;
  if (product.isApprox(
          Eigen::MatrixXcd::Identity(product.rows(), product.rows()),
          tolerance)) {
    return;
  }
  if (product.rows() == 1) {
    // Of course, for 0-qubit circuits there's no distinction between
    // state vectors and 1x1 unitaries! Don't worry about this.
    throw std::runtime_error("State vector is not normalised");
  }
  throw std::runtime_error("Matrix is not unitary");
}

static bool compare_statevectors_or_unitaries_may_throw(
    const Eigen::MatrixXcd& m1, const Eigen::MatrixXcd& m2,
    MatrixEquivalence equivalence, double tolerance) {
  throw_if_invalid_sizes(m1, m2);
  throw_if_not_unitary_or_unnormalised_statevector(m1, tolerance);
  throw_if_not_unitary_or_unnormalised_statevector(m2, tolerance);
  if (equivalence == MatrixEquivalence::EQUAL) {
    return m1.isApprox(m2, tolerance);
  }

  // We allow equivalence only up to global phase.
  // We now know that U,V are EITHER almost unitary,
  // OR almost norm one column vectors.
  // See the above mathematical discussion:
  // if A = cB for some |c|=1, then  (A adj)B = (c*)(B adj)B = (c*)Id,
  // where Id may be 1x1.
  //
  // Thus (U adj)V will be approximately diagonal, with diagonal entries
  // almost equal to each other.
  const Eigen::MatrixXcd product = m1.adjoint() * m2;
  auto entry = product(0, 0);
  const double entry_abs = std::abs(entry);
  if (!(std::abs(entry_abs - 1.0) < tolerance)) {
    // Catch NaNs also! Although they should already have been
    // caught above in the unitary/norm one checks.
    return false;
  }
  const auto size = product.rows();
  if (size == 1) {
    return true;
  }
  // Shouldn't make much difference but do it anyway.
  entry /= entry_abs;
  return product.isApprox(
      entry * Eigen::MatrixXcd::Identity(size, size), tolerance);
}

bool compare_statevectors_or_unitaries(
    const Eigen::MatrixXcd& m1, const Eigen::MatrixXcd& m2,
    MatrixEquivalence equivalence, double tolerance) {
  try {
    return compare_statevectors_or_unitaries_may_throw(
        m1, m2, equivalence, tolerance);
  } catch (const std::exception& e) {
    std::stringstream ss;
    ss << "Input matrices have sizes (" << m1.rows() << "," << m1.cols()
       << ") and (" << m2.rows() << "," << m2.cols() << "). tol=" << tolerance
       << " : " << e.what();
    throw std::runtime_error(ss.str());
  }
}

}  // namespace tket_sim
}  // namespace tket
