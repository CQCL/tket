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

// NOTE: this is used by both tket-proptests and tket-tests.
// Since it's test-only code, it is NOT placed in tket.
// However, we don't want to make proptests depend on tket-tests,
// or make another package for both to depend on.
// Therefore, we make copies of this here and in proptests.
// If you update one, please update the other!

// Numerical discussion on "compare_statevectors_or_unitaries":
//
// For complex column vectors a,b:
//
// Let  a,b  be the (unknown) exact answers, from the two circuits.
//
// (Thus they are state vectors, i.e.  a = Up  and  b = Vp
// where U,V are the unitary matrices of the two circuits,
// and  p  is some norm one column vector. We really want to
// test if U=V, or U=cV for some complex  c  with  |c|=1,  but we perhaps
// cannot, because U,V might be large and expensive to compute.
// However,  Up, Vp  are cheaper to compute. Often, we use the standard
// computational basis with  p = (1 0 0 0 ...) transposed,
// but this is not necessary).
//
// Let  da,db  be the vectors of numerical errors.
// Let's assume that  a,b  really are the statevectors
// of two circuits which are equivalent, but only up to global phase.
//
// Thus a = cb  for some complex number c with |c|=1, and ||a||=||b||=1.
//
// Let  u = a+da,  v = b+db.
//
// Thus,  u,v  are the known calculated state vectors which are passed
// into the function, since numerical errors have unavoidably been added.
// Thus  ||da||, ||db|| < eps  (but we don't know a,b,da,db).
//
// (Of course, in practice ||da||, ||db|| would grow with N, the dimension,
// but we'll ignore this since N never gets very large.
// In most cases N<2^10).
//
// (In practice eps is quite small, but quite a bit larger
// than the actual likely roundoff error. We don't really care much
// about the distinction between eps, 5.eps, etc. and write  O(eps)  roughly
// to mean a not-too-large multiple of eps).
//
// Now  <u,v> = c + O(eps),  and hence  |<u,v>| = 1 + O(eps).
//
// Notice that if we divide u by ||u|| = 1 + O(eps), etc. to normalise,
// this leads to a product of 1 + O(eps) terms and so we still get
// |<u,v>| = 1 + O(eps), so there is little to be gained by extra normalising.
//
// Thus, we simply test   | |<u,v>| - 1 | < eps.
// If it fails this test, we are highly confident that a,b weren't EXACTLY
// equivalent, because then it WOULD have passed the test with that eps
// (if we are confident that eps is a good upper bound for numerical errors).
//
// But if it passes this test, how confident should we be that a,b
// ARE equivalent? I.e., how different can two circuits be whilst still being
// within roundoff error of the known  u,v?
// This seems surprisingly tricky and subtle, but we ignore the problem here.
// It seems very unlikely that routines designed to preserve unitaries exactly
// would have bugs producing two inequivalent circuits, but with unitary
// matrices so close to each other that they falsely pass the test.
// (Of course we're testing  a,b  which are the products of the unitary
// matrices with vectors, rather than the matrices directly).
//
// For unitary matrices, the argument is similar:
//
// Assume that A,B are the (unknown) EXACT unitary matrices
// of the circuits, but we have been GIVEN matrices
//
//   U = A+dA,  V = B+dB  with ||dA||,||dB|| = O(eps).
//
// We also have  (U adj)U = I + O(eps),  and similarly for V
// (i.e., due to roundoff, U,V are not exactly unitary).
//
// If we had  A = cB  for some complex c with |c|=1, then we would have
//
//   (U adj)V = (A adj)B + O(eps) = (c*)I + O(eps),
//
// since  A,B  are exactly unitary.
// Thus (U adj)V is a nearly diagonal matrix, with diagonal entries
// nearly constant (i.e., all within O(eps) of each other).
//

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
