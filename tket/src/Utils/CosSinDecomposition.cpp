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

#include "CosSinDecomposition.hpp"

#include <cmath>

#include "Constants.hpp"
#include "EigenConfig.hpp"
#include "MatrixAnalysis.hpp"
#include "TketLog.hpp"

namespace tket {

csd_t CS_decomp(const Eigen::MatrixXcd &u) {
  if (!is_unitary(u)) {
    throw std::invalid_argument("Matrix for CS decomposition is not unitary");
  }
  unsigned N = u.rows();
  if (N % 2 != 0) {
    throw std::invalid_argument(
        "Matrix for CS decomposition has odd dimensions");
  }
  unsigned n = N / 2;
  Eigen::MatrixXcd u00 = u.topLeftCorner(n, n);
  Eigen::MatrixXcd u01 = u.topRightCorner(n, n);
  Eigen::MatrixXcd u10 = u.bottomLeftCorner(n, n);
  Eigen::MatrixXcd u11 = u.bottomRightCorner(n, n);

  Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::NoQRPreconditioner> svd(
      u00, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXcd l0 = svd.matrixU().rowwise().reverse();
  Eigen::MatrixXcd r0_dag = svd.matrixV().rowwise().reverse();
  Eigen::MatrixXd c = svd.singularValues().reverse().asDiagonal();
  Eigen::MatrixXcd r0 = r0_dag.adjoint();

  // Now u00 = l0 c r0; l0 and r0 are unitary, and c is diagonal with positive
  // non-decreasing entries. Because u00 is a submatrix of a unitary matrix, its
  // singular values (the entries of c) are all <= 1.

  Eigen::HouseholderQR<Eigen::MatrixXcd> qr = (u10 * r0_dag).householderQr();
  Eigen::MatrixXcd l1 = qr.householderQ();
  Eigen::MatrixXcd S = qr.matrixQR().triangularView<Eigen::Upper>();

  // Now u10 r0* = l1 S; l1 is unitary, and S is upper triangular.
  //
  // Claim: S is diagonal.
  // Proof: Since u is unitary, we have
  //     I = u00* u00 + u10* u10
  //       = (l0 c r0)* (l0 c r0) + (l1 S r0)* (l1 S r0)
  //       = r0* c l0* l0 c r0 + r0* S* l1* l1 S r0
  //       = r0* c^2 r0 + r0* S* S r0
  //       = r0* (c^2 + S* S) r0
  // So I = c^2 + (S* S), so (S* S) = I - c^2 is a diagonal matrix with non-
  // increasing entries in the range [0,1). As S is upper triangular, this
  // implies that S must be diagonal. (Proof by induction over the dimension n
  // of S: consider the two cases S_00 = 0 and S_00 != 0 and reduce to the n-1
  // case.)
  //
  // We want S to be real. This is not guaranteed, though since it is diagonal
  // it can be made so by adjusting l1. In fact, it seems that Eigen always does
  // give us a real S. We will not assume this, but will log an informational
  // message and do the adjustment ourselves if that behaviour does change.

  if (!S.imag().isZero()) {
    tket_log()->info(
        "Eigen surprisingly returned a non-real diagonal R in QR "
        "decomposition; adjusting Q and R to make it real.");
    for (unsigned j = 0; j < n; j++) {
      std::complex<double> z = S(j, j);
      double r = std::abs(z);
      if (r > EPS) {
        std::complex<double> w = std::conj(z) / r;
        S(j, j) *= w;
        l1.col(j) /= w;
      }
    }
  }

  // Now S is real and diagonal, and c^2 + S^2 = I.

  Eigen::MatrixXd s = S.real();

  // Make all entries in s non-negative.
  for (unsigned j = 0; j < n; j++) {
    if (s(j, j) < 0) {
      s(j, j) = -s(j, j);
      l1.col(j) = -l1.col(j);
    }
  }

  // Finally compute r1, being careful not to divide by small things.
  Eigen::MatrixXcd r1 = Eigen::MatrixXcd::Zero(n, n);
  for (unsigned i = 0; i < n; i++) {
    if (s(i, i) > c(i, i)) {
      r1.row(i) = -(l0.adjoint() * u01).row(i) / s(i, i);
    } else {
      r1.row(i) = (l1.adjoint() * u11).row(i) / c(i, i);
    }
  }
  return {l0, l1, r0, r1, c, s};
}

}  // namespace tket
