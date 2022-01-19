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

#include <unsupported/Eigen/KroneckerProduct>

#include "GateUnitaryMatrixImplementations.hpp"
#include "Utils/Constants.hpp"

// This file is for gates with unitary matrices which serve as
// building blocks for others ("composite gates").
// Of course, this choice is subjective and could change.

namespace tket {
namespace internal {

Eigen::Matrix2cd GateUnitaryMatrixImplementations::Rx(double value) {
  const double angle = 0.5 * PI * value;
  const double cc = cos(angle);
  const double ss = sin(angle);
  Eigen::Matrix2cd matr;
  matr << cc, -i_ * ss, -i_ * ss, cc;
  return matr;
}

Eigen::Matrix2cd GateUnitaryMatrixImplementations::Ry(double value) {
  const double angle = 0.5 * PI * value;
  const double cc = cos(angle);
  const double ss = sin(angle);
  Eigen::Matrix2cd matr;
  matr << cc, -ss, ss, cc;
  return matr;
}

Eigen::Matrix2cd GateUnitaryMatrixImplementations::Rz(double value) {
  const double angle = 0.5 * PI * value;
  const double cc = cos(angle);
  const double ss = sin(angle);
  Eigen::Matrix2cd matr;
  matr << cc - i_ * ss, 0, 0, cc + i_ * ss;
  return matr;
}

Eigen::Matrix2cd GateUnitaryMatrixImplementations::U1(double value) {
  Eigen::Matrix2cd matr;
  matr << 1, 0, 0, std::polar(1.0, PI * value);
  return matr;
}

Eigen::Matrix4cd GateUnitaryMatrixImplementations::ISWAP(double alpha) {
  Eigen::Matrix4cd matr = Eigen::Matrix4cd::Identity();
  const double angle = 0.5 * PI * alpha;
  matr(1, 1) = matr(2, 2) = cos(angle);
  matr(2, 1) = matr(1, 2) = i_ * sin(angle);
  return matr;
}

Eigen::Matrix4cd GateUnitaryMatrixImplementations::XXPhase(double alpha) {
  const double angle = 0.5 * PI * alpha;
  Eigen::Matrix4cd matr = cos(angle) * Eigen::Matrix4cd::Identity();
  matr(3, 0) = matr(2, 1) = matr(1, 2) = matr(0, 3) = i_ * (-sin(angle));
  return matr;
}

Eigen::Matrix4cd GateUnitaryMatrixImplementations::YYPhase(double alpha) {
  auto matr = XXPhase(alpha);
  matr(3, 0) = matr(0, 3) = std::conj(matr(3, 0));
  return matr;
}

Eigen::Matrix4cd GateUnitaryMatrixImplementations::ZZPhase(double alpha) {
  Eigen::Matrix4cd matr = Eigen::Matrix4cd::Zero();
  const auto exp_entry = std::polar(1.0, 0.5 * PI * alpha);
  matr(1, 1) = matr(2, 2) = exp_entry;
  matr(0, 0) = matr(3, 3) = std::conj(exp_entry);
  return matr;
}

Matrix8cd GateUnitaryMatrixImplementations::XXPhase3(double alpha) {
  Eigen::Matrix2cd PauliI, PauliX;
  PauliI = Eigen::Matrix2cd::Identity();
  PauliX << 0, 1, 1, 0;

  Eigen::Matrix4cd XI = Eigen::kroneckerProduct(PauliX, PauliI);
  Eigen::Matrix4cd IX = Eigen::kroneckerProduct(PauliI, PauliX);
  Matrix8cd XXI = Eigen::kroneckerProduct(PauliX, XI);
  Matrix8cd IXX = Eigen::kroneckerProduct(IX, PauliX);
  Matrix8cd XIX = Eigen::kroneckerProduct(XI, PauliX);

  Matrix8cd exponent = -0.5 * alpha * PI * i_ * (XXI + IXX + XIX);
  Matrix8cd matr = exponent.exp();

  return matr;
}

Eigen::Matrix4cd GateUnitaryMatrixImplementations::ESWAP(double alpha) {
  Eigen::Matrix4cd matr = Eigen::Matrix4cd::Identity();
  const double angle = 0.5 * PI * alpha;
  const double cc = cos(angle);
  const double ss = sin(angle);
  matr(0, 0) = matr(3, 3) = std::complex(cc, -ss);
  matr(1, 1) = matr(2, 2) = cc;
  matr(1, 2) = matr(2, 1) = -i_ * ss;
  return matr;
}

Eigen::Matrix4cd GateUnitaryMatrixImplementations::FSim(
    double alpha, double beta) {
  Eigen::Matrix4cd matr = Eigen::Matrix4cd::Identity();
  const double angle = PI * alpha;
  matr(1, 1) = matr(2, 2) = cos(angle);
  matr(1, 2) = matr(2, 1) = -i_ * sin(angle);
  matr(3, 3) = std::polar(1.0, -PI * beta);
  return matr;
}

Eigen::MatrixXcd GateUnitaryMatrixImplementations::PhaseGadget(
    unsigned int number_of_qubits, double alpha) {
  const auto diagonal_entries =
      PhaseGadget_diagonal_entries(number_of_qubits, alpha);
  return diagonal_entries.asDiagonal();
}

// Returns either 0 or 1.
static unsigned number_of_bits_mod_2_in_binary_expansion(unsigned ii) {
  unsigned result = 0;
  while (ii > 0) {
    // Standard trick: if the bits of i are  abc...z10...0,
    // then ~ii = ~(abc...z) 01...1
    unsigned next_bit = (~ii + 1) & ii;
    result = 1 - result;
    ii &= ~next_bit;
  }
  return result;
}

Eigen::VectorXcd GateUnitaryMatrixImplementations::PhaseGadget_diagonal_entries(
    unsigned int number_of_qubits, double alpha) {
  // We will count the bits in the binary expansion and reduce mod 2.
  const unsigned size = get_matrix_size(number_of_qubits);
  Eigen::VectorXcd result(size);

  std::array<std::complex<double>, 2> diagonal_terms;
  diagonal_terms[1] = std::polar(1.0, 0.5 * PI * alpha);
  diagonal_terms[0] = std::conj(diagonal_terms[1]);

  for (unsigned ii = 0; ii < size; ++ii) {
    // Could also do this iteratively/recursively,
    // but maybe not much faster (surely no more than 20 bits, so the
    // number_of_bits function is effectively constant time
    // with a not too large constant)
    result[ii] = diagonal_terms[number_of_bits_mod_2_in_binary_expansion(ii)];
  }
  return result;
}

}  // namespace internal
}  // namespace tket
