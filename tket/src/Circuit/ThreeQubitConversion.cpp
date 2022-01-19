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

#include "ThreeQubitConversion.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>

#include "CircUtils.hpp"
#include "Circuit.hpp"
#include "Circuit/Command.hpp"
#include "Gate/GatePtr.hpp"
#include "Gate/Rotation.hpp"
#include "OpType/OpType.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Constants.hpp"
#include "Utils/CosSinDecomposition.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

// Return a 3-qubit circuit implementing the unitary
//     [ D     ]
//     [    D* ]
// using 4 Rz and 4 CX operations, where D is a 4x4 diagonal unitary matrix.
// The circuit consists of Rz operations on qubit 0 and CX operations with
// target qubit 0.
static Circuit two_qubit_diag_adjoint_plex(const Eigen::Matrix4cd &D) {
  // Compute angles a_i such that D_ii = e^{-i pi/2 a_i}.
  static const double f = -2 / PI;
  double a0 = f * std::arg(D(0, 0));
  double a1 = f * std::arg(D(1, 1));
  double a2 = f * std::arg(D(2, 2));
  double a3 = f * std::arg(D(3, 3));
  double t0 = (a0 + a1 + a2 + a3) / 4;
  double t1 = (a0 + a1 - a2 - a3) / 4;
  double t2 = (a0 - a1 - a2 + a3) / 4;
  double t3 = (a0 - a1 + a2 - a3) / 4;
  Circuit circ(3);
  circ.add_op<unsigned>(OpType::Rz, t0, {0});
  circ.add_op<unsigned>(OpType::CX, {1, 0});
  circ.add_op<unsigned>(OpType::Rz, t1, {0});
  circ.add_op<unsigned>(OpType::CX, {2, 0});
  circ.add_op<unsigned>(OpType::Rz, t2, {0});
  circ.add_op<unsigned>(OpType::CX, {1, 0});
  circ.add_op<unsigned>(OpType::Rz, t3, {0});
  circ.add_op<unsigned>(OpType::CX, {2, 0});
  return circ;
}

// Return a 3-qubit circuit and a unit complex number z which together implement
// the unitary
//     [ U0    ]
//     [    U1 ]
// where U0 and U1 are 4x4 unitaries.
//
// The unitary is implemented by the circuit followed by the diagonal operator
// diag(z, z*, z*, z) on qubits 1 and 2.
//
// If `extract_final_diagonal` is false, then z=1 and the circuit implements the
// unitary exactly, using 9 CX gates.
//
// If `extract_final_diagonal` is true, then the circuit uses 8 CX gates.
static std::pair<Circuit, Complex> two_qubit_plex(
    const Eigen::Matrix4cd &U0, const Eigen::Matrix4cd &U1,
    bool extract_final_diagonal) {
  // 1. Decompose U0 U1* as L T L* where L and T are unitary and T is diagonal.
  Eigen::Matrix4cd U = U0 * U1.adjoint();
  Eigen::ComplexSchur<Eigen::Matrix4cd> schur(U);
  Eigen::Matrix4cd L = schur.matrixU();
  Eigen::Matrix4cd T = schur.matrixT();
  // By construction T is unitary and upper-triangular, hence diagonal.
  TKET_ASSERT(T.isDiagonal());
  // 2. Let D = sqrt(T)
  Eigen::Matrix4cd D = Eigen::Matrix4cd::Zero();
  for (unsigned i = 0; i < 4; i++) {
    D(i, i) = sqrt(T(i, i));
  }
  // 3. Compute R such that U0 = L D R and U1 = L D* R.
  Eigen::Matrix4cd R = D * L.adjoint() * U1;
  // 4. Decompose R into a 2-CX circuit followed by a diagonal.
  auto [R_circ, w0] = decompose_2cx_DV(R);
  Complex w1 = std::conj(w0);
  // 5. Construct the circuit.
  // The diagonal from R's decomposition commutes forward through the controls
  // on qubits 1 and 2, so can be absorbed into L.
  Circuit circ(3);
  unit_map_t qm;
  qm.insert({Qubit(0), Qubit(1)});
  qm.insert({Qubit(1), Qubit(2)});
  circ.append_with_map(R_circ, qm);
  circ.append(two_qubit_diag_adjoint_plex(D));
  L.col(0) *= w0;
  L.col(1) *= w1;
  L.col(2) *= w1;
  L.col(3) *= w0;
  if (extract_final_diagonal) {
    auto [L_circ, z0] = decompose_2cx_DV(L);
    circ.append_with_map(L_circ, qm);
    return {circ, z0};
  } else {
    circ.append_with_map(two_qubit_canonical(L), qm);
    return {circ, 1};
  }
}

// Return a 3-qubit circuit implementing the unitary
//     [ C  -S ]
//     [ DS DC ]
// using H, Ry and CX operations, where C and S are 4x4 real diagonal matrices,
// C^2 + S^2 = I, and D = diag(1,-1,1,-1).
//
// Note that
//     [ C  -S ] = U [ C -S ]
//     [ DS DC ]     [ S  C ]
// where U represents a CZ on qubits 0 and 2. We convert the CZ operations to CX
// by adding Hadamards and simplifying H Ry(t) H to Ry(-t).
static Circuit two_qubit_modified_cossin_circ(
    const Eigen::Matrix4d &C, const Eigen::Matrix4d &S) {
  // Compute angles a_i such that C_ii = cos(pi/2 a_i) and S_ii = sin(pi/2 a_i).
  static const double f = 2 / PI;
  double a0 = f * atan2(S(0, 0), C(0, 0));
  double a1 = f * atan2(S(1, 1), C(1, 1));
  double a2 = f * atan2(S(2, 2), C(2, 2));
  double a3 = f * atan2(S(3, 3), C(3, 3));
  double t0 = (a0 + a1 + a2 + a3) / 4;
  double t1 = (a0 + a1 - a2 - a3) / 4;
  double t2 = (a0 - a1 - a2 + a3) / 4;
  double t3 = (a0 - a1 + a2 - a3) / 4;
  Circuit circ(3);
  circ.add_op<unsigned>(OpType::Ry, t0, {0});
  circ.add_op<unsigned>(OpType::H, {0});
  circ.add_op<unsigned>(OpType::CX, {1, 0});
  circ.add_op<unsigned>(OpType::Ry, -t1, {0});
  circ.add_op<unsigned>(OpType::CX, {2, 0});
  circ.add_op<unsigned>(OpType::Ry, -t2, {0});
  circ.add_op<unsigned>(OpType::CX, {1, 0});
  circ.add_op<unsigned>(OpType::H, {0});
  circ.add_op<unsigned>(OpType::Ry, t3, {0});
  return circ;
}

// Given matrices U and V, if UV* is a scalar multiple of the identity, return
// the scalar.
static std::optional<Complex> id_coeff(
    const Eigen::Matrix4cd &U, const Eigen::Matrix4cd &V) {
  Eigen::Matrix4cd W = U * V.adjoint();
  Complex w = W(0, 0);
  if (W.isApprox(w * Eigen::Matrix4cd::Identity())) {
    return w;
  }
  if (W.isZero()) {
    // Eigen's isApprox compares multiplicatively, so doesn't catch this case.
    return 0.;
  }
  return std::nullopt;
}

// If the given 8x8 unitary represents a circuit with no entanglement between
// qubit 0 and the other qubits, return a decomposition into a circuit on qubit
// 0 and a circuit on qubits 1 and 2. (The qubits on the second circuit are
// indexed with 0 and 1.)
static std::optional<std::pair<Circuit, Circuit>> separate_0_12(
    const Eigen::MatrixXcd &U) {
  // We want to check whether the unitary is of the form
  // [ w_00 V  w_01 V ]
  // [ w_10 V  w_11 V ]
  // where W is a 2x2 unitary and V is a 4x4 unitary.
  // W.l.o.g. we will assume w_00 (or w_10) is real and positive, compute the
  // w_ij assuming the above form, and then check that the form is correct.
  Eigen::Matrix4cd U00 = U.topLeftCorner(4, 4);
  Eigen::Matrix4cd U01 = U.topRightCorner(4, 4);
  Eigen::Matrix4cd U10 = U.bottomLeftCorner(4, 4);
  Eigen::Matrix4cd U11 = U.bottomRightCorner(4, 4);
  // If U is of the desired form, then U_ij U_kl* = w_ij W_kl* I for all
  // i, j, k, l.
  std::optional<Complex> w0000 = id_coeff(U00, U00);  // |w_00|^2
  if (!w0000) return std::nullopt;
  std::optional<Complex> w0101 = id_coeff(U01, U01);  // |w_01|^2
  if (!w0101) return std::nullopt;
  double x0000 = w0000->real(), y0000 = w0000->imag();
  double x0101 = w0101->real(), y0101 = w0101->imag();
  if (std::abs(y0000) > EPS || std::abs(y0101) > EPS || x0000 < -EPS ||
      x0101 < -EPS) {
    return std::nullopt;
  }
  // If negative because of rounding errors, clamp to 0:
  if (x0000 < 0.) x0000 = 0.;
  if (x0101 < 0.) x0101 = 0.;
  // By unitarity of U, w0000 + w0101 = 1, and so x0000 + x0101 = 1.
  // Choose the larger one to work with.
  Complex w00, w01, w10, w11;
  Eigen::Matrix4cd V;
  if (x0000 >= x0101) {
    w00 = sqrt(x0000);
    V = U00 / w00;
    std::optional<Complex> w0001 = id_coeff(U00, U01);
    if (!w0001) return std::nullopt;
    std::optional<Complex> w0010 = id_coeff(U00, U10);
    if (!w0010) return std::nullopt;
    std::optional<Complex> w0011 = id_coeff(U00, U11);
    if (!w0011) return std::nullopt;
    w01 = std::conj(*w0001) / w00;
    w10 = std::conj(*w0010) / w00;
    w11 = std::conj(*w0011) / w00;
  } else {
    w01 = sqrt(x0101);
    V = U01 / w01;
    std::optional<Complex> w0100 = id_coeff(U01, U00);
    if (!w0100) return std::nullopt;
    std::optional<Complex> w0110 = id_coeff(U01, U10);
    if (!w0110) return std::nullopt;
    std::optional<Complex> w0111 = id_coeff(U01, U11);
    if (!w0111) return std::nullopt;
    w00 = std::conj(*w0100) / w01;
    w10 = std::conj(*w0110) / w01;
    w11 = std::conj(*w0111) / w01;
  }

  Eigen::Matrix2cd W;
  W << w00, w01, w10, w11;
  if (U.isApprox(Eigen::KroneckerProduct(W, V))) {
    std::vector<double> angs = tk1_angles_from_unitary(W);
    Circuit c_1q(1);
    c_1q.add_op<unsigned>(OpType::TK1, {angs[0], angs[1], angs[2]}, {0});
    c_1q.add_phase(angs[3]);
    Circuit c_2q = two_qubit_canonical(V);
    return std::pair<Circuit, Circuit>(c_1q, c_2q);
  }

  return std::nullopt;
}

// Special cases worth handling. This is not necessary for correctness, but
// allows us to obtain circuits that are more amenable to later optimization.
static std::optional<Circuit> special_3q_synth(const Eigen::MatrixXcd &U) {
  static const Eigen::PermutationMatrix<8> P1 = []() {
    Eigen::VectorXi V1(8);
    V1 << 0, 1, 4, 5, 2, 3, 6, 7;
    return Eigen::PermutationMatrix<8>(V1);
  }();
  static const Eigen::PermutationMatrix<8> P2 = []() {
    Eigen::VectorXi V2(8);
    V2 << 0, 4, 2, 6, 1, 5, 3, 7;
    return Eigen::PermutationMatrix<8>(V2);
  }();

  // Try separating qubit 0 from qubits 1 and 2:
  std::optional<std::pair<Circuit, Circuit>> c0 = separate_0_12(U);
  if (c0) {
    Circuit c_1q = c0->first;
    Circuit c_2q = c0->second;
    Circuit c(3);
    c.append(c_1q);
    c.append_with_map(c_2q, {{Qubit(0), Qubit(1)}, {Qubit(1), Qubit(2)}});
    return c;
  }

  // Try separating qubit 1 from qubits 0 and 2:
  std::optional<std::pair<Circuit, Circuit>> c1 = separate_0_12(P1 * U * P1);
  if (c1) {
    Circuit c_1q = c1->first;
    Circuit c_2q = c1->second;
    Circuit c(3);
    c.append_with_map(c_1q, {{Qubit(0), Qubit(1)}});
    c.append_with_map(c_2q, {{Qubit(1), Qubit(2)}});
    return c;
  }

  // Try separating qubit 2 from qubits 0 and 1:
  std::optional<std::pair<Circuit, Circuit>> c2 = separate_0_12(P2 * U * P2);
  if (c2) {
    Circuit c_1q = c2->first;
    Circuit c_2q = c2->second;
    Circuit c(3);
    c.append_with_map(c_1q, {{Qubit(0), Qubit(2)}});
    c.append_with_map(c_2q, {{Qubit(0), Qubit(1)}, {Qubit(1), Qubit(0)}});
    return c;
  }

  return std::nullopt;
}

Circuit three_qubit_synthesis(const Eigen::MatrixXcd &U) {
  if (U.rows() != 8 || U.cols() != 8) {
    throw std::invalid_argument("Wrong-size matrix for three-qubit synthesis");
  }

  std::optional<Circuit> c_special = special_3q_synth(U);
  if (c_special) return *c_special;

  auto [l0, l1, r0, r1, c, s] = CS_decomp(U);
  auto [R_circ, z0] = two_qubit_plex(r0, r1, true);
  Complex z1 = std::conj(z0);
  Circuit circ(3);
  circ.append(R_circ);
  circ.append(two_qubit_modified_cossin_circ(c, s));
  // We chopped off the last CZ (on qubits 0 and 2) from the circuit
  // implementing the CS decomposition. Account for this by changing the signs
  // of columns 1 and 3 of l1.
  //
  // We also carried a diagonal from the earlier subcircuit, which commutes
  // through the controls in the middle circuit and merges with l0 and l1.
  //
  // Together these imply the following adjustments:
  l0.col(0) *= z0;
  l0.col(1) *= z1;
  l0.col(2) *= z1;
  l0.col(3) *= z0;
  l1.col(0) *= z0;
  l1.col(1) *= -z1;
  l1.col(2) *= z1;
  l1.col(3) *= -z0;
  circ.append(two_qubit_plex(l0, l1, false).first);
  return circ;
}

Eigen::MatrixXcd get_3q_unitary(const Circuit &c) {
  if (c.n_qubits() != 3) {
    throw CircuitInvalidity("Circuit in get_3q_unitary must have 3 qubits");
  }

  // Construct map from qubits to indices {0,1,2}.
  qubit_vector_t all_qbs = c.all_qubits();
  std::map<Qubit, unsigned> idx;
  for (unsigned i = 0; i < 3; i++) {
    idx[all_qbs[i]] = i;
  }

  // Step through commands, building unitary as we go.
  Eigen::MatrixXcd U = Eigen::MatrixXcd::Identity(8, 8);
  for (const Command &cmd : c) {
    qubit_vector_t qbs = cmd.get_qubits();
    Op_ptr op = cmd.get_op_ptr();
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(8, 8);
    switch (qbs.size()) {
      case 1: {
        std::vector<Expr> angles = as_gate_ptr(op)->get_tk1_angles();
        Eigen::Matrix2cd u = get_matrix_from_tk1_angles(angles);
        // Construct the 8x8 matrix representing u.
        unsigned i = idx[qbs[0]];
        switch (i) {
          case 0:
            M(0, 0) = M(1, 1) = M(2, 2) = M(3, 3) = u(0, 0);
            M(0, 4) = M(1, 5) = M(2, 6) = M(3, 7) = u(0, 1);
            M(4, 0) = M(5, 1) = M(6, 2) = M(7, 3) = u(1, 0);
            M(4, 4) = M(5, 5) = M(6, 6) = M(7, 7) = u(1, 1);
            break;
          case 1:
            M(0, 0) = M(1, 1) = M(4, 4) = M(5, 5) = u(0, 0);
            M(0, 2) = M(1, 3) = M(4, 6) = M(5, 7) = u(0, 1);
            M(2, 0) = M(3, 1) = M(6, 4) = M(7, 5) = u(1, 0);
            M(2, 2) = M(3, 3) = M(6, 6) = M(7, 7) = u(1, 1);
            break;
          case 2:
            M(0, 0) = M(2, 2) = M(4, 4) = M(6, 6) = u(0, 0);
            M(0, 1) = M(2, 3) = M(4, 5) = M(6, 7) = u(0, 1);
            M(1, 0) = M(3, 2) = M(5, 4) = M(7, 6) = u(1, 0);
            M(1, 1) = M(3, 3) = M(5, 5) = M(7, 7) = u(1, 1);
            break;
          default:
            TKET_ASSERT(!"Invalid index");
        }
        break;
      }
      case 2: {
        if (op->get_type() != OpType::CX) {
          throw CircuitInvalidity(
              "Circuit in get_3q_unitary contains non-CX 2-qubit gates");
        }
        unsigned i = idx[qbs[0]];
        unsigned j = idx[qbs[1]];
        // Construct the 8x8 matrix representing CX.
        switch (i) {
          case 0:
            switch (j) {
              case 1:
                M(0, 0) = M(1, 1) = M(2, 2) = M(3, 3) = M(4, 6) = M(5, 7) =
                    M(6, 4) = M(7, 5) = 1;
                break;
              case 2:
                M(0, 0) = M(1, 1) = M(2, 2) = M(3, 3) = M(4, 5) = M(5, 4) =
                    M(6, 7) = M(7, 6) = 1;
                break;
              default:
                TKET_ASSERT(!"Invalid index");
            }
            break;
          case 1:
            switch (j) {
              case 0:
                M(0, 0) = M(1, 1) = M(2, 6) = M(3, 7) = M(4, 4) = M(5, 5) =
                    M(6, 2) = M(7, 3) = 1;
                break;
              case 2:
                M(0, 0) = M(1, 1) = M(2, 3) = M(3, 2) = M(4, 4) = M(5, 5) =
                    M(6, 7) = M(7, 6) = 1;
                break;
              default:
                TKET_ASSERT(!"Invalid index");
            }
            break;
          case 2:
            switch (j) {
              case 0:
                M(0, 0) = M(1, 5) = M(2, 2) = M(3, 7) = M(4, 4) = M(5, 1) =
                    M(6, 6) = M(7, 3) = 1;
                break;
              case 1:
                M(0, 0) = M(1, 3) = M(2, 2) = M(3, 1) = M(4, 4) = M(5, 7) =
                    M(6, 6) = M(7, 5) = 1;
                break;
              default:
                TKET_ASSERT(!"Invalid index");
            }
            break;
          default:
            TKET_ASSERT(!"Invalid index");
        }
        break;
      }
      default:
        throw CircuitInvalidity(
            "Circuit in get_3q_unitary contains gates with more than 2 qubits");
    }
    U = M * U;
  }

  return std::exp(i_ * PI * eval_expr(c.get_phase()).value()) * U;
}

}  // namespace tket
