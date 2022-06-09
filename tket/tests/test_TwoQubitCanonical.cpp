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

#include <catch2/catch_test_macros.hpp>

#include "Circuit/CircUtils.hpp"
#include "Circuit/Command.hpp"
#include "Gate/Rotation.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_TwoQubitCanonical {

static void check_get_information_content(const Eigen::Matrix4cd &U) {
  Eigen::Matrix2cd PauliX, PauliY, PauliZ;
  PauliX << 0, 1, 1, 0;
  PauliY << 0, -i_, i_, 0;
  PauliZ << 1, 0, 0, -1;
  const auto [K1, A, K2] = get_information_content(U);
  const auto [a, b, c] = A;
  const Eigen::Matrix4cd arg = -0.5 * PI * i_ *
                               (a * Eigen::kroneckerProduct(PauliX, PauliX) +
                                b * Eigen::kroneckerProduct(PauliY, PauliY) +
                                c * Eigen::kroneckerProduct(PauliZ, PauliZ));
  Eigen::Matrix4cd res = K1 * arg.exp() * K2;
  REQUIRE(res.isApprox(U));
}

SCENARIO("Testing get_matrix_from_2qb_circ") {
  Circuit c(2);
  GIVEN("A CX") { c.add_op<unsigned>(OpType::CX, {0, 1}); }
  GIVEN("A reverse CX") { c.add_op<unsigned>(OpType::CX, {1, 0}); }
  GIVEN("A Swap") { c.add_op<unsigned>(OpType::SWAP, {0, 1}); }
  GIVEN("A TK1") { c.add_op<unsigned>(OpType::TK1, {0.3, .2, -.6}, {0}); }
  GIVEN("A TK2") { c.add_op<unsigned>(OpType::TK2, {0.3, .2, -.6}, {0, 1}); }
  GIVEN("A reverse TK2") {
    c.add_op<unsigned>(OpType::TK2, {0.3, .2, -.6}, {1, 0});
  }
  GIVEN("A bunch of gates") {
    c.add_op<unsigned>(OpType::TK1, {0.3, .2, -.6}, {0});
    c.add_op<unsigned>(OpType::TK1, {0.3, 2.39, 1.6}, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Vdg, {0});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::Tdg, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
  }

  THEN("The unitaries from tket_sim and get_matrix are identical") {
    auto u1 = tket_sim::get_unitary(c);
    auto u2 = get_matrix_from_2qb_circ(c);
    REQUIRE(u1.isApprox(u2));
  }
}

SCENARIO("Testing two-qubit canonical forms") {
  GIVEN("Decomposing a kronecker product of matrices (0)") {
    Eigen::Matrix2cd testA, testB, resA, resB;
    testA << 1, 0, 0, exp(i_ * 2.4);
    testB << 0, exp(i_ * 3.01), exp(i_ * 0.45), 0;
    Eigen::Matrix4cd U = Eigen::kroneckerProduct(testA, testB);
    std::tie(resA, resB) = kronecker_decomposition(U);
    resA /= resA(0, 0);
    testB /= testB(1, 0);
    resB /= resB(1, 0);
    bool same = true;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        same &= std::abs(testA(i, j) - resA(i, j)) < ERR_EPS;
        same &= std::abs(testB(i, j) - resB(i, j)) < ERR_EPS;
      }
    }
    REQUIRE(same);
  }

  GIVEN("Decomposing a kronecker product of matrices (1)") {
    Eigen::Matrix2cd testA, testB, resA, resB;
    testA = get_matrix_from_tk1_angles({1.984, 4.480, 2.061, 0});
    testB = get_matrix_from_tk1_angles({0.165, 3.645, 1.062, 0});
    Eigen::Matrix4cd U = Eigen::kroneckerProduct(testA, testB);
    std::tie(resA, resB) = kronecker_decomposition(U);
    testA /= testA(0, 0);
    resA /= resA(0, 0);
    testB /= testB(0, 0);
    resB /= resB(0, 0);
    bool same = true;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        same &= std::abs(testA(i, j) - resA(i, j)) < ERR_EPS;
        same &= std::abs(testB(i, j) - resB(i, j)) < ERR_EPS;
      }
    }
    REQUIRE(same);
  }

  GIVEN("Identifying TK1 parameters from a matrix (0)") {
    Eigen::Matrix2cd test = get_matrix_from_tk1_angles({0, 2.061, 3.103, 0});
    std::vector<double> res = tk1_angles_from_unitary(test);
    Eigen::Matrix2cd res_mat =
        get_matrix_from_tk1_angles({res[0], res[1], res[2], res[3]});
    REQUIRE(test.isApprox(res_mat));
  }

  GIVEN("Identifying TK1 parameters from a matrix (1)") {
    Eigen::Matrix2cd test = get_matrix_from_tk1_angles({1., 1.054, 3.612, 0});
    std::vector<double> res = tk1_angles_from_unitary(test);
    Eigen::Matrix2cd res_mat =
        get_matrix_from_tk1_angles({res[0], res[1], res[2], res[3]});
    REQUIRE(test.isApprox(res_mat));
  }

  GIVEN("Performing a KAK decomposition (0)") {
    Eigen::Matrix4cd test;
    test << 1, 0, 0, 0, 0, 0, exp(i_), 0, 0, exp(i_), 0, 0, 0, 0, 0,
        exp(i_ * 2.814);
    check_get_information_content(test);

    Eigen::Matrix4cd CX;
    CX << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
    check_get_information_content(CX);

    for (unsigned i = 0; i < 100; i++) {
      Eigen::Matrix4cd U = random_unitary(4, i);
      check_get_information_content(U);
    }
  }

  GIVEN("Performing a KAK decomposition (1)") {
    Eigen::Matrix2cd PauliX, PauliY, PauliZ;
    PauliX << 0, 1, 1, 0;
    PauliY << 0, -i_, i_, 0;
    PauliZ << 1, 0, 0, -1;
    // some simple 1-CX circuit
    Circuit circ(2);
    circ.add_op<unsigned>(tket::OpType::Rz, -1.4, {0});
    circ.add_op<unsigned>(tket::OpType::Ry, 1., {1});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.8, {0});
    circ.add_op<unsigned>(tket::OpType::CX, {1, 0});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(tket::OpType::Rx, 1.2, {1});
    Eigen::Matrix4cd U = tket_sim::get_unitary(circ);

    const auto [K1, A, K2] = get_information_content(U);
    const auto [a, b, c] = A;
    const Eigen::Matrix4cd arg = -0.5 * PI * i_ *
                                 (a * Eigen::kroneckerProduct(PauliX, PauliX) +
                                  b * Eigen::kroneckerProduct(PauliY, PauliY) +
                                  c * Eigen::kroneckerProduct(PauliZ, PauliZ));
    Eigen::Matrix4cd res = K1 * arg.exp() * K2;
    REQUIRE(res.isApprox(U));
  }

  GIVEN("Performing a KAK decomposition (2)") {
    Eigen::Matrix2cd PauliX, PauliY, PauliZ;
    PauliX << 0, 1, 1, 0;
    PauliY << 0, -i_, i_, 0;
    PauliZ << 1, 0, 0, -1;
    // some simple 2-CX circuit
    Circuit circ(2);
    circ.add_op<unsigned>(tket::OpType::Rz, -1.4, {0});
    circ.add_op<unsigned>(tket::OpType::Ry, 1., {1});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.8, {0});
    circ.add_op<unsigned>(tket::OpType::CX, {1, 0});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(tket::OpType::Rx, -0.58, {1});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(tket::OpType::CX, {1, 0});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.2, {0});
    Eigen::Matrix4cd U = tket_sim::get_unitary(circ);

    const auto [K1, A, K2] = get_information_content(U);
    const auto [a, b, c] = A;
    const Eigen::Matrix4cd arg = -0.5 * PI * i_ *
                                 (a * Eigen::kroneckerProduct(PauliX, PauliX) +
                                  b * Eigen::kroneckerProduct(PauliY, PauliY) +
                                  c * Eigen::kroneckerProduct(PauliZ, PauliZ));
    Eigen::Matrix4cd res = K1 * arg.exp() * K2;
    REQUIRE(res.isApprox(U));
  }

  GIVEN("Performing a KAK decomposition (3)") {
    Eigen::Matrix2cd PauliX, PauliY, PauliZ;
    PauliX << 0, 1, 1, 0;
    PauliY << 0, -i_, i_, 0;
    PauliZ << 1, 0, 0, -1;
    // Arbitrary matrix
    Eigen::Matrix4cd B;
    B << 1. + 2. * i_, 2. + 3. * i_, 3. + 4. * i_, 4. + 5. * i_, 5. + 6. * i_,
        6. + 7. * i_, 7. + 8. * i_, 8. + 9. * i_, 9. + 1. * i_, 1. + 2. * i_,
        2. + 3. * i_, 3. + 4. * i_, 4. + 5. * i_, 5. + 6. * i_, 6. + 7. * i_,
        7. + 8. * i_;
    Eigen::Matrix4cd A = B + B.adjoint();  // hermitian
    Eigen::Matrix4cd I = Eigen::Matrix4cd::Identity();
    Eigen::Matrix4cd U = (I - i_ * A).inverse() * (I + i_ * A);  // unitary

    const auto [K1, AA, K2] = get_information_content(U);
    const auto [a, b, c] = AA;
    const Eigen::Matrix4cd arg = -0.5 * PI * i_ *
                                 (a * Eigen::kroneckerProduct(PauliX, PauliX) +
                                  b * Eigen::kroneckerProduct(PauliY, PauliY) +
                                  c * Eigen::kroneckerProduct(PauliZ, PauliZ));
    Eigen::Matrix4cd res = K1 * arg.exp() * K2;
    REQUIRE(res.isApprox(U));
  }

  GIVEN(
      "Decomposing information content of a fixed matrix "
      "deterministically") {
    Eigen::Matrix4cd X;
    double s = 1. / sqrt(2.);
    X << s + 0. * i_, 0. + 0. * i_, 0. + 0. * i_, s + -4.32978e-17 * i_,
        0. + 0. * i_, 0.5 + 0.5 * i_, -0.5 + -0.5 * i_, 0. + 0. * i_,
        0.5 + -0.5 * i_, 0. + 0. * i_, 0. + 0. * i_, -0.5 + 0.5 * i_,
        0. + 0. * i_, s + 0. * i_, s + -5.55112e-17 * i_, 0. + 0. * i_;
    auto [K1, A, K2] = get_information_content(X);
    bool all_deterministic = true;
    for (unsigned i = 0; i < 10; ++i) {
      auto [K1p, Ap, K2p] = get_information_content(X);
      if (!(matrices_are_equal(K1p, K1) && Ap == A &&
            matrices_are_equal(K2p, K2))) {
        all_deterministic = false;
        break;
      }
    }
    REQUIRE(all_deterministic);
  }

  GIVEN("Identifying a canonical circuit from a matrix (0)") {
    Eigen::Matrix4cd test;
    test << 1, 0, 0, 0, 0, 0, exp(i_), 0, 0, exp(i_), 0, 0, 0, 0, 0,
        exp(i_ * 2.814);
    Circuit result = two_qubit_canonical(test);
    Eigen::Matrix4cd res = tket_sim::get_unitary(result);
    REQUIRE(res.isApprox(test));
  }

  GIVEN("Identifying a canonical circuit from a matrix (1)") {
    Eigen::Matrix4cd test;
    test << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
    Circuit result = two_qubit_canonical(test);
    Eigen::Matrix4cd res = tket_sim::get_unitary(result);
    REQUIRE(res.isApprox(test));
  }

  GIVEN("Identifying a canonical circuit from a matrix (2)") {
    // Arbitrary matrix
    Eigen::Matrix4cd B;
    B << 1. + 2. * i_, 2. + 3. * i_, 3. + 4. * i_, 4. + 5. * i_, 5. + 6. * i_,
        6. + 7. * i_, 7. + 8. * i_, 8. + 9. * i_, 9. + 1. * i_, 1. + 2. * i_,
        2. + 3. * i_, 3. + 4. * i_, 4. + 5. * i_, 5. + 6. * i_, 6. + 7. * i_,
        7. + 8. * i_;
    Eigen::Matrix4cd A = B + B.adjoint();  // hermitian
    Eigen::Matrix4cd I = Eigen::Matrix4cd::Identity();
    Eigen::Matrix4cd U = (I - i_ * A).inverse() * (I + i_ * A);  // unitary
    Circuit result = two_qubit_canonical(U);
    Eigen::Matrix4cd res = tket_sim::get_unitary(result);
    REQUIRE(res.isApprox(U));
  }

  GIVEN("Identifying a canonical circuit from a matrix (3)") {
    Eigen::Matrix4cd test;
    test << -i_, 1, -i_, 1, -1, i_, 1, -i_, 1, -i_, 1, -i_, -i_, 1, i_, -1;
    test *= 0.5 * exp(i_ * PI * 0.25);
    Circuit result = two_qubit_canonical(test);
    Eigen::Matrix4cd res = tket_sim::get_unitary(result);
    REQUIRE(res.isApprox(test));
  }

  GIVEN("A two qubit circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Eigen::Matrix4cd mat = tket_sim::get_unitary(circ);
    bool success = Transforms::two_qubit_squash().apply(circ);
    REQUIRE(success);
    REQUIRE(circ.count_gates(OpType::CX) == 2);
    Eigen::Matrix4cd result = tket_sim::get_unitary(circ);
    mat /= mat(0, 0);
    result /= result(0, 0);
    bool same = true;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        same &= std::abs(mat(i, j) - result(i, j)) < ERR_EPS;
      }
    }
    REQUIRE(same);
  }

  GIVEN("A two qubit circuit with 0 CNOTs") {
    Circuit circ(2);
    circ.add_op<unsigned>(tket::OpType::Rz, -1.4, {0});
    circ.add_op<unsigned>(tket::OpType::Ry, 1., {1});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.8, {0});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(tket::OpType::Rx, 1.5, {0});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.2, {0});
    circ.add_op<unsigned>(tket::OpType::CX, {0, 1});
    circ.add_op<unsigned>(tket::OpType::CX, {0, 1});
    REQUIRE(Transforms::two_qubit_squash().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
  }

  GIVEN("A two qubit circuit that simplifies") {
    Circuit circ(2);
    circ.add_op<unsigned>(tket::OpType::Rz, -1.4, {0});
    circ.add_op<unsigned>(tket::OpType::Ry, 1., {1});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.8, {0});
    circ.add_op<unsigned>(tket::OpType::CX, {1, 0});
    circ.add_op<unsigned>(tket::OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(tket::OpType::Rx, 1.5, {0});
    circ.add_op<unsigned>(tket::OpType::CX, {0, 1});
    circ.add_op<unsigned>(tket::OpType::Rz, 1.2, {0});
    REQUIRE(Transforms::two_qubit_squash().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 1);
  }

  GIVEN("A swap cannot be simplified") {
    Circuit circ(2);
    add_2qb_gates(circ, OpType::CX, {{1, 0}, {0, 1}, {1, 0}});
    REQUIRE(!Transforms::two_qubit_squash().apply(circ));
  }

  GIVEN("A two qubit circuit with measures") {
    Circuit circ(4);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {0, 1}, {1, 0}});
    circ.add_op<unsigned>(OpType::Collapse, {0});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {0, 1}, {1, 0}});
    circ.add_op<unsigned>(OpType::Collapse, {0});
    circ.add_op<unsigned>(OpType::Collapse, {1});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {3, 2}});
    circ.add_op<unsigned>(OpType::Collapse, {2});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {3, 2}, {2, 3}, {3, 2}});
    REQUIRE(Transforms::two_qubit_squash().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 8);
  }

  GIVEN("An optimal circuit") {
    Circuit circ(3);
    add_2qb_gates(
        circ, OpType::CX,
        {{0, 1}, {1, 0}, {0, 1}, {0, 2}, {0, 1}, {1, 0}, {0, 1}});
    bool success = Transforms::two_qubit_squash().apply(circ);
    REQUIRE(!success);
  }

  GIVEN("Multiple subcircuits to optimise") {
    Circuit circ(4);
    add_2qb_gates(
        circ, OpType::CX,
        {{0, 1},
         {1, 0},
         {0, 1},
         {1, 0},
         {2, 3},
         {3, 2},
         {2, 3},
         {3, 2},
         {0, 2},
         {2, 0},
         {0, 2},
         {2, 0},
         {1, 3},
         {3, 1},
         {1, 3},
         {3, 1}});
    const StateVector s0 = tket_sim::get_statevector(circ);
    bool success = Transforms::two_qubit_squash().apply(circ);
    REQUIRE(success);
    REQUIRE(circ.count_gates(OpType::CX) == 8);
    const StateVector s1 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
}

SCENARIO("Testing two qubit decomposition with fidelity tradeoff") {
  GIVEN("Average fidelity must be greater than 3-CX perfect circuit") {
    // Arbitrary matrix
    Eigen::Matrix4cd B;
    B << 1. + 2. * i_, 2. + 3. * i_, 3. + 4. * i_, 4. + 5. * i_, 5. + 6. * i_,
        6. + 7. * i_, 7. + 8. * i_, 8. + 9. * i_, 9. + 1. * i_, 1. + 2. * i_,
        2. + 3. * i_, 3. + 4. * i_, 4. + 5. * i_, 5. + 6. * i_, 6. + 7. * i_,
        7. + 8. * i_;
    Eigen::Matrix4cd A = B + B.adjoint();  // hermitian
    Eigen::Matrix4cd I = Eigen::Matrix4cd::Identity();
    Eigen::Matrix4cd U = (I - i_ * A).inverse() * (I + i_ * A);  // unitary
    auto get_fid = [&U](const Eigen::Matrix4cd &Up) {
      return (4. + pow(abs((Up.adjoint() * U).trace()), 2)) / 20.;
    };
    bool same = true;
    Circuit circ_out = two_qubit_canonical(U);
    Transforms::TwoQbFidelities fid;
    for (double gate_fid = 0.; gate_fid < 1.; gate_fid += 0.01) {
      Circuit circ_approx = circ_out;
      fid.CX_fidelity = gate_fid;
      decompose_TK2(fid).apply(circ_out);
      Eigen::Matrix4cd out = tket_sim::get_unitary(circ_out);
      const int nb_cx = circ_out.count_gates(OpType::CX);
      const double fid_eff = get_fid(out) * pow(gate_fid, nb_cx);
      const double fid_theo = pow(gate_fid, 3);
      same &= fid_eff > fid_theo - ERR_EPS;
    }
    REQUIRE(same);
  }
}

SCENARIO("KAK Decomposition around symbolic gates") {
  GIVEN("Inefficient two-qubit circuit with symbolic gates") {
    Circuit circ(4);
    Sym a = SymEngine::symbol("alpha");
    Sym b = SymEngine::symbol("beta");
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {0, 1}, {1, 0}});
    circ.add_op<unsigned>(OpType::Rz, {Expr(a)}, {0});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {0, 1}, {1, 0}});
    circ.add_op<unsigned>(OpType::Rx, {-Expr(a)}, {0});
    circ.add_op<unsigned>(OpType::Ry, {Expr(b)}, {1});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {3, 2}});
    circ.add_op<unsigned>(OpType::U2, {0.5, -Expr(b)}, {2});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {3, 2}, {2, 3}, {3, 2}});
    REQUIRE(Transforms::two_qubit_squash().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 8);
  }
  GIVEN("Efficient two-qubit circuit with symbolic gates") {
    Circuit circ(4);
    Sym a = SymEngine::symbol("alpha");
    Sym b = SymEngine::symbol("beta");
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {0, 1}});
    circ.add_op<unsigned>(OpType::Rz, {Expr(a)}, {0});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}});
    circ.add_op<unsigned>(OpType::Rx, {-Expr(a)}, {0});
    circ.add_op<unsigned>(OpType::Ry, {Expr(b)}, {1});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::U2, {0.5, -Expr(b)}, {2});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {3, 2}, {2, 3}});
    REQUIRE_FALSE(Transforms::two_qubit_squash().apply(circ));
  }
}

SCENARIO("Test qubit reversal") {
  GIVEN("A 4x4 matrix") {
    Eigen::Matrix4cd test, correct;
    // clang-format off
        test << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;
        correct << 1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0,
            0, 1, 0, 0;
    // clang-format on
    REQUIRE(matrices_are_equal(reverse_indexing(test), correct));
  }

  GIVEN("An 8x8 matrix") {
    Eigen::Matrix<Complex, 8, 8> test, correct;
    // clang-format off
        test << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
        correct << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 1, 0, 0, 0, 0;
    // clang-format on
    Eigen::MatrixXcd Xtest = test;
    Eigen::MatrixXcd Xcorrect = correct;
    REQUIRE(matrices_are_equal(reverse_indexing(Xtest), Xcorrect));
  }

  GIVEN("An 8-vector") {
    Eigen::Matrix<Complex, 8, 1> test, correct;
    test << 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7;
    correct << 0, 0.4, 0.2, 0.6, 0.1, 0.5, 0.3, 0.7;
    Eigen::VectorXcd Xtest = test;
    Eigen::VectorXcd Xcorrect = correct;
    REQUIRE(matrices_are_equal(reverse_indexing(Xtest), Xcorrect));
  }
}

static void check_decompose_2cx_VD(const Eigen::Matrix4cd &U) {
  auto [circ, z0] = decompose_2cx_VD(U);
  unsigned n_cx = 0;
  for (const Command &cmd : circ) {
    OpType optype = cmd.get_op_ptr()->get_type();
    if (optype == OpType::CX) {
      n_cx++;
    } else {
      CHECK(optype == OpType::TK1);
    }
  }
  CHECK(n_cx <= 2);
  Eigen::Matrix4cd V = tket_sim::get_unitary(circ);
  Eigen::Matrix4cd D = Eigen::Matrix4cd::Zero();
  Complex z1 = std::conj(z0);
  D(0, 0) = z0;
  D(1, 1) = z1;
  D(2, 2) = z1;
  D(3, 3) = z0;
  CHECK(is_unitary(D));
  CHECK(U.isApprox(V * D));
}

static void check_decompose_2cx_DV(const Eigen::Matrix4cd &U) {
  auto [circ, z0] = decompose_2cx_DV(U);
  unsigned n_cx = 0;
  for (const Command &cmd : circ) {
    OpType optype = cmd.get_op_ptr()->get_type();
    if (optype == OpType::CX) {
      n_cx++;
    } else {
      CHECK(optype == OpType::TK1);
    }
  }
  CHECK(n_cx <= 2);
  Eigen::Matrix4cd V = tket_sim::get_unitary(circ);
  Eigen::Matrix4cd D = Eigen::Matrix4cd::Zero();
  Complex z1 = std::conj(z0);
  D(0, 0) = z0;
  D(1, 1) = z1;
  D(2, 2) = z1;
  D(3, 3) = z0;
  CHECK(is_unitary(D));
  CHECK(U.isApprox(D * V));
}

static void check_decompose_2cx_plus_diag(const Eigen::Matrix4cd &U) {
  check_decompose_2cx_VD(U);
  check_decompose_2cx_DV(U);
}

SCENARIO("Test decomposition into 2-CX circuit plus diagonal") {
  GIVEN("A fixed 4x4 unitary") {
    // Randomly generated with scipy.stats.unitary_group.rvs.
    Eigen::Matrix4cd U = Eigen::Matrix4cd::Zero();
    U(0, 0) = {-0.20152561587695295, 0.6507745766671906};
    U(0, 1) = {-0.4408881481052427, 0.27850972852126277};
    U(0, 2) = {0.35512207181773037, -0.27983369659344315};
    U(0, 3) = {0.23006105131436833, 0.08113678275144227};
    U(1, 0) = {0.5137659960929305, -0.039374703160842156};
    U(1, 1) = {-0.7012946739198794, 0.050511013385731204};
    U(1, 2) = {-0.14084755836866267, 0.40342398818925584};
    U(1, 3) = {-0.1880781494682805, 0.14888321804568522};
    U(2, 0) = {0.2840858425126659, -0.33809784885176974};
    U(2, 1) = {-0.15515861149283824, -0.3885892561931721};
    U(2, 2) = {0.1045319779935326, -0.48351730194381587};
    U(2, 3) = {0.49837718713122997, 0.36988314043954695};
    U(3, 0) = {-0.24596349093976072, 0.12190590768740035};
    U(3, 1) = {0.0912551074951825, 0.224234454187113};
    U(3, 2) = {-0.6068434390886989, -0.004194299289027856};
    U(3, 3) = {0.026106715046833248, 0.7050349022743666};
    check_decompose_2cx_plus_diag(U);
  }
  GIVEN("Random unitaries") {
    for (unsigned i = 0; i < 100; i++) {
      Eigen::Matrix4cd U = random_unitary(4, i);
      check_decompose_2cx_plus_diag(U);
    }
  }
  GIVEN("Some special matrices") {
    Eigen::Matrix4cd U;
    // clang-format off
    U << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;
    check_decompose_2cx_plus_diag(U);
    U << 0, 0, 1, 0,
         0, 0, 0, 1,
         1, 0, 0, 0,
         0, 1, 0, 0;
    check_decompose_2cx_plus_diag(U);
    U << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 0, 1,
         0, 0, 1, 0;
    check_decompose_2cx_plus_diag(U);
    U << 0, 0, 0, 1,
         0, 0, 1, 0,
         0, 1, 0, 0,
         1, 0, 0, 0;
    check_decompose_2cx_plus_diag(U);
    // clang-format on
  }
  GIVEN("unitary close to identity") {
    for (unsigned i = 0; i < 20; i++) {
      std::srand(i);
      Eigen::Matrix4d A = Eigen::Matrix4d::Zero();
      for (unsigned r = 0; r < 4; r++) {
        for (unsigned s = 0; s < 4; s++) {
          // If 0.01 is replaced with 0.001 we do get failures. See the
          // commentary for the function `decompose_VD`.
          A(r, s) = 0.01 * rand() / RAND_MAX;
        }
      }
      Eigen::Matrix4cd iH = i_ * (A + A.transpose());
      Eigen::Matrix4cd U = (iH).exp();
      check_decompose_2cx_plus_diag(U);
    }
  }
}

}  // namespace test_TwoQubitCanonical
}  // namespace tket
