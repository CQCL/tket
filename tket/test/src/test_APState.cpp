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

#include <catch2/catch_test_macros.hpp>

#include "Simulation/ComparisonFunctions.hpp"
#include "testutil.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Clifford/APState.hpp"
#include "tket/Converters/Converters.hpp"

namespace tket {
namespace test_APState {

void test_apply_gate(
    const MatrixXb& A, const VectorXb& B, const MatrixXb& E,
    const Eigen::VectorXi& P, OpType ot, const std::vector<unsigned>& args) {
  MatrixXb C = MatrixXb::Zero(A.cols(), A.cols());
  APState ap(A, B, C, E, P, 0);
  ap.verify();
  auto sv_before = ap.to_statevector();

  Circuit circ(A.cols());
  circ.add_op<unsigned>(ot, args);
  auto gate_u = tket_sim::get_unitary(circ);

  ap.apply_gate(ot, args);

  ap.verify();
  auto sv_after = ap.to_statevector();

  CHECK(tket_sim::compare_statevectors_or_unitaries(
      gate_u * sv_before, sv_after, tket_sim::MatrixEquivalence::EQUAL));
}

void test_apply_gate_dm(
    const MatrixXb& A, const VectorXb& B, const MatrixXb& C, const MatrixXb& E,
    const Eigen::VectorXi& P, OpType ot, const std::vector<unsigned>& args) {
  APState ap(A, B, C, E, P, 0);
  ap.verify();
  auto dm_before = ap.to_density_matrix();

  Circuit circ(A.cols());
  circ.add_op<unsigned>(ot, args);
  auto gate_u = tket_sim::get_unitary(circ);

  ap.apply_gate(ot, args);

  ap.verify();
  auto dm_after = ap.to_density_matrix();

  CHECK(dm_after.isApprox(gate_u * dm_before * gate_u.adjoint(), EPS));
}

void test_apply_gate_dm_input(
    const MatrixXb& A, const VectorXb& B, const MatrixXb& C, const MatrixXb& E,
    const Eigen::VectorXi& P, OpType ot, const qubit_vector_t& args) {
  ChoiAPState ap(A, B, C, E, P, 0, A.cols());
  ap.ap_.verify();
  auto dm_before = ap.ap_.to_density_matrix();

  Circuit circ(A.cols());
  circ.add_op<Qubit>(ot, args);
  auto gate_u = tket_sim::get_unitary(circ);

  ap.apply_gate(ot, args, ChoiAPState::TableauSegment::Input);

  ap.ap_.verify();
  auto dm_after = ap.ap_.to_density_matrix();

  CHECK(dm_after.isApprox(
      gate_u.transpose() * dm_before * gate_u.conjugate(), EPS));
}

SCENARIO("Normal form") {
  GIVEN("Make A reduced row-echelon form (pure state)") {
    MatrixXb A = MatrixXb::Zero(4, 4);
    VectorXb B = VectorXb::Zero(4);
    MatrixXb C = MatrixXb::Zero(4, 4);
    MatrixXb E = MatrixXb::Zero(4, 4);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
    A(0, 2) = A(0, 3) = true;
    A(1, 0) = A(1, 1) = A(1, 2) = true;
    A(2, 0) = A(2, 1) = A(2, 2) = true;
    A(3, 0) = A(3, 1) = true;
    APState ap(A, B, C, E, P, 0);
    auto sv_before = ap.to_statevector();
    ap.normal_form();
    auto sv_after = ap.to_statevector();
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_before, sv_after, tket_sim::MatrixEquivalence::EQUAL));
    MatrixXb corrA = MatrixXb::Zero(4, 4);
    corrA(0, 0) = corrA(0, 1) = true;
    corrA(1, 2) = true;
    corrA(2, 3) = true;
    CHECK(ap.A == corrA);
  }
  GIVEN("Make A and C reduced row-echelon form") {
    MatrixXb A = MatrixXb::Zero(4, 4);
    VectorXb B = VectorXb::Zero(4);
    MatrixXb C = MatrixXb::Zero(4, 4);
    MatrixXb E = MatrixXb::Zero(4, 4);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
    A(0, 2) = A(0, 3) = true;
    A(1, 0) = A(1, 1) = A(1, 2) = true;
    C(0, 0) = C(0, 1) = C(0, 2) = true;
    C(1, 0) = true;
    C(3, 2) = true;
    APState ap(A, B, C, E, P, 0);
    auto dm_before = ap.to_density_matrix();
    ap.normal_form();
    auto dm_after = ap.to_density_matrix();
    CHECK(dm_after.isApprox(dm_before, EPS));
    MatrixXb corrA = MatrixXb::Zero(4, 4);
    MatrixXb corrC = MatrixXb::Zero(4, 4);
    corrA(0, 0) = corrA(0, 1) = corrA(0, 3) = true;
    corrA(1, 2) = corrA(1, 3) = true;
    corrC(0, 1) = true;
    corrC(1, 3) = true;
    CHECK(ap.A == corrA);
    CHECK(ap.C == corrC);
  }
  GIVEN("Removing leaders from E and P (pure state)") {
    MatrixXb A = MatrixXb::Zero(5, 5);
    VectorXb B = VectorXb::Zero(5);
    MatrixXb C = MatrixXb::Zero(5, 5);
    MatrixXb E = MatrixXb::Zero(5, 5);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(5);
    A(0, 0) = A(0, 2) = A(0, 3) = true;
    A(1, 1) = A(1, 2) = true;
    E(0, 1) = E(1, 0) = true;
    E(0, 3) = E(3, 0) = true;
    E(0, 4) = E(4, 0) = true;
    for (unsigned b = 0; b < 2; ++b) {
      for (unsigned p = 0; p < 4; ++p) {
        B(0) = (b == 1);
        P(0) = p;
        APState ap(A, B, C, E, P, 0);
        auto sv_before = ap.to_statevector();
        ap.normal_form();
        auto sv_after = ap.to_statevector();
        // Just check using statevector; too much changes in each case to nicely
        // test the matrices
        CHECK(tket_sim::compare_statevectors_or_unitaries(
            sv_before, sv_after, tket_sim::MatrixEquivalence::EQUAL));
      }
    }
  }
  GIVEN("Removing mixed qubits from E and P") {
    MatrixXb A = MatrixXb::Zero(5, 5);
    VectorXb B = VectorXb::Zero(5);
    MatrixXb C = MatrixXb::Zero(5, 5);
    MatrixXb E = MatrixXb::Zero(5, 5);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(5);
    C(0, 0) = C(0, 2) = C(0, 3) = true;
    C(1, 1) = true;
    E(0, 1) = E(1, 0) = true;
    E(0, 2) = E(2, 0) = true;
    E(0, 4) = E(4, 0) = true;
    for (unsigned p = 0; p < 4; ++p) {
      P(0) = p;
      APState ap(A, B, C, E, P, 0);
      auto dm_before = ap.to_density_matrix();
      ap.normal_form();
      auto dm_after = ap.to_density_matrix();
      // Just check using statevector; too much changes in each case to nicely
      // test the matrices
      CHECK(dm_after.isApprox(dm_before, EPS));
    }
  }
}

SCENARIO("CZ cases") {
  GIVEN("CZ on free qubits") {
    MatrixXb A = MatrixXb::Zero(4, 4);
    VectorXb B = VectorXb::Zero(4);
    MatrixXb E = MatrixXb::Zero(4, 4);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
    A(0, 0) = A(0, 1) = A(0, 3) = true;
    E(2, 3) = E(3, 2) = true;
    test_apply_gate(A, B, E, P, OpType::CZ, {1, 2});
  }
  GIVEN("CZ on one leading qubit and connected free") {
    for (unsigned b = 0; b < 2; ++b) {
      MatrixXb A = MatrixXb::Zero(5, 5);
      VectorXb B = VectorXb::Zero(5);
      MatrixXb E = MatrixXb::Zero(5, 5);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(5);
      A(0, 0) = A(0, 2) = A(0, 3) = true;
      A(1, 1) = A(1, 3) = A(1, 4) = true;
      B(1) = (b == 1);
      test_apply_gate(A, B, E, P, OpType::CZ, {0, 3});
    }
  }
  GIVEN("CZ on one leading qubit and unconnected free") {
    for (unsigned b = 0; b < 2; ++b) {
      MatrixXb A = MatrixXb::Zero(5, 5);
      VectorXb B = VectorXb::Zero(5);
      MatrixXb E = MatrixXb::Zero(5, 5);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(5);
      A(0, 0) = A(0, 2) = A(0, 3) = true;
      A(1, 1) = A(1, 3) = A(1, 4) = true;
      B(1) = (b == 1);
      test_apply_gate(A, B, E, P, OpType::CZ, {0, 4});
    }
  }
  GIVEN("CZ on leading qubits") {
    for (unsigned b1 = 0; b1 < 2; ++b1) {
      for (unsigned b2 = 0; b2 < 2; ++b2) {
        MatrixXb A = MatrixXb::Zero(8, 8);
        VectorXb B = VectorXb::Zero(8);
        MatrixXb E = MatrixXb::Zero(8, 8);
        Eigen::VectorXi P = Eigen::VectorXi::Zero(8);
        A(0, 0) = A(0, 2) = A(0, 3) = A(0, 4) = A(0, 5) = true;
        A(1, 1) = A(1, 4) = A(1, 5) = A(1, 6) = A(1, 7) = true;
        B(0) = (b1 == 1);
        B(1) = (b2 == 1);
        test_apply_gate(A, B, E, P, OpType::CZ, {0, 1});
      }
    }
  }
  GIVEN("CZ on mixed state") {
    for (unsigned b1 = 0; b1 < 2; ++b1) {
      for (unsigned b2 = 0; b2 < 2; ++b2) {
        MatrixXb A = MatrixXb::Zero(5, 5);
        VectorXb B = VectorXb::Zero(5);
        MatrixXb C = MatrixXb::Zero(5, 5);
        MatrixXb E = MatrixXb::Zero(5, 5);
        Eigen::VectorXi P = Eigen::VectorXi::Zero(5);
        A(0, 1) = A(0, 3) = A(0, 4) = true;
        A(1, 0) = A(1, 2) = A(1, 3) = true;
        B(0) = (b1 == 1);
        B(1) = (b1 == 1);
        C(0, 0) = C(0, 2) = true;
        E(0, 3) = E(3, 0) = true;
        E(1, 3) = E(3, 1) = true;
        P(0) = 3;
        P(1) = 1;
        test_apply_gate_dm(A, B, C, E, P, OpType::CZ, {0, 1});
      }
    }
  }
}

SCENARIO("S cases") {
  GIVEN("S on free qubit") {
    MatrixXb A = MatrixXb::Zero(3, 3);
    VectorXb B = VectorXb::Zero(3);
    MatrixXb E = MatrixXb::Zero(3, 3);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
    A(0, 0) = A(0, 1) = A(0, 2) = true;
    E(1, 2) = E(2, 1) = true;
    test_apply_gate(A, B, E, P, OpType::S, {2});
  }
  GIVEN("S on leading qubit") {
    for (unsigned b = 0; b < 2; ++b) {
      MatrixXb A = MatrixXb::Zero(4, 4);
      VectorXb B = VectorXb::Zero(4);
      MatrixXb E = MatrixXb::Zero(4, 4);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
      A(0, 0) = A(0, 1) = A(0, 2) = true;
      B(0) = (b == 1);
      E(1, 3) = E(3, 1) = true;
      test_apply_gate(A, B, E, P, OpType::S, {0});
    }
  }
  GIVEN("S on disconnected leading qubit") {
    for (unsigned b = 0; b < 2; ++b) {
      MatrixXb A = MatrixXb::Zero(1, 1);
      VectorXb B = VectorXb::Zero(1);
      MatrixXb E = MatrixXb::Zero(1, 1);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(1);
      A(0, 0) = true;
      B(0) = (b == 1);
      test_apply_gate(A, B, E, P, OpType::S, {0});
    }
  }
  GIVEN("S on a mixed state") {
    MatrixXb A = MatrixXb::Zero(3, 3);
    VectorXb B = VectorXb::Zero(3);
    MatrixXb C = MatrixXb::Zero(3, 3);
    MatrixXb E = MatrixXb::Zero(3, 3);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
    A(0, 0) = A(0, 1) = A(0, 2) = true;
    C(0, 1) = C(0, 2) = true;
    E(1, 2) = E(2, 1) = true;
    test_apply_gate_dm(A, B, C, E, P, OpType::S, {2});
  }
}

SCENARIO("V cases") {
  GIVEN("V on leading qubit") {
    for (unsigned b = 0; b < 2; ++b) {
      MatrixXb A = MatrixXb::Zero(4, 4);
      VectorXb B = VectorXb::Zero(4);
      MatrixXb E = MatrixXb::Zero(4, 4);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
      A(0, 0) = A(0, 2) = A(0, 3) = true;
      A(1, 1) = A(1, 3) = true;
      B(0) = (b == 1);
      test_apply_gate(A, B, E, P, OpType::V, {0});
    }
  }
  GIVEN("V on free qubit with some leading") {
    for (unsigned b = 0; b < 2; ++b) {
      for (unsigned p = 0; p < 4; ++p) {
        MatrixXb A = MatrixXb::Zero(9, 9);
        VectorXb B = VectorXb::Zero(9);
        MatrixXb E = MatrixXb::Zero(9, 9);
        Eigen::VectorXi P = Eigen::VectorXi::Zero(9);
        A(0, 0) = A(0, 2) = A(0, 4) = A(0, 5) = A(0, 7) = true;
        A(1, 1) = A(1, 2) = A(1, 3) = A(1, 4) = A(1, 5) = A(1, 6) = true;
        B(1) = (b == 1);
        E(4, 5) = E(5, 4) = true;
        E(4, 6) = E(6, 4) = true;
        E(4, 7) = E(7, 4) = true;
        E(4, 8) = E(8, 4) = true;
        P(4) = p;
        test_apply_gate(A, B, E, P, OpType::V, {4});
      }
    }
  }
  GIVEN("V on free qubit with some earlier connected free") {
    for (unsigned p1 = 0; p1 < 4; ++p1) {
      for (unsigned p2 = 0; p2 < 4; ++p2) {
        MatrixXb A = MatrixXb::Zero(9, 9);
        VectorXb B = VectorXb::Zero(9);
        MatrixXb E = MatrixXb::Zero(9, 9);
        Eigen::VectorXi P = Eigen::VectorXi::Zero(9);
        A(0, 0) = true;
        A(1, 1) = A(1, 4) = A(1, 6) = true;
        A(2, 2) = A(2, 4) = A(2, 7) = true;
        A(3, 3) = A(3, 4) = A(3, 8) = A(3, 8) = true;
        E(4, 5) = E(5, 4) = true;
        E(4, 7) = E(7, 4) = true;
        E(4, 8) = E(8, 4) = true;
        E(5, 6) = E(6, 5) = true;
        E(5, 7) = E(7, 5) = true;
        E(5, 8) = E(8, 5) = true;
        P(4) = p1;
        P(5) = p2;
        test_apply_gate(A, B, E, P, OpType::V, {5});
      }
    }
  }
  GIVEN("V on free qubit with no earlier connected free") {
    for (unsigned p = 0; p < 4; ++p) {
      MatrixXb A = MatrixXb::Zero(4, 4);
      VectorXb B = VectorXb::Zero(4);
      MatrixXb E = MatrixXb::Zero(4, 4);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
      A(0, 0) = A(0, 2) = true;
      E(1, 2) = E(2, 1) = true;
      E(1, 3) = E(3, 1) = true;
      P(1) = p;
      test_apply_gate(A, B, E, P, OpType::V, {1});
    }
  }
  GIVEN("V on disconnected free qubit") {
    for (unsigned p = 0; p < 4; ++p) {
      MatrixXb A = MatrixXb::Zero(1, 1);
      VectorXb B = VectorXb::Zero(1);
      MatrixXb E = MatrixXb::Zero(1, 1);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(1, 1);
      P(0) = p;
      test_apply_gate(A, B, E, P, OpType::V, {0});
    }
  }
  GIVEN("V on a qubit involved in A and C") {
    for (unsigned p = 0; p < 4; ++p) {
      MatrixXb A = MatrixXb::Zero(7, 7);
      VectorXb B = VectorXb::Zero(7);
      MatrixXb C = MatrixXb::Zero(7, 7);
      MatrixXb E = MatrixXb::Zero(7, 7);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(7);
      A(0, 0) = A(0, 2) = true;
      A(1, 0) = A(1, 1) = A(1, 2) = A(1, 3) = A(1, 4) = true;
      C(0, 0) = C(0, 2) = C(0, 3) = true;
      C(1, 0) = C(1, 1) = C(1, 2) = C(1, 5) = C(1, 6) = true;
      E(0, 3) = E(3, 0) = true;
      E(0, 4) = E(4, 0) = true;
      E(0, 5) = E(5, 0) = true;
      E(0, 6) = E(6, 0) = true;
      P(0) = p;
      test_apply_gate_dm(A, B, C, E, P, OpType::V, {0});
    }
  }
  GIVEN("V on a mixed state with zero A") {
    for (unsigned p = 0; p < 4; ++p) {
      MatrixXb A = MatrixXb::Zero(7, 7);
      VectorXb B = VectorXb::Zero(7);
      MatrixXb C = MatrixXb::Zero(7, 7);
      MatrixXb E = MatrixXb::Zero(7, 7);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(7);
      C(0, 0) = C(0, 2) = true;
      C(1, 0) = C(1, 1) = C(1, 2) = C(1, 3) = C(1, 4) = true;
      C(2, 0) = C(2, 2) = C(0, 3) = true;
      C(3, 0) = C(3, 1) = C(3, 2) = C(3, 5) = C(3, 6) = true;
      E(0, 3) = E(3, 0) = true;
      E(0, 4) = E(4, 0) = true;
      E(0, 5) = E(5, 0) = true;
      E(0, 6) = E(6, 0) = true;
      P(0) = p;
      test_apply_gate_dm(A, B, C, E, P, OpType::V, {0});
    }
  }
}

SCENARIO("Qubit Reset") {
  GIVEN("Reset a qubit with a local state") {
    // Qubit 0 is in the |0> state
    MatrixXb A = MatrixXb::Zero(3, 3);
    VectorXb B = VectorXb::Zero(3);
    MatrixXb C = MatrixXb::Zero(3, 3);
    MatrixXb E = MatrixXb::Zero(3, 3);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
    A(0, 0) = true;
    A(1, 1) = A(1, 2) = true;
    B(1) = true;
    P(2) = 1;
    APState correct(A, B, C, E, P, 0);
    // correct is already in normal form
    for (unsigned s = 0; s < 6; ++s) {
      // s = 0,1,2,3: XY basis states
      // s = 4: |0>
      // s = 5: |1>
      A(0, 0) = (s < 4) ? false : true;
      P(0) = (s < 4) ? s : 0;
      B(0) = (s == 5);
      APState ap(A, B, C, E, P, 0);
      ap.apply_gate(OpType::Reset, {0});
      // Check up to global phase
      ap.normal_form();
      ap.phase = correct.phase;
      CHECK(ap == correct);
    }
  }
  GIVEN("Reset one side of a Bell state") {
    MatrixXb A = MatrixXb::Zero(3, 3);
    VectorXb B = VectorXb::Zero(3);
    MatrixXb C = MatrixXb::Zero(3, 3);
    MatrixXb E = MatrixXb::Zero(3, 3);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
    A(0, 0) = A(0, 1) = true;
    APState ap(A, B, C, E, P, 0);
    ap.apply_gate(OpType::Reset, {0});
    ap.normal_form();
    // Qubit 0 ends in |0>
    A(0, 1) = false;
    // Qubit 1 ends in maximally-mixed state
    C(0, 1) = true;
    APState correct(A, B, C, E, P, 0);
    CHECK(ap == correct);
  }
  GIVEN("Reset on a normal form state") {
    MatrixXb A = MatrixXb::Zero(4, 4);
    VectorXb B = VectorXb::Zero(4);
    MatrixXb C = MatrixXb::Zero(4, 4);
    MatrixXb E = MatrixXb::Zero(4, 4);
    Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
    A(0, 0) = A(0, 1) = A(0, 3) = true;
    B(0) = true;
    C(0, 1) = C(0, 2) = true;
    E(1, 3) = E(3, 1) = true;
    E(2, 3) = E(3, 2) = true;
    P(2) = 1;
    P(3) = 2;
    APState ap(A, B, C, E, P, 0);
    WHEN("Apply to Qubit 0") {
      ap.apply_gate(OpType::Reset, {0});
      ap.normal_form();
      // Qubit 0 in state |0>
      A(0, 1) = A(0, 3) = false;
      B(0) = false;
      // A row becomes a C row (combine with other row for gaussian form)
      C(1, 2) = C(1, 3) = true;
      // More gaussian steps
      C(0, 2) = false;
      C(0, 3) = true;
      // LC about C(1, -) to remove P(2)
      E(2, 3) = E(3, 2) = false;
      P(2) = 0;
      P(3) = 1;
      APState correct(A, B, C, E, P, {0});
      // Check equality up to global phase
      ap.phase = correct.phase;
      CHECK(ap == correct);
    }
    WHEN("Apply to Qubit 1") {
      ap.apply_gate(OpType::Reset, {1});
      ap.normal_form();
      // Correct form verified by hand
      MatrixXb A = MatrixXb::Zero(4, 4);
      VectorXb B = VectorXb::Zero(4);
      MatrixXb C = MatrixXb::Zero(4, 4);
      MatrixXb E = MatrixXb::Zero(4, 4);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
      A(0, 1) = true;
      C(0, 0) = C(0, 3) = true;
      C(1, 2) = true;
      E(0, 3) = E(3, 0) = true;
      E(2, 3) = E(3, 2) = true;
      P(3) = 2;
      APState correct(A, B, C, E, P, {0});
      // Check equality up to global phase
      ap.phase = correct.phase;
      CHECK(ap == correct);
    }
    WHEN("Apply to Qubit 2") {
      ap.apply_gate(OpType::Reset, {2});
      ap.normal_form();
      // Correct form verified by hand
      MatrixXb A = MatrixXb::Zero(4, 4);
      VectorXb B = VectorXb::Zero(4);
      MatrixXb C = MatrixXb::Zero(4, 4);
      MatrixXb E = MatrixXb::Zero(4, 4);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
      A(0, 0) = A(0, 1) = A(0, 3) = true;
      A(1, 2) = true;
      B(0) = true;
      C(0, 1) = true;
      C(1, 3) = true;
      APState correct(A, B, C, E, P, {0});
      // Check equality up to global phase
      ap.phase = correct.phase;
      CHECK(ap == correct);
    }
    WHEN("Apply to Qubit 3") {
      ap.apply_gate(OpType::Reset, {3});
      ap.normal_form();
      // Correct form verified by hand
      MatrixXb A = MatrixXb::Zero(4, 4);
      VectorXb B = VectorXb::Zero(4);
      MatrixXb C = MatrixXb::Zero(4, 4);
      MatrixXb E = MatrixXb::Zero(4, 4);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
      A(0, 3) = true;
      C(0, 0) = C(0, 2) = true;
      C(1, 1) = C(1, 2) = true;
      P(2) = 1;
      APState correct(A, B, C, E, P, {0});
      // Check equality up to global phase
      ap.phase = correct.phase;
      CHECK(ap == correct);
    }
  }
}

SCENARIO("Gate encodings") {
  std::list<std::pair<OpType, std::vector<unsigned>>> test_gates = {
      {OpType::Z, {0}},       {OpType::X, {0}},
      {OpType::Y, {0}},       {OpType::S, {0}},
      {OpType::Sdg, {0}},     {OpType::V, {0}},
      {OpType::Vdg, {0}},     {OpType::SX, {0}},
      {OpType::SXdg, {0}},    {OpType::H, {0}},
      {OpType::CX, {0, 1}},   {OpType::CY, {0, 1}},
      {OpType::CZ, {0, 1}},   {OpType::ZZMax, {0, 1}},
      {OpType::ECR, {0, 1}},  {OpType::ISWAPMax, {0, 1}},
      {OpType::SWAP, {0, 1}}, {OpType::BRIDGE, {0, 1, 2}},
      {OpType::noop, {0}},
  };
  GIVEN("Check Z actions") {
    for (const auto& com : test_gates) {
      MatrixXb A = MatrixXb::Identity(3, 3);
      VectorXb B = VectorXb::Zero(3);
      MatrixXb E = MatrixXb::Zero(3, 3);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
      test_apply_gate(A, B, E, P, com.first, com.second);
    }
  }
  GIVEN("Check X actions") {
    for (const auto& com : test_gates) {
      MatrixXb A = MatrixXb::Zero(3, 3);
      VectorXb B = VectorXb::Zero(3);
      MatrixXb E = MatrixXb::Zero(3, 3);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
      test_apply_gate(A, B, E, P, com.first, com.second);
    }
  }
  GIVEN("Check Z actions on inputs of ChoiAPState") {
    for (const auto& com : test_gates) {
      MatrixXb A = MatrixXb::Identity(3, 3);
      VectorXb B = VectorXb::Zero(3);
      MatrixXb C = MatrixXb::Zero(3, 3);
      MatrixXb E = MatrixXb::Zero(3, 3);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
      qubit_vector_t qbs;
      for (unsigned q : com.second) qbs.push_back(Qubit(q));
      test_apply_gate_dm_input(A, B, C, E, P, com.first, qbs);
    }
  }
  GIVEN("Check X actions on inputs of ChoiAPState") {
    for (const auto& com : test_gates) {
      MatrixXb A = MatrixXb::Zero(3, 3);
      VectorXb B = VectorXb::Zero(3);
      MatrixXb C = MatrixXb::Zero(3, 3);
      MatrixXb E = MatrixXb::Zero(3, 3);
      Eigen::VectorXi P = Eigen::VectorXi::Zero(3);
      qubit_vector_t qbs;
      for (unsigned q : com.second) qbs.push_back(Qubit(q));
      test_apply_gate_dm_input(A, B, C, E, P, com.first, qbs);
    }
  }
}

SCENARIO("Loading from a statevector") {
  MatrixXb A = MatrixXb::Zero(4, 4);
  VectorXb B = VectorXb::Zero(4);
  MatrixXb C = MatrixXb::Zero(4, 4);
  MatrixXb E = MatrixXb::Zero(4, 4);
  Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
  A(0, 0) = A(0, 2) = A(0, 3) = true;
  A(1, 1) = A(1, 2) = true;
  B(0) = true;
  E(2, 3) = E(3, 2) = true;
  P(2) = 1;
  P(3) = 2;
  APState ap(A, B, C, E, P, 0);
  auto sv = ap.to_statevector();
  APState reconstructed(sv);
  auto sv2 = reconstructed.to_statevector();
  CHECK(tket_sim::compare_statevectors_or_unitaries(
      sv, sv2, tket_sim::MatrixEquivalence::EQUAL));
}

SCENARIO("Loading from a density matrix") {
  MatrixXb A = MatrixXb::Zero(4, 4);
  VectorXb B = VectorXb::Zero(4);
  MatrixXb C = MatrixXb::Zero(4, 4);
  MatrixXb E = MatrixXb::Zero(4, 4);
  Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
  A(0, 0) = A(0, 2) = A(0, 3) = true;
  C(0, 1) = C(0, 2) = true;
  C(1, 0) = true;
  B(0) = true;
  E(2, 3) = E(3, 2) = true;
  P(0) = 3;
  P(1) = 2;
  P(2) = 1;
  P(3) = 2;
  APState ap(A, B, C, E, P, 0);
  auto dm = ap.to_density_matrix();
  APState reconstructed(dm);
  auto dm2 = reconstructed.to_density_matrix();
  CHECK(dm2.isApprox(dm, EPS));
  // This state is mixed, so only check equality of normal forms up to phase
  ap.normal_form();
  reconstructed.normal_form();
  reconstructed.phase = ap.phase;
  CHECK(ap == reconstructed);
  THEN("Test serialisation") {
    nlohmann::json j_ap = ap;
    APState ap2(0);
    j_ap.get_to(ap2);
    REQUIRE(ap == ap2);
  }
}

SCENARIO("Converting from/to a circuit") {
  GIVEN("A pure circuit in the standard AP form") {
    Circuit circ(4);
    circ.qubit_create_all();
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 0});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {3, 1});
    circ.add_op<unsigned>(OpType::S, {2});
    circ.add_op<unsigned>(OpType::Z, {3});
    APState ap = circuit_to_apstate(circ);
    auto sv_circ = tket_sim::get_statevector(circ);
    auto sv_ap = ap.to_statevector();
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_ap, tket_sim::MatrixEquivalence::EQUAL));
    ap.normal_form();
    sv_ap = ap.to_statevector();
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_ap, tket_sim::MatrixEquivalence::EQUAL));
    Circuit reconstructed = apstate_to_circuit(ap);
    CHECK(circ == reconstructed);
  }
  GIVEN("A generic pure circuit") {
    Circuit circ(4);
    circ.qubit_create_all();
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CY, {1, 3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::ZZMax, {2, 3});
    APState ap = circuit_to_apstate(circ);
    auto sv_circ = tket_sim::get_statevector(circ);
    auto sv_ap = ap.to_statevector();
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_ap, tket_sim::MatrixEquivalence::EQUAL));
    Circuit reconstructed = apstate_to_circuit(ap);
    auto sv_rec = tket_sim::get_statevector(reconstructed);
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_rec, tket_sim::MatrixEquivalence::EQUAL));
  }
  GIVEN("Initialisations, collapses, discards and post-selections") {
    Circuit circ(5);
    circ.qubit_create(Qubit(1));
    circ.qubit_create(Qubit(2));
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::Collapse, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 1});
    circ.add_op<unsigned>(OpType::CX, {4, 2});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.qubit_discard(Qubit(0));
    ChoiAPState ap = circuit_to_choi_apstate(circ);
    ap.post_select(Qubit(3), ChoiAPState::TableauSegment::Output);
    ap.canonical_column_order();
    ap.normal_form();
    // Define correct form from a hand calculation
    MatrixXb A = MatrixXb::Zero(6, 6);
    VectorXb B = VectorXb::Zero(6);
    MatrixXb C = MatrixXb::Zero(6, 6);
    C(0, 0) = C(0, 3) = true;
    C(1, 1) = true;
    MatrixXb E = MatrixXb::Zero(6, 6);
    E(1, 2) = E(2, 1) = true;
    E(1, 4) = E(4, 1) = true;
    E(1, 5) = E(5, 1) = true;
    E(3, 4) = E(4, 3) = true;
    Eigen::VectorXi P = Eigen::VectorXi::Zero(6);
    P(3) = P(4) = 3;
    // Ignore phase by setting them to match
    APState correct(A, B, C, E, P, ap.ap_.phase);
    CHECK(ap.ap_ == correct);
    std::pair<Circuit, qubit_map_t> res_uni =
        choi_apstate_to_unitary_extension_circuit(ap, {Qubit(1)}, {Qubit(0)});
    // Rebuild state by initialising, post-selecting, etc.
    ChoiAPState res_ap = circuit_to_choi_apstate(res_uni.first);
    qubit_map_t perm;
    for (const std::pair<const Qubit, Qubit>& p : res_uni.second)
      perm.insert({p.second, p.first});
    res_ap.rename_qubits(perm, ChoiAPState::TableauSegment::Output);
    // Post-select/initialise
    res_ap.post_select(Qubit(1), ChoiAPState::TableauSegment::Input);
    res_ap.post_select(Qubit(0), ChoiAPState::TableauSegment::Output);
    // Collapsing q[4] in X basis as per circ
    res_ap.apply_gate(
        OpType::H, {Qubit(4)}, ChoiAPState::TableauSegment::Output);
    res_ap.collapse_qubit(Qubit(4), ChoiAPState::TableauSegment::Output);
    res_ap.apply_gate(
        OpType::H, {Qubit(4)}, ChoiAPState::TableauSegment::Output);
    // Discarding q[0] also removes Z row for q[0], so recreate this by
    // XCollapse at input
    res_ap.apply_gate(
        OpType::H, {Qubit(0)}, ChoiAPState::TableauSegment::Input);
    res_ap.collapse_qubit(Qubit(0), ChoiAPState::TableauSegment::Input);
    res_ap.apply_gate(
        OpType::H, {Qubit(0)}, ChoiAPState::TableauSegment::Input);
    res_ap.canonical_column_order();
    res_ap.normal_form();
    // Mixed state, so only guaranteed up to phase
    res_ap.ap_.phase = ap.ap_.phase;
    CHECK(res_ap == ap);
    THEN("Serialize and deserialize") {
      nlohmann::json j_ap = ap;
      ChoiAPState ap2({});
      j_ap.get_to(ap2);
      CHECK(ap == ap2);
    }
    THEN("Check conversion to/from a tableau") {
      ChoiMixTableau tab = choi_apstate_to_cm_tableau(ap);
      ChoiMixTableau tab2 = circuit_to_cm_tableau(circ);
      tab2.post_select(Qubit(3), ChoiMixTableau::TableauSegment::Output);
      tab.canonical_column_order();
      tab.gaussian_form();
      tab2.canonical_column_order();
      tab2.gaussian_form();
      REQUIRE(tab == tab2);
      ChoiAPState ap2 = cm_tableau_to_choi_apstate(tab);
      ap2.canonical_column_order();
      ap2.normal_form();
      // Converting to a tableau drops phase, so ignore this in equivalence
      // check
      ap2.ap_.phase = ap.ap_.phase;
      REQUIRE(ap == ap2);
    }
  }
}

SCENARIO("Converting from/to a tableau") {
  GIVEN("Check pure state up to global phase using circuit") {
    Circuit circ(8);
    circ.qubit_create_all();
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {5});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {6});
    circ.add_op<unsigned>(OpType::H, {7});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {4, 0});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::CX, {6, 0});
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {7, 5});
    circ.add_op<unsigned>(OpType::CZ, {2, 6});
    circ.add_op<unsigned>(OpType::CZ, {4, 6});
    circ.add_op<unsigned>(OpType::CZ, {4, 7});
    circ.add_op<unsigned>(OpType::CZ, {6, 7});
    circ.add_op<unsigned>(OpType::S, {2});
    circ.add_op<unsigned>(OpType::Sdg, {4});
    circ.add_op<unsigned>(OpType::Z, {7});
    ChoiMixTableau cmt = circuit_to_cm_tableau(circ);
    ChoiAPState ap = cm_tableau_to_choi_apstate(cmt);
    auto sv_circ = tket_sim::get_statevector(circ);
    auto sv_ap = ap.ap_.to_statevector();
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_ap, tket_sim::MatrixEquivalence::EQUAL));
    ChoiMixTableau cmt2 = choi_apstate_to_cm_tableau(ap);
    auto [circ2, perm] = cm_tableau_to_exact_circuit(cmt2);
    qubit_map_t inv;
    for (const std::pair<const Qubit, Qubit>& qp : perm)
      inv.insert({qp.second, qp.first});
    circ2.permute_boundary_output(inv);
    auto sv_circ2 = tket_sim::get_statevector(circ2);
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_circ2,
        tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE));
  }
}

}  // namespace test_APState
}  // namespace tket