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
  APState ap(A, B, E, P, 0);
  ap.verify();
  auto sv_before = ap.to_statevector();

  Circuit circ(A.rows());
  circ.add_op<unsigned>(ot, args);
  auto gate_u = tket_sim::get_unitary(circ);

  ap.apply_gate(ot, args);

  ap.verify();
  auto sv_after = ap.to_statevector();

  CHECK(tket_sim::compare_statevectors_or_unitaries(
      gate_u * sv_before, sv_after, tket_sim::MatrixEquivalence::EQUAL));
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
}

SCENARIO("Loading from a statevector") {
  MatrixXb A = MatrixXb::Zero(4, 4);
  VectorXb B = VectorXb::Zero(4);
  MatrixXb E = MatrixXb::Zero(4, 4);
  Eigen::VectorXi P = Eigen::VectorXi::Zero(4);
  A(0, 0) = A(0, 2) = A(0, 3) = true;
  A(1, 1) = A(1, 2) = true;
  B(0) = true;
  E(2, 3) = E(3, 2) = true;
  P(2) = 1;
  P(3) = 2;
  APState ap(A, B, E, P, 0);
  auto sv = ap.to_statevector();
  APState reconstructed(sv);
  auto sv2 = reconstructed.to_statevector();
  CHECK(tket_sim::compare_statevectors_or_unitaries(
      sv, sv2, tket_sim::MatrixEquivalence::EQUAL));
}

SCENARIO("Converting from/to a circuit") {
  GIVEN("A circuit in the standard AP form") {
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
    Circuit reconstructed = apstate_to_circuit(ap);
    CHECK(circ == reconstructed);
  }
  GIVEN("A generic circuit") {
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
}

SCENARIO("Converting from/to a tableau") {
  GIVEN("Check up to global phase using circuit") {
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
    APState ap = tableau_to_apstate(cmt.tab_);
    auto sv_circ = tket_sim::get_statevector(circ);
    auto sv_ap = ap.to_statevector();
    CHECK(tket_sim::compare_statevectors_or_unitaries(
        sv_circ, sv_ap, tket_sim::MatrixEquivalence::EQUAL));
    SymplecticTableau tab2 = apstate_to_tableau(ap);
    ChoiMixTableau cmt2(tab2.xmat, tab2.zmat, tab2.phase);
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
