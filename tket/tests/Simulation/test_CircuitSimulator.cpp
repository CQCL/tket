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

// This file is specifically for testing of tket-sim functions
// (which, of course, are used in other tests for correctness).
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../Gate/GatesData.hpp"
#include "../testutil.hpp"
#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "ComparisonFunctions.hpp"
#include "Gate/GateUnitaryMatrix.hpp"
#include "Gate/GateUnitaryMatrixUtils.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {
namespace test_CircuitSimulator {
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

SCENARIO("Simple circuits produce the correct statevectors") {
  GIVEN("A 1 qubit circ with X-gate") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::X, {0});
    const StateVector statevector = tket_sim::get_statevector(circ);
    CHECK(statevector(0) == 0.);
    CHECK(statevector(1) == 1.);
  }
  GIVEN("An 8 qubit circ with 8 X-gates") {
    Circuit circ(8);
    for (unsigned i = 0; i < circ.n_qubits(); ++i) {
      circ.add_op<unsigned>(OpType::X, {i});
    }
    const StateVector sv = tket_sim::get_statevector(circ);
    for (unsigned i = 0; i < sv.size() - 1; ++i) {
      CHECK(sv(i) == 0.);  // All elements are 0 except the |111....1> state
    }
    REQUIRE(sv(sv.size() - 1) == 1.);
  }
  GIVEN(
      "An N-qubit circuit with N Z-gates and then N more Z-gates (ie "
      "identity)") {
    unsigned N = 11;
    Circuit circ(N);
    for (unsigned i = 0; i < circ.n_qubits(); ++i) {
      circ.add_op<unsigned>(OpType::Z, {i});
      circ.add_op<unsigned>(OpType::Z, {i});
    }
    const StateVector sv = tket_sim::get_statevector(circ);
    REQUIRE(sv(0) == 1.);
    for (unsigned i = 1; i < sv.size(); ++i) {
      CHECK(sv(i) == 0.);
    }
  }

  GIVEN("A 2-qubit circuit with one hadamard") {
    unsigned N = 2;
    Circuit circ(N);
    circ.add_op<unsigned>(OpType::H, {1});
    WHEN("Statevector is calculated") {
      StateVector sv = tket_sim::get_statevector(circ);

      THEN("Statevector matches") {
        CHECK(sv(0).real() == Approx(0.70710678));
        CHECK(sv(0).imag() == 0.0);
        CHECK(sv(1).real() == Approx(0.70710678));
        CHECK(sv(2) == 0.0);
        CHECK(sv(3) == 0.0);
      }
    }
  }

  GIVEN("A circuit with all statevector elements non-zero.") {
    unsigned N = 3;
    Circuit circ(N);
    circ.add_op<unsigned>(OpType::Rx, 0.1, {0});
    circ.add_op<unsigned>(OpType::Ry, 2.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});

    circ.add_op<unsigned>(OpType::Rx, 1.22, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.4, {1});
    circ.add_op<unsigned>(OpType::Rz, 12, {2});

    WHEN("StateVector is calculated") {
      const StateVector sv = tket_sim::get_statevector(circ);
      THEN("Values match expected.") {
        CHECK(sv(0).real() == Approx(0.43583703));
        CHECK(sv(2).imag() == Approx(0.14799435));
        CHECK(sv(1).real() == Approx(-0.01194328));
        CHECK(sv(3).imag() == Approx(0.03517244));
      }
    }
    WHEN("Statevectors compared before and after some minor changes") {
      const StateVector sv1 = tket_sim::get_statevector(circ);
      circ.add_op<unsigned>(OpType::H, {2});
      circ.add_op<unsigned>(OpType::H, {2});
      circ.add_op<unsigned>(OpType::X, {1});
      circ.add_op<unsigned>(OpType::X, {1});
      circ.add_op<unsigned>(OpType::Rz, 4. / 3, {0});
      circ.add_op<unsigned>(OpType::Rz, 4. / 3, {0});
      circ.add_op<unsigned>(OpType::Rz, 4. / 3, {0});

      const StateVector sv2 = tket_sim::get_statevector(circ);

      REQUIRE(tket_sim::compare_statevectors_or_unitaries(sv1, sv2));

      circ.add_op<unsigned>(OpType::H, {2});
      StateVector sv3 = tket_sim::get_statevector(circ);
      REQUIRE(!tket_sim::compare_statevectors_or_unitaries(sv1, sv3));
    }
  }

  GIVEN("A circuit where the only difference is a swapped round CZ") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});

    Circuit circ2(2);
    circ2.add_op<unsigned>(OpType::H, {0});
    circ2.add_op<unsigned>(OpType::CZ, {1, 0});
    circ2.add_op<unsigned>(OpType::H, {0});
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        tket_sim::get_statevector(circ), tket_sim::get_statevector(circ2)));
  }
}

SCENARIO("Ignored op types don't affect get unitary") {
  Circuit circ1(3);
  // circ2 will add the same ops as circ1, but with extra ops
  Circuit circ2(3, 2);

  circ1.add_op<unsigned>(OpType::H, {0});
  circ2.add_op<unsigned>(OpType::Measure, {0, 1});
  circ2.add_op<unsigned>(OpType::H, {0});

  circ1.add_op<unsigned>(OpType::CZ, {0, 1});
  circ2.add_op<unsigned>(OpType::CZ, {0, 1});

  circ2.add_barrier({1, 2});

  circ1.add_op<unsigned>(OpType::Ry, 2.1, {2});
  circ2.add_op<unsigned>(OpType::noop, {2});
  circ2.add_op<unsigned>(OpType::Ry, 2.1, {2});

  circ2.add_op<unsigned>(OpType::Measure, {2, 0});
  REQUIRE(matrices_are_equal(
      tket_sim::get_statevector(circ1), tket_sim::get_statevector(circ2)));
  REQUIRE(matrices_are_equal(
      tket_sim::get_unitary(circ1), tket_sim::get_unitary(circ2)));
}

SCENARIO("Circuits without gates") {
  // Include the zero qubit case!
  for (unsigned nn = 0; nn < 4; ++nn) {
    const Circuit circ(nn);
    const auto size = get_matrix_size(nn);
    const auto u = tket_sim::get_unitary(circ);
    CHECK(matrices_are_equal(Eigen::MatrixXcd::Identity(size, size), u));

    const auto sv = tket_sim::get_statevector(circ);
    CHECK(matrices_are_equal(Eigen::MatrixXcd::Identity(size, 1), sv));

    // A "random" matrix...
    Eigen::MatrixXcd rectangular_matrix(size, 3);
    auto entry = std::polar(1.0, 0.12345);
    for (unsigned jj = 0; jj < rectangular_matrix.cols(); ++jj) {
      for (unsigned ii = 0; ii < size; ++ii) {
        rectangular_matrix(ii, jj) = entry;
        entry *= std::polar(1.0123456, 0.9876543);
      }
    }
    const auto original_copy = rectangular_matrix;
    tket_sim::apply_unitary(circ, rectangular_matrix);

    // Multiplication by 1.0 or 0.0 should give the EXACT answer,
    // so it's OK to demand exact equality here.
    CHECK(matrices_are_equal(original_copy, rectangular_matrix));
  }
}

SCENARIO("Directly simulate circuits with >= 3 qubit gates") {
  Circuit circ(4);
  circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
  circ.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
  circ.add_op<unsigned>(OpType::CSWAP, {0, 1, 2});
  circ.add_op<unsigned>(OpType::CnRy, 0.1234, {0, 1, 2, 3});
  circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
  circ.add_op<unsigned>(OpType::PhaseGadget, 0.1, {0, 1, 2, 3});
  const auto u = tket_sim::get_unitary(circ);
  REQUIRE(is_unitary(u));
}

SCENARIO("Directly simulate circuits with CircBox") {
  Circuit w(3);
  w.add_op<unsigned>(OpType::Rx, 0.5, {0});
  w.add_op<unsigned>(OpType::CX, {0, 1});
  auto w_copy = w;
  {
    Circuit temp_circ(2);
    temp_circ.add_op<unsigned>(OpType::Ry, 0.75, {0});
    temp_circ.add_op<unsigned>(OpType::CX, {1, 0});
    const CircBox temp_box(temp_circ);
    w.add_box(temp_box, {2, 0});
  }
  // Add some junk at the end.
  w.add_op<unsigned>(OpType::X, {1});
  w.add_op<unsigned>(OpType::X, {2});
  const auto unitary = tket_sim::get_unitary(w);

  // "manually" recreate and add the circ box,
  // remembering the altered qubit indices!
  w_copy.add_op<unsigned>(OpType::Ry, 0.75, {2});
  w_copy.add_op<unsigned>(OpType::CX, {0, 2});
  // Add the same junk at the end.
  w_copy.add_op<unsigned>(OpType::X, {1});
  w_copy.add_op<unsigned>(OpType::X, {2});
  const auto recreated_unitary = tket_sim::get_unitary(w_copy);
  REQUIRE(unitary.isApprox(recreated_unitary));
}

// Whenever a specific test involving circuit simulation fails
// (e.g., proptests), copy it here to verify that it is fixed.
SCENARIO("Specific previous failures") {
  GIVEN("Circuit with PhasedISWAP") {
    // Directly simulate a circuit containing PhasedISWAP;
    // copied from a specific failing proptest.
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::PhasedISWAP, {0.509675, 1.34623}, {0, 1});
    const auto u1 = tket_sim::get_unitary(circ);
    Transforms::synthesise_tket().apply(circ);
    const auto u2 = tket_sim::get_unitary(circ);
    CHECK(u1.isApprox(u2));

    // Let's recreate the transformed circuit "manually".
    Circuit manual_circ(2);
    // These numbers were copied from printing out "circ",
    // however they are not given to many significant figures.
    manual_circ.add_op<unsigned>(OpType::U3, {0.5, 1.5, 1.00968}, {0});
    manual_circ.add_op<unsigned>(OpType::U3, {0.5, 1.5, 1.99032}, {1});
    manual_circ.add_op<unsigned>(OpType::CX, {0, 1});
    manual_circ.add_op<unsigned>(OpType::U3, {3.32688, 1.5, 0.5}, {0});
    manual_circ.add_op<unsigned>(OpType::U1, 1.32688, {1});
    manual_circ.add_op<unsigned>(OpType::CX, {0, 1});
    manual_circ.add_op<unsigned>(OpType::U3, {3.5, 0.990325, 0.5}, {0});
    manual_circ.add_op<unsigned>(OpType::U3, {3.5, 0.009675, 0.5}, {1});
    manual_circ.add_phase(0.3365575);
    const Eigen::MatrixXcd u3 = tket_sim::get_unitary(manual_circ);
    // Because the numerical values above are not very accurate,
    // the matrices don't match as closely as usual.
    CHECK(u1.isApprox(u3, 1e-4));
  }
}

// For testing a unitary via an equivalent circuit,
// requiring that the type doesn't occur in the circuit.
static Eigen::MatrixXcd get_unitary_without_op_type(
    const Circuit& circ, OpType type) {
  for (const auto& command : circ.get_commands()) {
    const auto& op_ptr = command.get_op_ptr();
    REQUIRE(op_ptr);
    REQUIRE(op_ptr->get_type() != type);
  }
  return tket_sim::get_unitary(circ);
}

// Gates which were not simulated much or at all in older tests.
SCENARIO("Check single gates") {
  GIVEN("BRIDGE") {
    const auto dense_unitary1 =
        GateUnitaryMatrix::get_unitary(OpType::BRIDGE, 3, {});
    const auto dense_unitary2 = get_unitary_without_op_type(
        CircPool::BRIDGE_using_CX_0(), OpType::BRIDGE);
    CHECK(dense_unitary1.isApprox(dense_unitary2));
  }
  GIVEN("CSWAP") {
    const auto dense_unitary1 =
        GateUnitaryMatrix::get_unitary(OpType::CSWAP, 3, {});
    const auto swap = GateUnitaryMatrix::get_unitary(OpType::SWAP, 2, {});
    const auto dense_unitary2 = internal::GateUnitaryMatrixUtils::
        get_multi_controlled_gate_dense_unitary(swap, 3);
    CHECK(dense_unitary1.isApprox(dense_unitary2));
  }
  GIVEN("PhaseGadget") {
    const double t = -1.23456789;
    const Expr t_expr = t;
    for (unsigned n_qubits = 1; n_qubits <= 4; ++n_qubits) {
      const auto dense_unitary1 =
          GateUnitaryMatrix::get_unitary(OpType::PhaseGadget, n_qubits, {t});
      const Circuit circ = phase_gadget(n_qubits, t_expr);
      const auto dense_unitary2 =
          get_unitary_without_op_type(circ, OpType::PhaseGadget);
      CHECK(dense_unitary1.isApprox(dense_unitary2));
    }
  }
}

SCENARIO("Match single gate unitaries against circuit simulator") {
  const auto& data = internal::GatesData::get().input_data;
  std::vector<Expr> parameters;
  std::vector<double> parameter_doubles;
  std::vector<unsigned> qubits;
  Eigen::MatrixXcd gate_unitary;
  Eigen::MatrixXcd circuit_sim_unitary;

  for (const auto& outer_entry : data) {
    const auto number_of_qubits = outer_entry.first;
    qubits.resize(number_of_qubits);
    std::iota(qubits.begin(), qubits.end(), 0);

    for (const auto& inner_entry : outer_entry.second) {
      parameters.resize(inner_entry.first);
      parameter_doubles.resize(inner_entry.first);
      for (unsigned nn = 0; nn < parameters.size(); ++nn) {
        parameter_doubles[nn] = 0.123456789 + nn * 0.23456789;
        parameters[nn] = parameter_doubles[nn];
      }
      for (OpType type : inner_entry.second) {
        gate_unitary = GateUnitaryMatrix::get_unitary(
            type, number_of_qubits, parameter_doubles);
        Circuit circ(number_of_qubits);
        circ.add_op<unsigned>(type, parameters, qubits);
        circuit_sim_unitary = tket_sim::get_unitary(circ);
        CHECK(gate_unitary.isApprox(circuit_sim_unitary));
      }
    }
  }
}

SCENARIO("Unitaries for controlled operations") {
  GIVEN("CRz") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CRz, 0.5, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = exp(-0.25 * i_ * PI);
    V(3, 3) = exp(0.25 * i_ * PI);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CRx") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CRx, 0.5, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = std::cos(0.25 * PI);
    V(2, 3) = i_ * std::sin(-0.25 * PI);
    V(3, 2) = i_ * std::sin(-0.25 * PI);
    V(3, 3) = std::cos(0.25 * PI);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CRy") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CRy, 0.5, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = std::cos(0.25 * PI);
    V(2, 3) = std::sin(-0.25 * PI);
    V(3, 2) = std::sin(0.25 * PI);
    V(3, 3) = std::cos(0.25 * PI);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CV") {
    const double sq = 1 / std::sqrt(2.);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CV, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = sq;
    V(2, 3) = sq * -i_;
    V(3, 2) = sq * -i_;
    V(3, 3) = sq;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CVdg") {
    const double sq = 1 / std::sqrt(2.);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CVdg, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = sq;
    V(2, 3) = sq * i_;
    V(3, 2) = sq * i_;
    V(3, 3) = sq;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CSX") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CSX, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = 0.5 * (1. + i_);
    V(2, 3) = 0.5 * (1. - i_);
    V(3, 2) = 0.5 * (1. - i_);
    V(3, 3) = 0.5 * (1. + i_);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CSXdg") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CSXdg, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = 0.5 * (1. - i_);
    V(2, 3) = 0.5 * (1. + i_);
    V(3, 2) = 0.5 * (1. + i_);
    V(3, 3) = 0.5 * (1. - i_);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CU1") {
    double a = 0.125;
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CU1, a, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(3, 3) = exp(i_ * PI * a);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CU3") {
    double a = 0.125, b = 0.375, c = 0.75;
    Circuit circ0(1);
    circ0.add_op<unsigned>(OpType::U3, {a, b, c}, {0});
    const Eigen::MatrixXcd U0 = tket_sim::get_unitary(circ0);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CU3, {a, b, c}, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = U0(0, 0);
    V(2, 3) = U0(0, 1);
    V(3, 2) = U0(1, 0);
    V(3, 3) = U0(1, 1);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CCX") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(8, 8);
    V(6, 6) = V(7, 7) = 0.;
    V(6, 7) = V(7, 6) = 1.;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CnX") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(16, 16);
    V(14, 14) = V(15, 15) = 0.;
    V(14, 15) = V(15, 14) = 1.;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("CnRy") {
    double a = 0.125;
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CnRy, a, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    double x = cos(0.5 * PI * a), y = sin(0.5 * PI * a);
    V(2, 2) = V(3, 3) = x;
    V(2, 3) = -y;
    V(3, 2) = y;
    REQUIRE(U.isApprox(V));
  }
}

SCENARIO("Handling internal qubit permutations") {
  GIVEN("A Clifford reduction introducing a wireswap") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit circ2(3);
    circ2.add_op<unsigned>(OpType::CX, {1, 0});
    circ2.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ2.replace_SWAPs();
    Eigen::MatrixXcd m1 = tket_sim::get_unitary(circ);
    Eigen::MatrixXcd m2 = tket_sim::get_unitary(circ2);
    REQUIRE(m1.isApprox(m2, ERR_EPS));
  }
}

SCENARIO("compare_statevectors_or_unitaries gives expected errors") {
  const std::array<tket_sim::MatrixEquivalence, 2> equivalences{
      tket_sim::MatrixEquivalence::EQUAL,
      tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE};
  const auto check_error = [&equivalences](
                               const Eigen::MatrixXcd& matr1,
                               const Eigen::MatrixXcd& matr2,
                               const std::string& message) {
    for (auto equiv : equivalences) {
      try {
        const bool equivalent =
            tket_sim::compare_statevectors_or_unitaries(matr1, matr2, equiv);
        INFO("matrices compared equivalent: " << equivalent);
        CHECK(false);
      } catch (const std::exception& e) {
        CHECK_THAT(e.what(), ContainsSubstring(message));
      }
    }
  };
  GIVEN("Non-square, non-statevector inputs") {
    const Eigen::MatrixXcd matr = Eigen::MatrixXcd::Identity(4, 2);
    check_error(matr, matr, "Not square, and also not column vectors");
  }
  GIVEN("Wrongly sized unitary inputs") {
    const Eigen::MatrixXcd matr = Eigen::MatrixXcd::Identity(3, 3);
    check_error(matr, matr, "matrix size 3 is not a power of two");
  }
  GIVEN("Wrongly sized statevector inputs") {
    const Eigen::MatrixXcd matr = Eigen::MatrixXcd::Identity(5, 1);
    check_error(matr, matr, "matrix size 5 is not a power of two");
  }
  GIVEN("Different sized unitary inputs") {
    const Eigen::MatrixXcd matr1 = Eigen::MatrixXcd::Identity(2, 2);
    const Eigen::MatrixXcd matr2 = Eigen::MatrixXcd::Identity(4, 4);
    check_error(matr1, matr2, "Different sized matrices");
  }
  GIVEN("Different sized statevector inputs") {
    const Eigen::MatrixXcd matr1 = Eigen::MatrixXcd::Identity(2, 1);
    const Eigen::MatrixXcd matr2 = Eigen::MatrixXcd::Identity(4, 1);
    check_error(matr1, matr2, "Different sized matrices");
  }
  GIVEN("Not norm 1 statevectors") {
    Eigen::MatrixXcd matr(4, 1);
    matr << 1, 2, 3, 4;
    check_error(matr, matr, "State vector is not normalised");
  }
  GIVEN("Non unitary inputs") {
    Eigen::MatrixXcd matr(2, 2);
    matr << 1, 2, 3, 4;
    check_error(matr, matr, "Matrix is not unitary");
  }
}

// We just want to avoid any simple pattern; doesn't matter exactly what.
static Eigen::MatrixXcd get_random_matrix(unsigned rows, unsigned cols) {
  Eigen::MatrixXcd matr(rows, cols);
  // Not very important, it just makes the Frobenius norm equal 1.
  const double rr = 1.0 / std::sqrt(rows * cols);
  for (unsigned jj = 0; jj < cols; ++jj) {
    for (unsigned ii = 0; ii < rows; ++ii) {
      matr(ii, jj) = std::polar(
          rr, 0.15 + 0.1 * cols * ii + 0.2 * rows * jj + 0.03 * ii * ii +
                  0.04 * jj * jj);
    }
  }
  return matr;
}

// Returns almost, but not exactly,  cM  for some complex number
// c with |c| = 1.
static Eigen::MatrixXcd get_almost_phase_equivalent_matrix(
    const Eigen::MatrixXcd& matr) {
  double theta = 1.23456789;
  const double dtheta = 1e-12 / (matr.rows() * matr.cols());
  Eigen::MatrixXcd new_matr(matr);
  for (unsigned jj = 0; jj < matr.cols(); ++jj) {
    for (unsigned ii = 0; ii < matr.rows(); ++ii) {
      new_matr(ii, jj) *= std::polar(1.0, theta);
      theta += dtheta;
    }
  }
  return new_matr;
}

SCENARIO(
    "compare_statevectors_or_unitaries works as expected for valid "
    "inputs") {
  Eigen::VectorXcd norm_one_vect(4);
  for (unsigned ii = 0; ii < norm_one_vect.rows(); ++ii) {
    // "random" entries.
    norm_one_vect(ii) = std::polar(0.5, 0.2 * (ii + 1) * (ii + 2));
  }
  // Householder matrix - more interesting than just a diagonal matrix
  const Eigen::MatrixXcd unitary =
      Eigen::MatrixXcd::Identity(4, 4) -
      2.0 * norm_one_vect * norm_one_vect.adjoint();

  const Eigen::VectorXcd vect_entries = get_random_matrix(4, 1);
  const auto matr_entries = get_random_matrix(4, 4);

  const double small_eps = 1e-12;
  const double large_eps = 1e-4;

  const Eigen::VectorXcd almost_equal_vect =
      norm_one_vect + small_eps * vect_entries;
  const auto almost_equal_unitary = unitary + small_eps * matr_entries;

  // Fairly close, but definitely different
  Eigen::VectorXcd different_vect = norm_one_vect + large_eps * vect_entries;
  different_vect /= different_vect.blueNorm();

  const Eigen::MatrixXcd different_matr =
      Eigen::MatrixXcd::Identity(4, 4) -
      2.0 * different_vect * different_vect.adjoint();

  const std::array<tket_sim::MatrixEquivalence, 2> equivalences{
      tket_sim::MatrixEquivalence::EQUAL,
      tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE};

  const std::array<std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd>, 2>
      almost_equal_pairs{
          std::make_pair<Eigen::MatrixXcd, Eigen::MatrixXcd>(
              norm_one_vect, almost_equal_vect),
          {unitary, almost_equal_unitary}};

  GIVEN("Equal up to roundoff") {
    for (const auto& entry : almost_equal_pairs) {
      for (auto equiv : equivalences) {
        CHECK(tket_sim::compare_statevectors_or_unitaries(
            entry.first, entry.second, equiv));
      }
      const auto almost_phase_equivalent_matrix =
          get_almost_phase_equivalent_matrix(entry.second);
      CHECK(!tket_sim::compare_statevectors_or_unitaries(
          entry.first, almost_phase_equivalent_matrix));
      CHECK(tket_sim::compare_statevectors_or_unitaries(
          entry.first, almost_phase_equivalent_matrix,
          tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE));
    }
  }
  GIVEN("Different inputs") {
    for (auto equiv : equivalences) {
      CHECK(!tket_sim::compare_statevectors_or_unitaries(
          norm_one_vect, different_vect, equiv));
      CHECK(!tket_sim::compare_statevectors_or_unitaries(
          unitary, different_matr, equiv));
      CHECK(!tket_sim::compare_statevectors_or_unitaries(
          norm_one_vect, get_almost_phase_equivalent_matrix(different_vect),
          equiv));
      CHECK(!tket_sim::compare_statevectors_or_unitaries(
          unitary, get_almost_phase_equivalent_matrix(different_matr), equiv));
    }
  }
}

}  // namespace test_CircuitSimulator
}  // namespace tket
