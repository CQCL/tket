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

#include "MeasurementSetup/MeasurementSetup.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_MeasurementSetup {

SCENARIO("verify_measurement_setup") {
  GIVEN("An empty setup") {
    MeasurementSetup ms;
    REQUIRE(ms.verify());
  }
  GIVEN("A basic Z-measure circuit") {
    MeasurementSetup ms;
    Circuit mc(2, 2);
    mc.add_measure(0, 0);
    mc.add_measure(1, 1);
    ms.add_measurement_circuit(mc);
    Qubit q0(q_default_reg(), 0);
    Qubit q1(q_default_reg(), 1);
    QubitPauliString ii;
    QubitPauliString zi({{q0, Pauli::Z}});
    QubitPauliString iz({{q1, Pauli::Z}});
    QubitPauliString zz({{q0, Pauli::Z}, {q1, Pauli::Z}});
    ms.add_result_for_term(ii, {0, {}, false});
    ms.add_result_for_term(zi, {0, {0}, false});
    ms.add_result_for_term(iz, {0, {1}, false});
    ms.add_result_for_term(zz, {0, {0, 1}, false});
    REQUIRE(ms.verify());
  }
  GIVEN("Multiple circuits with same Pauli") {
    MeasurementSetup ms;
    Circuit mc0(2, 2);
    mc0.add_measure(0, 0);
    mc0.add_measure(1, 1);
    ms.add_measurement_circuit(mc0);
    Circuit mc1(2, 2);
    mc1.add_measure(0, 0);
    mc1.add_op<unsigned>(OpType::H, {1});
    mc1.add_measure(1, 1);
    ms.add_measurement_circuit(mc1);
    Qubit q0(q_default_reg(), 0);
    Qubit q1(q_default_reg(), 1);
    QubitPauliString zi({{q0, Pauli::Z}});
    ms.add_result_for_term(zi, {0, {0}, false});
    ms.add_result_for_term(zi, {1, {0}, false});
    REQUIRE(ms.verify());
  }
  GIVEN("Parity flips") {
    MeasurementSetup ms;
    Circuit mc(1, 1);
    mc.add_op<unsigned>(OpType::X, {0});
    mc.add_measure(0, 0);
    ms.add_measurement_circuit(mc);
    Qubit q0(q_default_reg(), 0);
    QubitPauliString z({{q0, Pauli::Z}});
    ms.add_result_for_term(z, {0, {0}, true});
    REQUIRE(ms.verify());
  }
  GIVEN("Wrong ones") {
    MeasurementSetup ms;
    Circuit mc(2, 2);
    mc.add_op<unsigned>(OpType::X, {0});
    mc.add_measure(0, 0);
    mc.add_op<unsigned>(OpType::V, {1});
    mc.add_measure(1, 1);
    ms.add_measurement_circuit(mc);
    Qubit q0(q_default_reg(), 0);
    Qubit q1(q_default_reg(), 1);
    WHEN("Wrong parity") {
      QubitPauliString zi({{q0, Pauli::Z}});
      ms.add_result_for_term(zi, {0, {0}, false});
      REQUIRE_FALSE(ms.verify());
    }
    WHEN("Wrong string") {
      QubitPauliString ix({{q1, Pauli::X}});
      ms.add_result_for_term(ix, {0, {1}, false});
      REQUIRE_FALSE(ms.verify());
    }
    WHEN("Wrong bit set") {
      QubitPauliString iy({{q1, Pauli::Y}});
      ms.add_result_for_term(iy, {0, {0, 1}, false});
      REQUIRE_FALSE(ms.verify());
    }
  }
  GIVEN("HQS experiment") {
    MeasurementSetup ms;
    const auto add_meas_circ = [&ms](Circuit& circ) {
      for (unsigned nn = 0; nn <= 3; ++nn) {
        circ.add_measure(nn, nn);
      }
      ms.add_measurement_circuit(circ);
    };

    Circuit mc0(4, 4);
    mc0.add_op<unsigned>(OpType::CX, {0, 2});
    mc0.add_op<unsigned>(OpType::H, {1});
    mc0.add_op<unsigned>(OpType::CX, {1, 0});
    mc0.add_op<unsigned>(OpType::H, {1});
    add_meas_circ(mc0);

    Circuit mc1(4, 4);
    add_meas_circ(mc1);

    Circuit mc2(4, 4);
    mc2.add_op<unsigned>(OpType::V, {0});
    mc2.add_op<unsigned>(OpType::H, {1});
    mc2.add_op<unsigned>(OpType::V, {2});
    mc2.add_op<unsigned>(OpType::V, {3});
    add_meas_circ(mc2);

    Circuit mc3(4, 4);
    add_1qb_gates(mc3, OpType::V, {0, 2, 3});
    add_2qb_gates(mc3, OpType::CX, {{0, 1}, {0, 2}});
    mc3.add_op<unsigned>(OpType::H, {0});
    add_meas_circ(mc3);

    Circuit mc4(4, 4);
    add_1qb_gates(mc4, OpType::H, {0, 1, 2, 3});
    add_meas_circ(mc4);

    Circuit mc5(4, 4);
    mc5.add_op<unsigned>(OpType::CX, {0, 2});
    mc5.add_op<unsigned>(OpType::H, {0});
    mc5.add_op<unsigned>(OpType::V, {1});
    mc5.add_op<unsigned>(OpType::V, {2});
    mc5.add_op<unsigned>(OpType::CX, {1, 2});
    mc5.add_op<unsigned>(OpType::V, {1});
    mc5.add_op<unsigned>(OpType::H, {3});
    add_meas_circ(mc5);

    Qubit q0(q_default_reg(), 0);
    Qubit q1(q_default_reg(), 1);
    Qubit q2(q_default_reg(), 2);
    Qubit q3(q_default_reg(), 3);

    QubitPauliTensor x0(q0, Pauli::X);
    QubitPauliTensor y0(q0, Pauli::Y);
    QubitPauliTensor z0(q0, Pauli::Z);
    QubitPauliTensor x1(q1, Pauli::X);
    QubitPauliTensor y1(q1, Pauli::Y);
    QubitPauliTensor z1(q1, Pauli::Z);
    QubitPauliTensor x2(q2, Pauli::X);
    QubitPauliTensor y2(q2, Pauli::Y);
    QubitPauliTensor z2(q2, Pauli::Z);
    QubitPauliTensor x3(q3, Pauli::X);
    QubitPauliTensor y3(q3, Pauli::Y);
    QubitPauliTensor z3(q3, Pauli::Z);

    ms.add_result_for_term(z0, {1, {0}, false});
    ms.add_result_for_term(z0 * z1, {1, {0, 1}, false});
    ms.add_result_for_term(z1, {1, {1}, false});
    ms.add_result_for_term(x0 * y1 * y2, {0, {0, 1, 2}, false});
    ms.add_result_for_term(y0 * x1 * y2, {2, {0, 1, 2}, false});
    ms.add_result_for_term(y0 * y1 * x2, {0, {0, 1}, false});
    ms.add_result_for_term(x0 * x1 * x2, {3, {0}, false});
    ms.add_result_for_term(z0 * x1, {0, {0}, false});
    ms.add_result_for_term(x1, {2, {1}, false});
    ms.add_result_for_term(x1 * z2, {0, {0, 2}, false});
    ms.add_result_for_term(z0 * x1 * z2, {3, {0, 2}, true});
    ms.add_result_for_term(y0 * z1 * y2, {0, {1, 2}, true});
    ms.add_result_for_term(x0 * z1 * x2, {0, {1}, false});
    ms.add_result_for_term(y0 * y3, {2, {0, 3}, false});
    ms.add_result_for_term(x0 * z1 * x3, {5, {0, 1, 3}, true});
    ms.add_result_for_term(x0 * x3, {4, {0, 3}, false});
    ms.add_result_for_term(y0 * z1 * y3, {3, {1, 3}, false});
    ms.add_result_for_term(z2, {1, {2}, false});
    ms.add_result_for_term(z0 * z2, {0, {2}, false});
    ms.add_result_for_term(z1 * z2, {1, {1, 2}, false});
    ms.add_result_for_term(z0 * z1 * z2, {1, {0, 1, 2}, false});
    ms.add_result_for_term(x0 * y1 * z2 * y3, {3, {0, 1, 2, 3}, false});
    ms.add_result_for_term(x0 * x1 * x3, {4, {0, 1, 3}, false});
    ms.add_result_for_term(y0 * y1 * z2 * x3, {5, {0, 2, 3}, false});
    ms.add_result_for_term(y0 * x1 * y3, {2, {0, 1, 3}, false});
    ms.add_result_for_term(z3, {0, {3}, false});
    ms.add_result_for_term(z0 * z3, {1, {0, 3}, false});
    ms.add_result_for_term(x0 * y1 * y2 * z3, {0, {0, 1, 2, 3}, false});
    ms.add_result_for_term(y0 * y1 * x2 * z3, {0, {0, 1, 3}, false});
    ms.add_result_for_term(z0 * y1 * x2 * y3, {3, {0, 1, 3}, false});
    ms.add_result_for_term(z0 * y1 * y2 * x3, {5, {2, 3}, false});
    ms.add_result_for_term(x1 * x2 * x3, {4, {1, 2, 3}, false});
    ms.add_result_for_term(x1 * y2 * y3, {2, {1, 2, 3}, false});
    ms.add_result_for_term(z0 * z1 * z3, {1, {0, 1, 3}, false});
    ms.add_result_for_term(z0 * x1 * z3, {0, {0, 3}, false});
    ms.add_result_for_term(x1 * z2 * z3, {0, {0, 2, 3}, false});
    ms.add_result_for_term(y2 * y3, {2, {2, 3}, false});
    ms.add_result_for_term(z1 * x2 * x3, {5, {1, 3}, true});
    ms.add_result_for_term(x2 * x3, {4, {2, 3}, false});
    ms.add_result_for_term(z1 * y2 * y3, {3, {1, 2, 3}, false});
    ms.add_result_for_term(z2 * z3, {1, {2, 3}, false});
    ms.add_result_for_term(z1 * z2 * z3, {1, {1, 2, 3}, false});

    REQUIRE(ms.verify());
  }
  GIVEN("H3 Singlet JW") {
    MeasurementSetup ms;
    const auto add_meas_circ = [&ms](Circuit& circ) {
      for (unsigned nn = 0; nn <= 5; ++nn) {
        circ.add_measure(nn, nn);
      }
      ms.add_measurement_circuit(circ);
    };

    Circuit mc0(6, 6);
    mc0.add_op<unsigned>(OpType::H, {0});
    mc0.add_op<unsigned>(OpType::V, {1});
    mc0.add_op<unsigned>(OpType::H, {2});
    mc0.add_op<unsigned>(OpType::V, {3});
    add_2qb_gates(mc0, OpType::CX, {{0, 1}, {2, 3}, {4, 5}, {0, 2}});
    mc0.add_op<unsigned>(OpType::H, {0});
    mc0.add_op<unsigned>(OpType::V, {4});
    add_meas_circ(mc0);

    Circuit mc1(6, 6);
    add_2qb_gates(mc1, OpType::CX, {{0, 2}, {3, 5}, {0, 3}});
    mc1.add_op<unsigned>(OpType::H, {0});
    add_meas_circ(mc1);

    Circuit mc2(6, 6);
    mc2.add_op<unsigned>(OpType::V, {0});
    mc2.add_op<unsigned>(OpType::V, {4});
    add_2qb_gates(mc2, OpType::CX, {{2, 1}, {2, 3}, {0, 4}});
    mc2.add_op<unsigned>(OpType::V, {2});
    mc2.add_op<unsigned>(OpType::CX, {0, 2});
    mc2.add_op<unsigned>(OpType::V, {0});
    add_meas_circ(mc2);

    Circuit mc3(6, 6);
    mc3.add_op<unsigned>(OpType::V, {0});
    mc3.add_op<unsigned>(OpType::H, {1});
    mc3.add_op<unsigned>(OpType::V, {2});
    mc3.add_op<unsigned>(OpType::H, {3});
    add_2qb_gates(mc3, OpType::CX, {{0, 1}, {2, 3}, {4, 5}, {0, 2}});
    mc3.add_op<unsigned>(OpType::H, {0});
    mc3.add_op<unsigned>(OpType::V, {4});
    add_meas_circ(mc3);

    Circuit mc4(6, 6);
    add_2qb_gates(mc4, OpType::CX, {{0, 1}, {2, 3}, {4, 5}});
    add_1qb_gates(mc4, OpType::H, {0, 2, 4});
    add_meas_circ(mc4);

    Circuit mc5(6, 6);
    add_1qb_gates(mc5, OpType::V, {0, 3});
    add_2qb_gates(mc5, OpType::CX, {{1, 2}, {4, 5}});
    mc5.add_op<unsigned>(OpType::V, {1});
    add_1qb_gates(mc5, OpType::H, {4, 5});
    add_2qb_gates(mc5, OpType::CX, {{0, 3}, {0, 5}, {1, 5}});
    add_1qb_gates(mc5, OpType::H, {0, 1});
    add_meas_circ(mc5);

    Circuit mc6(6, 6);
    add_2qb_gates(mc6, OpType::CX, {{0, 4}, {1, 2}, {3, 5}, {4, 5}});
    add_1qb_gates(mc6, OpType::H, {0, 1, 3, 4});
    add_meas_circ(mc6);

    Circuit mc7(6, 6);
    mc7.add_op<unsigned>(OpType::CX, {1, 5});
    mc7.add_op<unsigned>(OpType::H, {1});
    add_meas_circ(mc7);

    Circuit mc8(6, 6);
    mc8.add_op<unsigned>(OpType::CX, {0, 2});
    mc8.add_op<unsigned>(OpType::CX, {3, 1});
    mc8.add_op<unsigned>(OpType::H, {0});
    mc8.add_op<unsigned>(OpType::CX, {0, 3});
    mc8.add_op<unsigned>(OpType::H, {0});
    add_meas_circ(mc8);

    Qubit q0(q_default_reg(), 0);
    Qubit q1(q_default_reg(), 1);
    Qubit q2(q_default_reg(), 2);
    Qubit q3(q_default_reg(), 3);
    Qubit q4(q_default_reg(), 4);
    Qubit q5(q_default_reg(), 5);

    QubitPauliTensor x0(q0, Pauli::X);
    QubitPauliTensor y0(q0, Pauli::Y);
    QubitPauliTensor z0(q0, Pauli::Z);
    QubitPauliTensor x1(q1, Pauli::X);
    QubitPauliTensor y1(q1, Pauli::Y);
    QubitPauliTensor z1(q1, Pauli::Z);
    QubitPauliTensor x2(q2, Pauli::X);
    QubitPauliTensor y2(q2, Pauli::Y);
    QubitPauliTensor z2(q2, Pauli::Z);
    QubitPauliTensor x3(q3, Pauli::X);
    QubitPauliTensor y3(q3, Pauli::Y);
    QubitPauliTensor z3(q3, Pauli::Z);
    QubitPauliTensor x4(q4, Pauli::X);
    QubitPauliTensor y4(q4, Pauli::Y);
    QubitPauliTensor z4(q4, Pauli::Z);
    QubitPauliTensor x5(q5, Pauli::X);
    QubitPauliTensor y5(q5, Pauli::Y);
    QubitPauliTensor z5(q5, Pauli::Z);

    ms.add_result_for_term(y0 * z1 * y2 * z3, {0, {0, 1, 3}, false});
    ms.add_result_for_term(
        y0 * z1 * z2 * x3 * x4 * y5, {0, {0, 1, 4, 5}, true});
    ms.add_result_for_term(y1 * y3, {0, {1, 2, 3}, false});
    ms.add_result_for_term(y1 * x2 * x4 * y5, {0, {1, 2, 4, 5}, false});
    ms.add_result_for_term(x0 * y1 * y4 * x5, {0, {1, 4}, false});
    ms.add_result_for_term(x2 * y3 * y4 * x5, {0, {3, 4}, false});
    ms.add_result_for_term(z4 * z5, {0, {5}, false});
    ms.add_result_for_term(z4 * z5, {3, {5}, false});
    ms.add_result_for_term(z4 * z5, {4, {5}, false});
    ms.add_result_for_term(z4 * z5, {8, {4, 5}, false});
    ms.add_result_for_term(
        y0 * z1 * y2 * y3 * z4 * y5, {1, {0, 1, 2, 4, 5}, false});
    ms.add_result_for_term(x0 * z1 * x2 * x3 * z4 * x5, {1, {0, 1, 4}, false});
    ms.add_result_for_term(
        y0 * z1 * y2 * x3 * z4 * x5, {1, {0, 1, 2, 4}, true});
    ms.add_result_for_term(z1, {1, {1}, false});
    ms.add_result_for_term(z3 * z5, {1, {5}, false});
    ms.add_result_for_term(z1 * z4, {1, {1, 4}, false});
    ms.add_result_for_term(
        x0 * z1 * x2 * y3 * z4 * y5, {1, {0, 1, 4, 5}, true});
    ms.add_result_for_term(z0 * z5, {1, {3, 5}, false});
    ms.add_result_for_term(z2 * z5, {1, {2, 3, 5}, false});
    ms.add_result_for_term(z0 * z3, {1, {3}, false});
    ms.add_result_for_term(z0 * z2, {1, {2}, false});
    ms.add_result_for_term(z4, {1, {4}, false});
    ms.add_result_for_term(z4, {7, {4}, false});
    ms.add_result_for_term(z4, {8, {4}, false});
    ms.add_result_for_term(z2 * z3, {1, {2, 3}, false});
    ms.add_result_for_term(x1 * x2 * y3 * y4, {2, {2, 3, 4}, false});
    ms.add_result_for_term(
        x0 * z1 * z2 * z3 * x4 * z5, {2, {0, 1, 2, 3, 5}, true});
    ms.add_result_for_term(y0 * x1 * x2 * y3, {2, {2, 3}, false});
    ms.add_result_for_term(y1 * x2 * x3 * y4, {2, {1, 2, 4}, false});
    ms.add_result_for_term(y0 * z1 * z2 * y4, {2, {1, 4}, false});
    ms.add_result_for_term(z5, {2, {5}, false});
    ms.add_result_for_term(z5, {8, {5}, false});
    ms.add_result_for_term(z1 * z3, {2, {1, 3}, false});
    ms.add_result_for_term(y0 * y1 * x2 * x3, {2, {1, 2}, false});
    ms.add_result_for_term(z1 * z2, {2, {1}, false});
    ms.add_result_for_term(x1 * y2 * y4 * x5, {3, {1, 2, 4}, false});
    ms.add_result_for_term(x0 * z1 * x2 * z3, {3, {0}, false});
    ms.add_result_for_term(x0 * z1 * z2 * y3 * y4 * x5, {3, {0, 3, 4}, true});
    ms.add_result_for_term(y0 * x1 * x4 * y5, {3, {1, 4, 5}, false});
    ms.add_result_for_term(x1 * x3, {3, {1, 2, 3}, false});
    ms.add_result_for_term(y2 * x3 * x4 * y5, {3, {3, 4, 5}, false});
    ms.add_result_for_term(z0 * z1, {4, {1}, false});
    ms.add_result_for_term(y0 * y1 * x4 * x5, {4, {0, 1, 4}, true});
    ms.add_result_for_term(x0 * x1 * y2 * y3, {4, {0, 2, 3}, true});
    ms.add_result_for_term(x2 * x3 * y4 * y5, {4, {2, 4, 5}, true});
    ms.add_result_for_term(x0 * x1 * y4 * y5, {4, {0, 4, 5}, true});
    ms.add_result_for_term(y2 * y3 * x4 * x5, {4, {2, 3, 4}, true});
    ms.add_result_for_term(x0 * z1 * z2 * x3 * y4 * y5, {5, {0, 2, 4}, true});
    ms.add_result_for_term(x1 * y2 * y3 * x4, {5, {2, 3, 4, 5}, false});
    ms.add_result_for_term(y0 * z1 * z2 * y3 * x4 * x5, {5, {2, 3, 4}, false});
    ms.add_result_for_term(y1 * y2 * y4 * y5, {5, {1, 2, 4}, false});
    ms.add_result_for_term(x0 * y1 * y2 * x3, {5, {0, 1, 2}, true});
    ms.add_result_for_term(y0 * z1 * z2 * z3 * y4 * z5, {6, {0, 2, 5}, true});
    ms.add_result_for_term(x1 * x2 * x4 * x5, {6, {1, 4}, false});
    ms.add_result_for_term(y1 * y2 * x3 * x4, {6, {1, 2, 3, 4}, true});
    ms.add_result_for_term(x0 * z1 * z2 * x4, {6, {0, 2}, false});
    ms.add_result_for_term(y1 * z3 * z4 * y5, {7, {1, 3, 4, 5}, true});
    ms.add_result_for_term(y1 * z2 * z3 * y5, {7, {1, 2, 3, 5}, true});
    ms.add_result_for_term(z3, {7, {3}, false});
    ms.add_result_for_term(x1 * z3 * z4 * x5, {7, {1, 3, 4}, false});
    ms.add_result_for_term(z1 * z5, {7, {5}, false});
    ms.add_result_for_term(x1 * z2 * z3 * x5, {7, {1, 2, 3}, false});
    ms.add_result_for_term(z3 * z4, {7, {3, 4}, false});
    ms.add_result_for_term(z0, {7, {0}, false});
    ms.add_result_for_term(z2, {7, {2}, false});
    ms.add_result_for_term(z0 * z4, {7, {0, 4}, false});
    ms.add_result_for_term(z2 * z4, {7, {2, 4}, false});
    ms.add_result_for_term(y1 * z2 * y3 * z4, {8, {0, 1, 2, 4}, true});
    ms.add_result_for_term(x1 * z2 * x3 * z4, {8, {0, 2, 4}, false});
    ms.add_result_for_term(y0 * z1 * y2 * z5, {8, {1, 2, 3, 5}, true});
    ms.add_result_for_term(x0 * z1 * x2 * z5, {8, {1, 3, 5}, false});

    REQUIRE(ms.verify());
  }
}

}  // namespace test_MeasurementSetup
}  // namespace tket
