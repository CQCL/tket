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

#include "MeasurementSetup/MeasurementReduction.hpp"

namespace tket {
namespace test_MeasurementReduction {

SCENARIO("Some QubitOperators") {
  GIVEN("4 string QubitOperator") {
    QubitPauliString qp_map0(Qubit(0), Pauli::I);
    QubitPauliString qp_map1(Qubit(0), Pauli::X);
    QubitPauliString qp_map2(Qubit(0), Pauli::Y);
    QubitPauliString qp_map3(Qubit(0), Pauli::Z);
    std::list<QubitPauliString> pts{qp_map0, qp_map1, qp_map2, qp_map3};

    WHEN("Commuting sets") {
      MeasurementSetup measurements =
          measurement_reduction(pts, PauliPartitionStrat::CommutingSets);
      REQUIRE(measurements.get_circs().size() == 3);
      for (const Circuit& circ : measurements.get_circs()) {
        REQUIRE(circ.all_bits()[0] == Bit(0));
      }
      REQUIRE(measurements.verify());
    }
    WHEN("Nonconflicting sets") {
      MeasurementSetup measurements =
          measurement_reduction(pts, PauliPartitionStrat::NonConflictingSets);
      REQUIRE(measurements.get_circs().size() == 3);
      REQUIRE(measurements.verify());
    }
  }
  GIVEN("7 string QubitOperator") {
    QubitPauliString qp_map0(Qubit(0), Pauli::Z);
    QubitPauliString qp_map1(Qubit(1), Pauli::Z);
    QubitPauliString qp_map2(Qubit(2), Pauli::Z);
    QubitPauliString qp_map3(Qubit(3), Pauli::Z);
    QubitPauliString qp_map4({Pauli::Z, Pauli::Z, Pauli::Z, Pauli::Z});
    QubitPauliString qp_map5({Pauli::X, Pauli::X, Pauli::Y, Pauli::Y});
    QubitPauliString qp_map6({Pauli::Y, Pauli::Y, Pauli::X, Pauli::X});
    std::list<QubitPauliString> pts{qp_map0, qp_map1, qp_map2, qp_map3,
                                    qp_map4, qp_map5, qp_map6};

    WHEN("Commuting sets") {
      MeasurementSetup measurements =
          measurement_reduction(pts, PauliPartitionStrat::CommutingSets);
      REQUIRE(measurements.get_circs().size() == 2);
      REQUIRE(measurements.verify());
    }
    WHEN("Nonconflicting sets") {
      MeasurementSetup measurements =
          measurement_reduction(pts, PauliPartitionStrat::NonConflictingSets);
      REQUIRE(measurements.get_circs().size() == 3);
      REQUIRE(measurements.verify());
    }
  }
  GIVEN("8 strings over 4 qubits") {
    QubitPauliString qp_map0({Pauli::X, Pauli::X, Pauli::X, Pauli::Y});
    QubitPauliString qp_map1({Pauli::X, Pauli::X, Pauli::Y, Pauli::X});
    QubitPauliString qp_map2({Pauli::X, Pauli::Y, Pauli::X, Pauli::X});
    QubitPauliString qp_map3({Pauli::X, Pauli::Y, Pauli::Y, Pauli::Y});
    QubitPauliString qp_map4({Pauli::Y, Pauli::X, Pauli::X, Pauli::X});
    QubitPauliString qp_map5({Pauli::Y, Pauli::X, Pauli::Y, Pauli::Y});
    QubitPauliString qp_map6({Pauli::Y, Pauli::Y, Pauli::X, Pauli::Y});
    QubitPauliString qp_map7({Pauli::Y, Pauli::Y, Pauli::Y, Pauli::X});
    QubitPauliString qp_map8({Pauli::I, Pauli::I, Pauli::I, Pauli::I});
    std::list<QubitPauliString> pts{qp_map0, qp_map1, qp_map2, qp_map3, qp_map4,
                                    qp_map5, qp_map6, qp_map7, qp_map8};
    WHEN("Commuting sets") {
      MeasurementSetup measurements =
          measurement_reduction(pts, PauliPartitionStrat::CommutingSets);
      REQUIRE(measurements.get_circs().size() == 1);
      REQUIRE(measurements.get_circs()[0].count_gates(OpType::CX) == 3);
      REQUIRE(measurements.verify());
    }
    WHEN("Nonconflicting sets") {
      MeasurementSetup measurements =
          measurement_reduction(pts, PauliPartitionStrat::NonConflictingSets);
      REQUIRE(measurements.get_circs().size() == 8);
      for (const Circuit& circ : measurements.get_circs()) {
        REQUIRE(circ.count_gates(OpType::CX) == 0);
      }
      REQUIRE(measurements.verify());
    }
  }
}

}  // namespace test_MeasurementReduction
}  // namespace tket
