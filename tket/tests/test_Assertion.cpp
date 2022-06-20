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

#include "Circuit/AssertionSynthesis.hpp"
#include "Circuit/Circuit.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/PauliStrings.hpp"
#include "testutil.hpp"
namespace tket {

SCENARIO("Testing projector based assertion synthesis") {
  GIVEN("A projector with rank < 2 ^ n-1 and rank is a power of 2") {
    WHEN("The projector is 1q") {
      Eigen::MatrixXcd P(2, 2);
      P << 1, 0, 0, 0;
      auto [c, _] = projector_assertion_synthesis(P);
      REQUIRE(c.n_qubits() == 1);
      REQUIRE(c.count_gates(OpType::Unitary1qBox) == 2);
    }
    WHEN("The projector is 2q") {
      Eigen::MatrixXcd bell(4, 4);
      bell << 0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0.5;
      auto [c, _] = projector_assertion_synthesis(bell);
      REQUIRE(c.n_qubits() == 2);
      REQUIRE(c.count_gates(OpType::Unitary2qBox) == 2);
    }
    WHEN("The projector is 3q") {
      Eigen::MatrixXcd P = Eigen::MatrixXcd::Zero(8, 8);
      P(0, 0) = 1;
      P(1, 1) = 1;
      P(2, 2) = 1;
      P(7, 7) = 1;
      auto [c, _] = projector_assertion_synthesis(P);
      REQUIRE(c.n_qubits() == 3);
      REQUIRE(c.count_gates(OpType::Unitary3qBox) == 2);
    }
  }

  GIVEN("A projector with rank < 2 ^ n-1 and rank is not a power of 2") {
    WHEN("The projector is 3q") {
      Eigen::MatrixXcd P = Eigen::MatrixXcd::Zero(8, 8);
      P(0, 0) = 1;
      P(1, 1) = 1;
      P(7, 7) = 1;
      auto [c, _] = projector_assertion_synthesis(P);
      REQUIRE(c.n_qubits() == 3);
      REQUIRE(c.count_gates(OpType::Unitary3qBox) == 4);
    }
  }

  GIVEN("A projector with rank > 2 ^ n-1") {
    WHEN("The projector is 2q") {
      Eigen::MatrixXcd P = Eigen::MatrixXcd::Zero(4, 4);
      P(0, 0) = 1;
      P(1, 1) = 1;
      P(2, 2) = 1;
      auto [c, _] = projector_assertion_synthesis(P);
      REQUIRE(c.n_qubits() == 3);
      REQUIRE(c.count_gates(OpType::Unitary3qBox) == 4);
    }
    WHEN("The projector is 3q") {
      Eigen::MatrixXcd P = Eigen::MatrixXcd::Zero(8, 8);
      P(0, 0) = 1;
      P(1, 1) = 1;
      P(2, 2) = 1;
      P(3, 3) = 1;
      P(4, 4) = 1;
      REQUIRE_THROWS_AS(projector_assertion_synthesis(P), CircuitInvalidity);
    }
  }
}

SCENARIO(
    "Testing adding a projector based assertion and decomposing the circuit") {
  GIVEN("A 2q projector") {
    Circuit circ(2);
    Eigen::MatrixXcd bell(4, 4);
    bell << 0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0.5;
    ProjectorAssertionBox box(bell);
    circ.add_assertion(
        box, {Qubit(0), Qubit(1)}, std::nullopt, "bell projector");
    circ.add_assertion(
        box, {Qubit(1), Qubit(0)}, std::nullopt, "bell projector");
    circ.add_assertion(box, {Qubit(0), Qubit(1)});
    circ.add_assertion(box, {Qubit(1), Qubit(0)});
    REQUIRE(
        circ.get_reg_info(
            c_debug_zero_prefix() + "_" + c_debug_default_name()) !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(
            c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_zero_prefix() + "_bell projector") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_zero_prefix() + "_bell projector(1)") !=
        std::nullopt);
    CompilationUnit cu(circ);
    REQUIRE(DecomposeBoxes()->apply(cu));
  }
  GIVEN("A 3q projector") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 1.5, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Eigen::MatrixXcd P = Eigen::MatrixXcd::Zero(8, 8);
    P(0, 0) = 1;
    P(1, 1) = 1;
    P(7, 7) = 1;
    ProjectorAssertionBox box(P);
    circ.add_assertion(
        box, {Qubit(0), Qubit(1), Qubit(2)}, std::nullopt, "random projector");
    circ.add_assertion(
        box, {Qubit(1), Qubit(0), Qubit(2)}, std::nullopt, "random projector");
    circ.add_assertion(box, {Qubit(0), Qubit(1), Qubit(2)});
    circ.add_assertion(box, {Qubit(1), Qubit(0), Qubit(2)});
    REQUIRE(
        circ.get_reg_info(
            c_debug_zero_prefix() + "_" + c_debug_default_name()) !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(
            c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_zero_prefix() + "_random projector") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_zero_prefix() + "_random projector(1)") !=
        std::nullopt);
    CompilationUnit cu(circ);
    REQUIRE(DecomposeBoxes()->apply(cu));
  }
}

SCENARIO("Testing stabiliser based assertion") {
  GIVEN("Random stabilisers") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 1.5, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliStabiliser pauli1 = {{Pauli::X, Pauli::X}, true};
    PauliStabiliser pauli2 = {{Pauli::Z, Pauli::Z}, true};
    PauliStabiliser pauli3 = {{Pauli::Z, Pauli::Z}, false};
    PauliStabiliserList stabilisers = {pauli1, pauli2, pauli3};
    StabiliserAssertionBox box(stabilisers);
    circ.add_assertion(
        box, {Qubit(0), Qubit(2)}, Qubit(1), "random stabiliser");
    circ.add_assertion(
        box, {Qubit(0), Qubit(2)}, Qubit(1), "random stabiliser");
    circ.add_assertion(box, {Qubit(0), Qubit(2)}, Qubit(1));
    circ.add_assertion(box, {Qubit(0), Qubit(2)}, Qubit(1));
    REQUIRE(
        circ.get_reg_info(
            c_debug_zero_prefix() + "_" + c_debug_default_name()) !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(
            c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_zero_prefix() + "_random stabiliser") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_zero_prefix() + "_random stabiliser(1)") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(
            c_debug_one_prefix() + "_" + c_debug_default_name()) !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(
            c_debug_one_prefix() + "_" + c_debug_default_name() + "(1)") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_one_prefix() + "_random stabiliser") !=
        std::nullopt);
    REQUIRE(
        circ.get_reg_info(c_debug_one_prefix() + "_random stabiliser(1)") !=
        std::nullopt);
    CompilationUnit cu(circ);
    REQUIRE(DecomposeBoxes()->apply(cu));
  }

  GIVEN("Invalid input") {
    WHEN("Empty input") {
      PauliStabiliserList stabilisers = {};
      REQUIRE_THROWS_AS(StabiliserAssertionBox(stabilisers), CircuitInvalidity);
    }
    WHEN("Unequal lengths") {
      PauliStabiliser pauli1 = {{Pauli::X}, true};
      PauliStabiliser pauli2 = {{Pauli::Z, Pauli::Z}, true};
      PauliStabiliserList stabilisers = {pauli1, pauli2};
      REQUIRE_THROWS_AS(StabiliserAssertionBox(stabilisers), CircuitInvalidity);
    }
    WHEN("Identity") {
      REQUIRE_THROWS_AS(
          PauliStabiliser({Pauli::I, Pauli::I, Pauli::I}, true), NotValid);
      REQUIRE_THROWS_AS(
          PauliStabiliser({Pauli::I, Pauli::I, Pauli::I}, false), NotValid);
    }
  }
}

SCENARIO("Testing stibiliser based assertion serialization") {
  GIVEN("Serialise a stabiliser box") {
    PauliStabiliser pauli1 = {{Pauli::X, Pauli::X}, true};
    PauliStabiliser pauli2 = {{Pauli::Z, Pauli::Z}, true};
    nlohmann::json j_pauli1 = pauli1;
    PauliStabiliser new_pauli1 = j_pauli1.get<PauliStabiliser>();
    REQUIRE(new_pauli1 == pauli1);
    PauliStabiliserList bell = {pauli1, pauli2};
    nlohmann::json j_bell = bell;
    PauliStabiliserList new_bell = j_bell.get<PauliStabiliserList>();
    REQUIRE(new_bell == bell);
    StabiliserAssertionBox bell_box(new_bell);
    Circuit circ(3);
    circ.add_assertion(
        bell_box, {Qubit(0), Qubit(2)}, Qubit(1), "bell stabiliser");
    nlohmann::json j_box = circ;
    const Circuit new_c = j_box.get<Circuit>();
    const auto& new_box = static_cast<const StabiliserAssertionBox&>(
        *new_c.get_commands()[0].get_op_ptr());
    REQUIRE(bell_box.get_stabilisers() == new_box.get_stabilisers());
  }
}
}  // namespace tket
