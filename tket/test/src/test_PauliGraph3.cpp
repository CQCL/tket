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
#include <fstream>

#include "testutil.hpp"
#include "tket/Converters3/Converters.hpp"
#include "tket/PauliGraph3/PauliGraph.hpp"

namespace tket {
namespace test_PauliGraph3 {

using namespace pg;

SCENARIO("Correct creation of refactored PauliGraphs") {
  GIVEN("A Clifford circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Vdg, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 1qb circuit") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with no interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rx, 1.6, {1});
    circ.add_op<unsigned>(OpType::Rz, 1.3, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with some anti-commuting interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {1});
    circ.add_op<unsigned>(OpType::XXPhase, 1.1, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with some commuting interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.1, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit a Clifford-angled ZZPhase") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 0.5, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 1qb circuit with stuff to merge") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with stuff to merge") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.1, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.8, {0});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.6, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A circuit with Cliffords and non-Cliffords") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {0});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.8, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A dense example") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {1});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {1});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {3});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A more interesting example (tof_3)") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::T, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Tdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::Tdg, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::Tdg, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::T, {2});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::T, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Tdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A circuit with a PauliExpBox") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::ZZPhase, 0.2, {0, 1});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    PauliExpBox peb({Pauli::Y, Pauli::X}, 0.333);
    circ.add_box(peb, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    THEN("Print diagram to file") {
      std::ofstream dot_file("pauligraph.dot");
      pg.to_graphviz(dot_file);
      dot_file.close();
      remove("pauligraph.dot");
    }
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("Teleportation") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_measure(0, 0);
    circ.add_measure(1, 1);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{2}, {1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Z, {}, uvec{2}, {0}, 1);
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
  }
}

}  // namespace test_PauliGraph3
}  // namespace tket
