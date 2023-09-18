// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Utils/PauliStrings2.hpp"

namespace tket {
namespace test_PauliString {

SCENARIO("Testing equality of sparse PauliTensor variants") {
  Qubit q0 = Qubit("q", 0);
  Qubit q1 = Qubit("q", 1);
  Qubit q2 = Qubit("r", 0);
  Qubit q3 = Qubit("s");
  Qubit q4 = Qubit("t", 0, 1);
  Qubit q5 = Qubit("t", 0, 0);
  GIVEN("Two exactly identical Pauli strings") {
    QubitPauliMap map = {
        {q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}};
    SpPauliTensor a(map, i_);
    SpPauliTensor b(map, i_);
    REQUIRE(a == b);
    THEN("We add some extra Is on each one") {
      a.set(q4, Pauli::I);
      b.set(q5, Pauli::I);
      REQUIRE(a == b);
    }
  }
  GIVEN("Two Pauli strings with different Paulis but same coefficient") {
    SpPauliString a({{q0, Pauli::X}});
    SpPauliString b({{q0, Pauli::Y}});
    REQUIRE(a != b);
  }
  GIVEN("Two Pauli strings with disjoint Paulis but same coefficient") {
    SpPauliString a(q0, Pauli::X);
    SpPauliString b(q1, Pauli::X);
    REQUIRE(a != b);
  }
  GIVEN("Two Pauli strings with same Paulis but different coefficient") {
    SpPauliTensor a(q0, Pauli::X, 1.);
    SpPauliTensor b(q0, Pauli::X, i_);
    REQUIRE(a != b);
  }
  GIVEN("Two completely different Pauli strings") {
    QubitPauliMap qpm_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliMap qpm_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    SpPauliTensor a(qpm_a, 1.);
    SpPauliTensor b(qpm_b, i_);
    REQUIRE(a != b);
  }
}

}  // namespace test_PauliString
}  // namespace tket
