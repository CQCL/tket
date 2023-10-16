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
    SpCxPauliTensor a(map, i_);
    SpCxPauliTensor b(map, i_);
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
    REQUIRE(a < b);
  }
  GIVEN("Two Pauli strings with disjoint Paulis but same coefficient") {
    SpPauliString a(q0, Pauli::X);
    SpPauliString b(q1, Pauli::X);
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two Pauli strings with same Paulis but different coefficient") {
    SpCxPauliTensor a(q0, Pauli::X, 1.);
    SpCxPauliTensor b(q0, Pauli::X, i_);
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two completely different Pauli strings") {
    QubitPauliMap qpm_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliMap qpm_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    SpCxPauliTensor a(qpm_a, 1.);
    SpCxPauliTensor b(qpm_b, i_);
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
}

SCENARIO("Testing equality of dense PauliTensor variants") {
  GIVEN("Two exactly identical Pauli strings") {
    DensePauliMap map = {Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
    CxPauliTensor a(map, i_);
    CxPauliTensor b(map, i_);
    REQUIRE(a == b);
    THEN("We add some extra Is on each one") {
      a.set(4, Pauli::I);
      b.set(5, Pauli::I);
      REQUIRE(a == b);
    }
  }
  GIVEN("Two Pauli strings with different Paulis but same coefficient") {
    PauliString a({Pauli::X});
    PauliString b({Pauli::Y});
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
  GIVEN("Two Pauli strings with disjoint Paulis but same coefficient") {
    PauliString a({Pauli::X});
    PauliString b({Pauli::I, Pauli::X});
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two Pauli strings with same Paulis but different coefficient") {
    CxPauliTensor a({Pauli::X}, 1.);
    CxPauliTensor b({Pauli::X}, i_);
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two completely different Pauli strings") {
    DensePauliMap qpm_a({Pauli::I, Pauli::X, Pauli::Y, Pauli::Z});
    DensePauliMap qpm_b({Pauli::X, Pauli::I, Pauli::Z, Pauli::Y});
    CxPauliTensor a(qpm_a, 1.);
    CxPauliTensor b(qpm_b, i_);
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
}

SCENARIO("Testing multiplication of sparse PauliTensor") {
  Qubit q0 = Qubit("q", 0);
  Qubit q1 = Qubit("q", 1);
  Qubit q2 = Qubit("r", 0);
  Qubit q3 = Qubit("s");
  Qubit q4 = Qubit("t", 0, 1);
  GIVEN("Two Pauli strings with disjoint non-trivial components") {
    SpPauliString a(q0, Pauli::X);
    SpPauliString b(q1, Pauli::Y);
    SpPauliString c({{q0, Pauli::X}, {q1, Pauli::Y}});
    REQUIRE((a * b) == c);
  }
  GIVEN("Multiplying by a trivial Pauli string") {
    SpCxPauliTensor a(q0, Pauli::X, 2.);
    SpCxPauliTensor b(q0, Pauli::X, 3. * i_);
    REQUIRE((a * SpCxPauliTensor({}, 1.5 * i_)) == b);
  }
  GIVEN("Two exactly identical Pauli strings") {
    QubitPauliMap map = {
        {q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}};
    SpPauliStabiliser a(map, 3);
    SpPauliStabiliser b({}, 2);
    REQUIRE((a * a) == b);
  }
  GIVEN("Each individual Pauli combination") {
    SpPauliStabiliser I(q0, Pauli::I);
    SpPauliStabiliser X(q0, Pauli::X);
    SpPauliStabiliser Y(q0, Pauli::Y);
    SpPauliStabiliser Z(q0, Pauli::Z);
    SpPauliStabiliser i({}, 1);
    SpPauliStabiliser mi({}, 3);
    REQUIRE((I * I) == I);
    REQUIRE((I * X) == X);
    REQUIRE((I * Y) == Y);
    REQUIRE((I * Z) == Z);
    REQUIRE((X * I) == X);
    REQUIRE((X * X) == I);
    REQUIRE((X * Y) == (i * Z));
    REQUIRE((X * Z) == (mi * Y));
    REQUIRE((Y * I) == Y);
    REQUIRE((Y * X) == (mi * Z));
    REQUIRE((Y * Y) == I);
    REQUIRE((Y * Z) == (i * X));
    REQUIRE((Z * I) == Z);
    REQUIRE((Z * X) == (i * Y));
    REQUIRE((Z * Y) == (mi * X));
    REQUIRE((Z * Z) == I);
  }
  GIVEN("2*IXYZ(I) * -1.5i*XIZ(I)Y") {
    QubitPauliMap tensor_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliMap tensor_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    SpCxPauliTensor a(tensor_a, 2.);
    SpCxPauliTensor b(tensor_b, -1.5 * i_);
    QubitPauliMap tensor_c(
        {{q0, Pauli::X},
         {q1, Pauli::X},
         {q2, Pauli::X},
         {q3, Pauli::Z},
         {q4, Pauli::Y}});
    SpCxPauliTensor c(tensor_c, 3.);
    REQUIRE((a * b) == c);
  }
}

SCENARIO("Testing multiplication of dense PauliTensor") {
  GIVEN("Two Pauli strings with disjoint non-trivial components") {
    PauliString a({Pauli::X});
    PauliString b({Pauli::I, Pauli::Y});
    PauliString c({Pauli::X, Pauli::Y});
    REQUIRE((a * b) == c);
  }
  GIVEN("Multiplying by a trivial Pauli string") {
    CxPauliTensor a({Pauli::X}, 2.);
    CxPauliTensor b({Pauli::X}, 3. * i_);
    REQUIRE((a * CxPauliTensor({}, 1.5 * i_)) == b);
  }
  GIVEN("Two exactly identical Pauli strings") {
    DensePauliMap map = {Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
    PauliStabiliser a(map, 3);
    PauliStabiliser b({}, 2);
    REQUIRE((a * a) == b);
  }
  GIVEN("Each individual Pauli combination") {
    PauliStabiliser I({Pauli::I});
    PauliStabiliser X({Pauli::X});
    PauliStabiliser Y({Pauli::Y});
    PauliStabiliser Z({Pauli::Z});
    PauliStabiliser i({}, 1);
    PauliStabiliser mi({}, 3);
    REQUIRE((I * I) == I);
    REQUIRE((I * X) == X);
    REQUIRE((I * Y) == Y);
    REQUIRE((I * Z) == Z);
    REQUIRE((X * I) == X);
    REQUIRE((X * X) == I);
    REQUIRE((X * Y) == (i * Z));
    REQUIRE((X * Z) == (mi * Y));
    REQUIRE((Y * I) == Y);
    REQUIRE((Y * X) == (mi * Z));
    REQUIRE((Y * Y) == I);
    REQUIRE((Y * Z) == (i * X));
    REQUIRE((Z * I) == Z);
    REQUIRE((Z * X) == (i * Y));
    REQUIRE((Z * Y) == (mi * X));
    REQUIRE((Z * Z) == I);
  }
  GIVEN("2*IXYZ(I) * -1.5i*XIZ(I)Y") {
    DensePauliMap tensor_a({Pauli::I, Pauli::X, Pauli::Y, Pauli::Z});
    DensePauliMap tensor_b({Pauli::X, Pauli::I, Pauli::Z, Pauli::I, Pauli::Y});
    CxPauliTensor a(tensor_a, 2.);
    CxPauliTensor b(tensor_b, -1.5 * i_);
    DensePauliMap tensor_c({Pauli::X, Pauli::X, Pauli::X, Pauli::Z, Pauli::Y});
    CxPauliTensor c(tensor_c, 3.);
    REQUIRE((a * b) == c);
  }
}

SCENARIO("Test hashing for sparse PauliTensor") {
  GIVEN("Trivial strings") {
    SpPauliString qps1;
    SpPauliString qps2;
    REQUIRE(qps1.hash_value() == qps2.hash_value());
    WHEN("Add I Pauli") {
      qps1.set(Qubit(0), Pauli::I);
      REQUIRE(qps1.hash_value() == qps2.hash_value());
    }
  }
  GIVEN("Nontrivial strings") {
    QubitPauliMap qpm{
        {Qubit(0), Pauli::Z},
        {Qubit(1), Pauli::Y},
        {Qubit(2), Pauli::X},
        {Qubit(3), Pauli::I}};
    SpPauliString qps1(qpm);
    SpPauliString qps2(qpm);
    qps1.set(Qubit(4), Pauli::X);
    qps2.set(Qubit(4), Pauli::X);
    qps2.set(Qubit(5), Pauli::I);
    REQUIRE(qps1.hash_value() == qps2.hash_value());
  }
  GIVEN("Trivial tensor") {
    SpCxPauliTensor qpt1;
    SpCxPauliTensor qpt2;
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    WHEN("Add I Pauli") {
      qpt1.set(Qubit(0), Pauli::I);
      REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    }
  }
  GIVEN("Nontrivial tensors") {
    QubitPauliMap qpm{
        {Qubit(0), Pauli::Z},
        {Qubit(1), Pauli::Y},
        {Qubit(2), Pauli::X},
        {Qubit(3), Pauli::I}};
    SpSymPauliTensor qpt1(qpm, .5 * i_);
    SpSymPauliTensor qpt2(qpm, .5 * i_);
    qpt1.set(Qubit(4), Pauli::X);
    qpt2.set(Qubit(4), Pauli::X);
    qpt2.set(Qubit(5), Pauli::I);
    qpt2.set(Qubit(6), Pauli::I);
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
  }
}

SCENARIO("Test hashing for dense PauliTensor") {
  GIVEN("Trivial strings") {
    PauliString qps1;
    PauliString qps2;
    REQUIRE(qps1.hash_value() == qps2.hash_value());
    WHEN("Add I Pauli") {
      qps1.set(0, Pauli::I);
      REQUIRE(qps1.hash_value() == qps2.hash_value());
    }
  }
  GIVEN("Nontrivial strings") {
    DensePauliMap qpm{Pauli::Z, Pauli::Y, Pauli::X, Pauli::I};
    PauliString qps1(qpm);
    PauliString qps2(qpm);
    qps1.set(4, Pauli::X);
    qps2.set(4, Pauli::X);
    qps2.set(5, Pauli::I);
    REQUIRE(qps1.hash_value() == qps2.hash_value());
  }
  GIVEN("Trivial tensor") {
    CxPauliTensor qpt1;
    CxPauliTensor qpt2;
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    WHEN("Add I Pauli") {
      qpt1.set(0, Pauli::I);
      REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    }
  }
  GIVEN("Nontrivial tensors") {
    DensePauliMap qpm{Pauli::Z, Pauli::Y, Pauli::X, Pauli::I};
    SymPauliTensor qpt1(qpm, .5 * i_);
    SymPauliTensor qpt2(qpm, .5 * i_);
    qpt1.set(4, Pauli::X);
    qpt2.set(4, Pauli::X);
    qpt2.set(5, Pauli::I);
    qpt2.set(6, Pauli::I);
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
  }
}

}  // namespace test_PauliString
}  // namespace tket
