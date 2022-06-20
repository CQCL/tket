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
#include <iostream>

#include "CircuitsForTesting.hpp"
#include "Converters/PauliGadget.hpp"
#include "PauliGraph/ConjugatePauliFunctions.hpp"
#include "Utils/PauliStrings.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_PauliString {

SCENARIO("Testing equality of QubitPauliTensor") {
  Qubit q0 = Qubit("q", 0);
  Qubit q1 = Qubit("q", 1);
  Qubit q2 = Qubit("r", 0);
  Qubit q3 = Qubit("s");
  Qubit q4 = Qubit("t", 0, 1);
  Qubit q5 = Qubit("t", 0, 0);
  GIVEN("Two exactly identical Pauli strings") {
    QubitPauliMap map = {
        {q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}};
    QubitPauliTensor a(map, i_);
    QubitPauliTensor b(map, i_);
    REQUIRE(a == b);
    THEN("We add some extra Is on each one") {
      a.string.map.insert({q4, Pauli::I});
      b.string.map.insert({q5, Pauli::I});
      REQUIRE(a == b);
    }
  }
  GIVEN("Two Pauli strings with different Paulis but same coefficient") {
    QubitPauliTensor a(q0, Pauli::X);
    QubitPauliTensor b(q0, Pauli::Y);
    REQUIRE(!(a == b));
  }
  GIVEN("Two Pauli strings with disjoint Paulis but same coefficient") {
    QubitPauliTensor a(q0, Pauli::X);
    QubitPauliTensor b(q1, Pauli::X);
    REQUIRE(!(a == b));
  }
  GIVEN("Two Pauli strings with same Paulis but different coefficient") {
    QubitPauliTensor a(q0, Pauli::X, 1.);
    QubitPauliTensor b(q0, Pauli::X, i_);
    REQUIRE(!(a == b));
  }
  GIVEN("Two completely different Pauli strings") {
    QubitPauliString tensor_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliString tensor_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    QubitPauliTensor a(tensor_a, 1.);
    QubitPauliTensor b(tensor_b, i_);
    REQUIRE(!(a == b));
  }
}

SCENARIO("Testing multiplication of QubitPauliTensor") {
  Qubit q0 = Qubit("q", 0);
  Qubit q1 = Qubit("q", 1);
  Qubit q2 = Qubit("r", 0);
  Qubit q3 = Qubit("s");
  Qubit q4 = Qubit("t", 0, 1);
  GIVEN("Two Pauli strings with disjoint non-trivial components") {
    QubitPauliTensor a(q0, Pauli::X, 2.);
    QubitPauliTensor b(q1, Pauli::Y, i_);
    QubitPauliTensor c({{q0, Pauli::X}, {q1, Pauli::Y}}, 2. * i_);
    REQUIRE((a * b) == c);
  }
  GIVEN("Multiplying by a trivial Pauli string") {
    QubitPauliTensor a(q0, Pauli::X, 2.);
    QubitPauliTensor b(q0, Pauli::X, 3. * i_);
    REQUIRE((a * QubitPauliTensor(1.5 * i_)) == b);
  }
  GIVEN("Two exactly identical Pauli strings") {
    QubitPauliMap map = {
        {q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}};
    QubitPauliTensor a(map, i_);
    QubitPauliTensor b(-1.);
    REQUIRE((a * a) == b);
  }
  GIVEN("Each individual Pauli combination") {
    QubitPauliTensor I(q0, Pauli::I);
    QubitPauliTensor X(q0, Pauli::X);
    QubitPauliTensor Y(q0, Pauli::Y);
    QubitPauliTensor Z(q0, Pauli::Z);
    QubitPauliTensor i(i_);
    QubitPauliTensor mi(-i_);
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
    QubitPauliString tensor_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliString tensor_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    QubitPauliTensor a(tensor_a, 2.);
    QubitPauliTensor b(tensor_b, -1.5 * i_);
    QubitPauliString tensor_c(
        {{q0, Pauli::X},
         {q1, Pauli::X},
         {q2, Pauli::X},
         {q3, Pauli::Z},
         {q4, Pauli::Y}});
    QubitPauliTensor c(tensor_c, 3.);
    REQUIRE((a * b) == c);
  }
}

SCENARIO("Test basic conjugations") {
  const Qubit q0 = Qubit("q", 0);
  const Qubit q1 = Qubit("q", 1);
  Circuit circ;
  circ.add_qubit(q0);
  circ.add_qubit(q1);
  // add some arbitrary rotations to get away from |00> state
  const auto& prepend = CircuitsForTesting::get().prepend_2qb_circuit;
  const double angle = 0.845;
  // generate all different 2-qb pauli strings
  std::vector<QubitPauliTensor> qps_vec;
  for (const auto& map_entry : QubitPauliTensor::get_mult_matrix()) {
    const std::pair<Pauli, Pauli>& paulis = map_entry.first;
    QubitPauliMap map = {{q0, paulis.first}, {q1, paulis.second}};
    QubitPauliTensor qps(map);
    qps_vec.push_back(qps);
  }

  const auto perform_test = [&q0, &qps_vec, angle, &prepend](
                                OpType op_type, OpType op_type_dag,
                                OpType tensor_op_type, bool reverse = false) {
    for (QubitPauliTensor qps : qps_vec) {
      Circuit test(2);
      test.add_op<unsigned>(op_type, {0});
      append_single_pauli_gadget(test, qps, angle);
      test.add_op<unsigned>(op_type_dag, {0});
      test = prepend >> test;
      conjugate_PauliTensor(qps, tensor_op_type, q0, reverse);
      Circuit test1 = prepend;
      append_single_pauli_gadget(test1, qps, angle);
      REQUIRE(test_statevector_comparison(test, test1));
    }
  };

  WHEN("Test Hs") { perform_test(OpType::H, OpType::H, OpType::H); }
  WHEN("Test Ss") {
    perform_test(OpType::S, OpType::Sdg, OpType::S);
    perform_test(OpType::Sdg, OpType::S, OpType::S, true);
  }
  WHEN("Test Vs") {
    perform_test(OpType::V, OpType::Vdg, OpType::V);
    perform_test(OpType::Vdg, OpType::V, OpType::V, true);
  }
  WHEN("Test Xs") { perform_test(OpType::X, OpType::X, OpType::X); }
  WHEN("Test Zs") { perform_test(OpType::Z, OpType::Z, OpType::Z); }
  WHEN("Test CXs") {
    for (QubitPauliTensor qps : qps_vec) {
      Circuit test(2);
      test.add_op<unsigned>(OpType::CX, {0, 1});
      append_single_pauli_gadget(test, qps, angle);
      test.add_op<unsigned>(OpType::CX, {0, 1});
      test = prepend >> test;
      conjugate_PauliTensor(qps, OpType::CX, q0, q1);
      Circuit test2 = prepend;
      append_single_pauli_gadget(test2, qps, angle);
      REQUIRE(test_statevector_comparison(test, test2));
    }
  }
}

SCENARIO("Test hashing") {
  GIVEN("Trivial strings") {
    QubitPauliString qps1;
    QubitPauliString qps2;
    REQUIRE(hash_value(qps1) == hash_value(qps2));
    WHEN("Add I Pauli") {
      qps1.map[Qubit(0)] = Pauli::I;
      REQUIRE(hash_value(qps1) == hash_value(qps2));
    }
  }
  GIVEN("Nontrivial strings") {
    QubitPauliMap qpm1{
        {Qubit(0), Pauli::Z},
        {Qubit(1), Pauli::Y},
        {Qubit(2), Pauli::X},
        {Qubit(3), Pauli::I}};
    QubitPauliString qps1(qpm1);
    QubitPauliString qps2(qpm1);
    qps1.map[Qubit(4)] = Pauli::X;
    qps2.map[Qubit(4)] = Pauli::X;
    qps2.map[Qubit(5)] = Pauli::I;
    REQUIRE(hash_value(qps1) == hash_value(qps2));
  }
  GIVEN("Trivial tensor") {
    QubitPauliTensor qpt1;
    QubitPauliTensor qpt2;
    REQUIRE(hash_value(qpt1) == hash_value(qpt2));
    WHEN("Add I Pauli") {
      qpt1.string.map[Qubit(0)] = Pauli::I;
      REQUIRE(hash_value(qpt1) == hash_value(qpt2));
    }
  }
  GIVEN("Nontrivial tensors") {
    QubitPauliMap qpm1{
        {Qubit(0), Pauli::Z},
        {Qubit(1), Pauli::Y},
        {Qubit(2), Pauli::X},
        {Qubit(3), Pauli::I}};
    QubitPauliString qps1(qpm1);
    QubitPauliString qps2(qpm1);
    qps1.map[Qubit(4)] = Pauli::X;
    qps2.map[Qubit(4)] = Pauli::X;
    qps2.map[Qubit(5)] = Pauli::I;
    QubitPauliTensor qpt1(qps1, .5 * i_);
    QubitPauliTensor qpt2(qps2, .5 * i_);
    qpt2.string.map[Qubit(6)] = Pauli::I;
    REQUIRE(hash_value(qpt1) == hash_value(qpt2));
  }
}

SCENARIO("Test matrix product utilities") {
  GIVEN("Simple operator and its +1 eigenvector") {
    const QubitPauliString op({{Qubit(0), Pauli::X}, {Qubit(1), Pauli::Y}});
    Eigen::VectorXcd state(4);

    state << Complex(0.5, 0), Complex(0, 0.5), Complex(0.5, 0), Complex(0, 0.5);
    WHEN("Operator acts on statevector") {
      Eigen::VectorXcd dotproduct = op.dot_state(state);
      THEN("The eigenvector is unchanged") {
        REQUIRE(dotproduct.isApprox(state));
      }
    }

    WHEN("Expectation value is calculated w.r.t. +1 eigenvector") {
      Complex eigenval = op.state_expectation(state);
      THEN("The eigenvalue is correct") { REQUIRE(eigenval == Complex(1, 0)); }
    }
  }
  GIVEN("A qubit list with repeats") {
    const QubitPauliString op({{Qubit(0), Pauli::X}, {Qubit(1), Pauli::Y}});
    REQUIRE_THROWS(op.to_sparse_matrix({Qubit(0), Qubit(0), Qubit(1)}));
  }
  GIVEN("A qubit list that doesn't contain all qubits in the string") {
    const QubitPauliString op({{Qubit(0), Pauli::X}, {Qubit(1), Pauli::Y}});
    REQUIRE_THROWS(op.to_sparse_matrix({Qubit(0), Qubit(2)}));
  }
}

}  // namespace test_PauliString
}  // namespace tket
