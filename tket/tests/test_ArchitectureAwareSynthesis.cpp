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

#include "Architecture/Architecture.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "testutil.hpp"

namespace tket {
using Connection = Architecture::Connection;
SCENARIO("Routing of aas example") {
  GIVEN("aas routing - simple example") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example II") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example III") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example IV") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example V") {
    Architecture arc(std::vector<Connection>{{Node(0), Node(1)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example VI") {
    Architecture arc(std::vector<Connection>{{Node(0), Node(2)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));

    const auto s = tket_sim::get_unitary(circ);
    const auto s1 = tket_sim::get_unitary(result);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        s, s1, tket_sim::MatrixEquivalence::EQUAL));
  }
  GIVEN("aas routing - simple example VII") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(2)}, {Node(2), Node(4)}, {Node(4), Node(6)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));

    const auto s = tket_sim::get_unitary(circ);
    const auto s1 = tket_sim::get_unitary(result);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        s, s1, tket_sim::MatrixEquivalence::EQUAL));
  }
  GIVEN("aas routing - simple example VIII") {
    Architecture arc(std::vector<Connection>{
        {Node(1000), Node(10)}, {Node(10), Node(100)}, {Node(100), Node(1)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example IX, other gate set") {
    Architecture arc(std::vector<Connection>{
        {Node(1000), Node(10)}, {Node(10), Node(100)}, {Node(100), Node(1)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::X, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::X, {3});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing with measure") {
    Architecture arc(std::vector<Connection>{{Node(0), Node(2)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(2, 2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    for (unsigned mes = 0; mes < 2; ++mes) {
      circ.add_measure(mes, mes);
    }

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
  }
  GIVEN("aas routing - circuit with fewer qubits then nodes in the arch") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - circuit with fewer qubits then nodes in the arch II") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)},
        {Node(1), Node(2)},
        {Node(2), Node(3)},
        {Node(3), Node(4)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
}
}  // namespace tket
