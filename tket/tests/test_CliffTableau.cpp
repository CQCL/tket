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

#include "Clifford/CliffTableau.hpp"
#include "Converters/Converters.hpp"

namespace tket {
namespace test_CliffTableau {

static void add_ops_list_one_to_circuit(Circuit& circ) {
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  circ.add_op<unsigned>(OpType::S, {1});
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  circ.add_op<unsigned>(OpType::CX, {1, 2});
  circ.add_op<unsigned>(OpType::Vdg, {1});
  circ.add_op<unsigned>(OpType::CX, {1, 0});
}
static Circuit get_test_circ() {
  Circuit circ(3);
  add_ops_list_one_to_circuit(circ);
  return circ;
}
static void add_ops_list_two_to_circuit(
    Circuit& circ, OpType middle_op = OpType::S) {
  circ.add_op<unsigned>(OpType::H, {0});
  circ.add_op<unsigned>(OpType::V, {1});
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  circ.add_op<unsigned>(OpType::CX, {1, 2});
  circ.add_op<unsigned>(middle_op, {2});
  circ.add_op<unsigned>(OpType::CX, {1, 2});
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  circ.add_op<unsigned>(OpType::H, {0});
  circ.add_op<unsigned>(OpType::Vdg, {1});
}
static CliffTableau get_tableau_with_gates_applied_at_front() {
  CliffTableau tab(3);
  tab.apply_gate_at_front(OpType::CX, std::vector<unsigned>({1, 0}));
  tab.apply_gate_at_front(OpType::Vdg, std::vector<unsigned>({1}));
  tab.apply_gate_at_front(OpType::CX, std::vector<unsigned>({1, 2}));
  tab.apply_gate_at_front(OpType::CX, std::vector<unsigned>({0, 1}));
  tab.apply_gate_at_front(OpType::S, std::vector<unsigned>({1}));
  tab.apply_gate_at_front(OpType::CX, std::vector<unsigned>({0, 1}));
  return tab;
}

SCENARIO("Correct creation of Tableaus") {
  GIVEN("An identity circuit") {
    Circuit circ(3);
    CliffTableau tab = circuit_to_tableau(circ);
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Z, 1.));
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 1)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::Z, 1.));
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 2)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 2), Pauli::Z, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 1)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::X, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 2)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 2), Pauli::X, 1.));
  }
  GIVEN("A single S gate") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::S, {0});
    CliffTableau tab = circuit_to_tableau(circ);
    CliffTableau s_tab(3);
    s_tab.apply_S_at_front(0);
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Z, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Y, -1.));
    REQUIRE(tab == s_tab);
  }
  GIVEN("A single V gate") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::V, {0});
    CliffTableau tab = circuit_to_tableau(circ);
    CliffTableau v_tab(3);
    v_tab.apply_V_at_front(0);
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Y, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X, 1.));
    REQUIRE(tab == v_tab);
  }
  GIVEN("A single H gate") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    CliffTableau tab = circuit_to_tableau(circ);
    CliffTableau h_tab(3);
    h_tab.apply_gate_at_front(OpType::H, std::vector<unsigned>({0}));
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Z, 1.));
    REQUIRE(tab == h_tab);
  }
  GIVEN("A single CX gate") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    CliffTableau tab = circuit_to_tableau(circ);
    CliffTableau cx_tab(3);
    cx_tab.apply_CX_at_front(0, 1);
    REQUIRE(
        tab.get_zpauli(Qubit(q_default_reg(), 0)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Z, 1.));
    REQUIRE(
        tab.get_xpauli(Qubit(q_default_reg(), 1)) ==
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::X, 1.));
    QubitPauliTensor correct =
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::Z) *
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::Z);
    REQUIRE(tab.get_zpauli(Qubit(q_default_reg(), 1)) == correct);
    correct = QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X) *
              QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::X);
    REQUIRE(tab.get_xpauli(Qubit(q_default_reg(), 0)) == correct);
    REQUIRE(tab == cx_tab);
  }
  GIVEN("A Clifford circuit") {
    const auto circ = get_test_circ();
    CliffTableau tab = circuit_to_tableau(circ);
    const auto rev_tab = get_tableau_with_gates_applied_at_front();
    REQUIRE(tab == rev_tab);
  }
  GIVEN("A PI/2 rotation") {
    auto circ = get_test_circ();
    CliffTableau tab = circuit_to_tableau(circ);
    QubitPauliTensor pauli =
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X) *
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::Y) *
        QubitPauliTensor(Qubit(q_default_reg(), 2), Pauli::Z);
    tab.apply_pauli_at_end(pauli, 3);

    add_ops_list_two_to_circuit(circ, OpType::Sdg);
    CliffTableau correct_tab = circuit_to_tableau(circ);
    REQUIRE(tab == correct_tab);
  }
  GIVEN("A PI/2 rotation at front") {
    auto tab = get_tableau_with_gates_applied_at_front();
    QubitPauliTensor pauli =
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X) *
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::Y) *
        QubitPauliTensor(Qubit(q_default_reg(), 2), Pauli::Z);
    tab.apply_pauli_at_front(pauli, 1);

    Circuit circ(3);
    add_ops_list_two_to_circuit(circ);
    add_ops_list_one_to_circuit(circ);
    CliffTableau correct_tab = circuit_to_tableau(circ);
    REQUIRE(tab == correct_tab);
  }
  GIVEN("Combining two circuits via tableau compose") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    CliffTableau first = circuit_to_tableau(circ);

    Circuit circ1(3);
    add_ops_list_two_to_circuit(circ1);
    CliffTableau second = circuit_to_tableau(circ1);
    CliffTableau correct = circuit_to_tableau(circ >> circ1);
    CliffTableau result = CliffTableau::compose(first, second);
    REQUIRE(result == correct);
  }
}

SCENARIO("Error handling in Tableau generation") {
  GIVEN("Add a non-clifford gate at end") {
    CliffTableau tab(2);
    REQUIRE_THROWS_AS(
        tab.apply_gate_at_end(OpType::T, std::vector<unsigned>({0})), NotValid);
  }
  GIVEN("Add a non-clifford gate at front") {
    CliffTableau tab(2);
    REQUIRE_THROWS_AS(
        tab.apply_gate_at_front(OpType::Tdg, std::vector<unsigned>({0})),
        NotValid);
  }
  GIVEN("Tableau from a non-Clifford circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CH, {1, 0});
    REQUIRE_THROWS_AS(circuit_to_tableau(circ), NotValid);
  }
}

SCENARIO("Synthesis of circuits from Tableaus") {
  GIVEN("An identity circuit") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    CliffTableau tab = circuit_to_tableau(circ);
    Circuit res = tableau_to_circuit(tab);
    CliffTableau res_tab = circuit_to_tableau(res);
    REQUIRE(res_tab == tab);
  }
}

}  // namespace test_CliffTableau
}  // namespace tket
