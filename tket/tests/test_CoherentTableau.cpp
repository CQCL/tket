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

#include "Converters/Converters.hpp"

namespace tket {
namespace test_CoherentTableau {

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
static CoherentTableau get_tableau_with_gates_applied_at_front() {
  CoherentTableau tab(3);
  tab.apply_gate(
      OpType::CX, {Qubit(1), Qubit(0)}, CoherentTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::Vdg, {Qubit(1)}, CoherentTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::CX, {Qubit(1), Qubit(2)}, CoherentTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::CX, {Qubit(0), Qubit(1)}, CoherentTableau::TableauSegment::Input);
  tab.apply_gate(OpType::S, {Qubit(1)}, CoherentTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::CX, {Qubit(0), Qubit(1)}, CoherentTableau::TableauSegment::Input);
  return tab;
}

SCENARIO("Correct creation of CoherentTableau") {
  GIVEN(
      "A circuit with an identity, a discarded input, and an initialised "
      "output") {
    Circuit circ(3);
    circ.qubit_discard(Qubit(1));
    circ.qubit_create(Qubit(2));
    CoherentTableau tab = circuit_to_coherent_tableau(circ);
    REQUIRE(tab.get_n_rows() == 3);
    REQUIRE(tab.get_n_boundaries() == 4);
    REQUIRE(tab.get_n_inputs() == 2);
    REQUIRE(tab.get_n_outputs() == 2);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::X),
                              QubitPauliTensor(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(2), Pauli::Z)});
  }
  GIVEN("Applying S gates") {
    CoherentTableau tab(3);
    tab.post_select(Qubit(1), CoherentTableau::TableauSegment::Output);
    tab.post_select(Qubit(2), CoherentTableau::TableauSegment::Input);
    // Check S on initialised/post-selected qubits does nothing
    CoherentTableau orig = tab;
    tab.apply_S(Qubit(1), CoherentTableau::TableauSegment::Input);
    tab.apply_S(Qubit(2), CoherentTableau::TableauSegment::Output);
    tab.gaussian_form();
    REQUIRE(tab == orig);
    // Check S on identity
    tab.apply_S(Qubit(0));
    REQUIRE(
        tab.get_row(0) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::X),
                              QubitPauliTensor(Qubit(0), Pauli::Y, -1.)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::Z), {}});
    REQUIRE(
        tab.get_row(3) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(2), Pauli::Z)});
    // Applying an S at the input end adds up to a net Z
    tab.apply_S(Qubit(0), CoherentTableau::TableauSegment::Input);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::X),
                              QubitPauliTensor(Qubit(0), Pauli::X, -1.)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::Z), {}});
    REQUIRE(
        tab.get_row(3) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(2), Pauli::Z)});
  }
  GIVEN("Applying V gates") {
    CoherentTableau tab(3);
    tab.post_select(Qubit(1), CoherentTableau::TableauSegment::Output);
    tab.post_select(Qubit(2), CoherentTableau::TableauSegment::Input);
    // V on initialised/post-selected qubits has non-trivial effect
    tab.apply_V(Qubit(1), CoherentTableau::TableauSegment::Input);
    tab.apply_V(Qubit(2), CoherentTableau::TableauSegment::Output);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::X),
                              QubitPauliTensor(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::Y), {}});
    REQUIRE(
        tab.get_row(3) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(2), Pauli::Y)});
    // Check V on identity
    tab.apply_V(Qubit(0), CoherentTableau::TableauSegment::Output);
    REQUIRE(
        tab.get_row(0) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::X),
                              QubitPauliTensor(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Y)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::Y), {}});
    REQUIRE(
        tab.get_row(3) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(2), Pauli::Y)});
    // Applying a V at the input end adds up to a net X
    tab.apply_V(Qubit(0), CoherentTableau::TableauSegment::Input);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::X),
                              QubitPauliTensor(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z, -1.)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::Y), {}});
    REQUIRE(
        tab.get_row(3) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(2), Pauli::Y)});
  }
  GIVEN("Applying CX gates") {
    CoherentTableau tab(4);
    tab.post_select(Qubit(2), CoherentTableau::TableauSegment::Output);
    tab.post_select(Qubit(3), CoherentTableau::TableauSegment::Input);
    // CX with control on initialised/post-selected qubits does nothing
    CoherentTableau orig = tab;
    tab.apply_CX(Qubit(2), Qubit(0), CoherentTableau::TableauSegment::Input);
    tab.apply_CX(Qubit(3), Qubit(1), CoherentTableau::TableauSegment::Output);
    tab.gaussian_form();
    REQUIRE(tab == orig);
    // Check CX on identity
    tab.apply_CX(Qubit(0), Qubit(1));
    REQUIRE(
        tab.get_row(0) ==
        CoherentTableau::row_tensor_t{
            QubitPauliTensor(Qubit(0), Pauli::X),
            QubitPauliTensor({{Qubit(0), Pauli::X}, {Qubit(1), Pauli::X}})});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::X),
                              QubitPauliTensor(Qubit(1), Pauli::X)});
    REQUIRE(
        tab.get_row(3) ==
        CoherentTableau::row_tensor_t{
            QubitPauliTensor(Qubit(1), Pauli::Z),
            QubitPauliTensor({{Qubit(0), Pauli::Z}, {Qubit(1), Pauli::Z}})});
    REQUIRE(
        tab.get_row(4) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(2), Pauli::Z), {}});
    REQUIRE(
        tab.get_row(5) == CoherentTableau::row_tensor_t{
                              {}, QubitPauliTensor(Qubit(3), Pauli::Z)});
    // CX on input cancels back to original
    tab.apply_CX(Qubit(0), Qubit(1), CoherentTableau::TableauSegment::Input);
    tab.gaussian_form();
    REQUIRE(tab == orig);
    // CX with target on initialised/post-selected qubit still entangles
    tab.apply_CX(Qubit(0), Qubit(2), CoherentTableau::TableauSegment::Input);
    tab.apply_CX(Qubit(1), Qubit(3));
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) ==
        CoherentTableau::row_tensor_t{
            QubitPauliTensor({{Qubit(0), Pauli::X}, {Qubit(2), Pauli::X}}),
            QubitPauliTensor(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(0), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) ==
        CoherentTableau::row_tensor_t{
            QubitPauliTensor(Qubit(1), Pauli::X),
            QubitPauliTensor({{Qubit(1), Pauli::X}, {Qubit(3), Pauli::X}})});
    REQUIRE(
        tab.get_row(3) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(1), Pauli::Z),
                              QubitPauliTensor(Qubit(1), Pauli::Z)});
    REQUIRE(
        tab.get_row(4) == CoherentTableau::row_tensor_t{
                              QubitPauliTensor(Qubit(2), Pauli::Z),
                              QubitPauliTensor(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(5) ==
        CoherentTableau::row_tensor_t{
            {},
            QubitPauliTensor({{Qubit(1), Pauli::Z}, {Qubit(3), Pauli::Z}})});
  }
  GIVEN("A full circuit") {
    Circuit circ = get_test_circ();
    CoherentTableau tab = circuit_to_coherent_tableau(circ);
    CoherentTableau rev_tab = get_tableau_with_gates_applied_at_front();
    tab.gaussian_form();
    rev_tab.gaussian_form();
    REQUIRE(tab == rev_tab);
  }
  GIVEN("A PI/2 rotation at front") {
    Circuit circ = get_test_circ();
    CoherentTableau tab = circuit_to_coherent_tableau(circ);
    QubitPauliTensor pauli{
        {{Qubit(0), Pauli::X}, {Qubit(1), Pauli::Y}, {Qubit(2), Pauli::Z}}};
    tab.apply_pauli(pauli, 3);
    tab.gaussian_form();

    add_ops_list_two_to_circuit(circ, OpType::Sdg);
    CoherentTableau correct_tab = circuit_to_coherent_tableau(circ);
    correct_tab.gaussian_form();
    REQUIRE(tab == correct_tab);
  }
  GIVEN("A PI/2 rotation at front") {
    CoherentTableau tab = get_tableau_with_gates_applied_at_front();
    QubitPauliTensor pauli =
        QubitPauliTensor(Qubit(q_default_reg(), 0), Pauli::X) *
        QubitPauliTensor(Qubit(q_default_reg(), 1), Pauli::Y) *
        QubitPauliTensor(Qubit(q_default_reg(), 2), Pauli::Z);
    tab.apply_pauli(pauli, 1, CoherentTableau::TableauSegment::Input);
    tab.gaussian_form();

    Circuit circ(3);
    add_ops_list_two_to_circuit(circ);
    add_ops_list_one_to_circuit(circ);
    CoherentTableau correct_tab = circuit_to_coherent_tableau(circ);
    correct_tab.gaussian_form();
    REQUIRE(tab == correct_tab);
  }
  GIVEN("Combining two non-unitary circuits via tableau compose") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    circ.qubit_discard(Qubit(2));
    CoherentTableau first = circuit_to_coherent_tableau(circ);

    Circuit circ1(3);
    add_ops_list_two_to_circuit(circ1);
    circ1.qubit_create(Qubit(2));
    CoherentTableau second = circuit_to_coherent_tableau(circ1);
    CoherentTableau correct = circuit_to_coherent_tableau(circ >> circ1);
    CoherentTableau result = CoherentTableau::compose(first, second);
    result.canonical_column_order();
    result.gaussian_form();
    correct.canonical_column_order();
    correct.gaussian_form();
    REQUIRE(result == correct);
  }
}

SCENARIO("Error handling in CoherentTableau generation") {
  GIVEN("Add a non-clifford gate at end") {
    CoherentTableau tab(2);
    REQUIRE_THROWS_AS(tab.apply_gate(OpType::T, {Qubit(0)}), BadOpType);
  }
  GIVEN("Add a non-clifford gate at front") {
    CoherentTableau tab(2);
    REQUIRE_THROWS_AS(
        tab.apply_gate(
            OpType::Tdg, {Qubit(0)}, CoherentTableau::TableauSegment::Input),
        BadOpType);
  }
  GIVEN("Tableau from a non-Clifford circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CH, {1, 0});
    REQUIRE_THROWS_AS(circuit_to_coherent_tableau(circ), BadOpType);
  }
}

SCENARIO("Synthesis of circuits from CoherentTableaus") {
  GIVEN("A unitary circuit") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    CliffTableau tab = circuit_to_tableau(circ);
    Circuit res = tableau_to_circuit(tab);
    CliffTableau res_tab = circuit_to_tableau(res);
    REQUIRE(res_tab == tab);
  }
}

}  // namespace test_CoherentTableau
}  // namespace tket
