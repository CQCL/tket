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

#include "testutil.hpp"
#include "tket/Converters/Converters.hpp"

namespace tket {
namespace test_ChoiMixTableau {

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
static ChoiMixTableau get_tableau_with_gates_applied_at_front() {
  ChoiMixTableau tab(3);
  tab.apply_gate(
      OpType::CX, {Qubit(1), Qubit(0)}, ChoiMixTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::Vdg, {Qubit(1)}, ChoiMixTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::CX, {Qubit(1), Qubit(2)}, ChoiMixTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::CX, {Qubit(0), Qubit(1)}, ChoiMixTableau::TableauSegment::Input);
  tab.apply_gate(OpType::S, {Qubit(1)}, ChoiMixTableau::TableauSegment::Input);
  tab.apply_gate(
      OpType::CX, {Qubit(0), Qubit(1)}, ChoiMixTableau::TableauSegment::Input);
  return tab;
}

SCENARIO("Correct creation of ChoiMixTableau") {
  GIVEN(
      "A circuit with an identity, a discarded input, and an initialised "
      "output") {
    Circuit circ(3);
    circ.qubit_discard(Qubit(1));
    circ.qubit_create(Qubit(2));
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    REQUIRE(tab.get_n_rows() == 3);
    REQUIRE(tab.get_n_boundaries() == 4);
    REQUIRE(tab.get_n_inputs() == 2);
    REQUIRE(tab.get_n_outputs() == 2);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::X),
                              SpPauliStabiliser(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(2), Pauli::Z)});
    REQUIRE(
        tab.get_row_product({0, 1}) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(0), Pauli::Y),
            SpPauliStabiliser(Qubit(0), Pauli::Y, 2)});
    THEN("Serialize and deserialize") {
      nlohmann::json j_tab = tab;
      ChoiMixTableau tab2{{}};
      j_tab.get_to(tab2);
      REQUIRE(tab == tab2);
    }
  }
  GIVEN("Applying S gates") {
    ChoiMixTableau tab(3);
    tab.post_select(Qubit(1), ChoiMixTableau::TableauSegment::Output);
    tab.post_select(Qubit(2), ChoiMixTableau::TableauSegment::Input);
    // Check S on initialised/post-selected qubits does nothing
    ChoiMixTableau orig = tab;
    tab.apply_S(Qubit(1), ChoiMixTableau::TableauSegment::Input);
    tab.apply_S(Qubit(2), ChoiMixTableau::TableauSegment::Output);
    tab.gaussian_form();
    REQUIRE(tab == orig);
    // Check S on identity
    tab.apply_S(Qubit(0));
    // e^{-i Z pi/4} X id X = (-iZX) e^{-i Z pi/4} X = +Y S X
    REQUIRE(
        tab.get_row(0) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::X),
                              SpPauliStabiliser(Qubit(0), Pauli::Y)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Z), {}});
    REQUIRE(
        tab.get_row(3) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(2), Pauli::Z)});
    // Applying an S at the input end adds up to a net Z
    tab.apply_S(Qubit(0), ChoiMixTableau::TableauSegment::Input);
    tab.canonical_column_order();
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::X),
                              SpPauliStabiliser(Qubit(0), Pauli::X, 2)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Z), {}});
    REQUIRE(
        tab.get_row(3) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(2), Pauli::Z)});
    THEN("Compare to explicitly generated tableau") {
      std::list<ChoiMixTableau::row_tensor_t> rows;
      rows.push_back(
          {SpPauliStabiliser(Qubit(0), Pauli::X),
           SpPauliStabiliser(Qubit(0), Pauli::X, 2)});
      rows.push_back(
          {SpPauliStabiliser(Qubit(0), Pauli::Z),
           SpPauliStabiliser(Qubit(0), Pauli::Z)});
      rows.push_back({SpPauliStabiliser(Qubit(1), Pauli::Z), {}});
      rows.push_back({{}, SpPauliStabiliser(Qubit(2), Pauli::Z)});
      ChoiMixTableau tab2(rows);
      tab2.canonical_column_order();
      REQUIRE(tab == tab2);
    }
  }
  GIVEN("Applying V gates") {
    ChoiMixTableau tab(3);
    tab.post_select(Qubit(1), ChoiMixTableau::TableauSegment::Output);
    tab.post_select(Qubit(2), ChoiMixTableau::TableauSegment::Input);
    // V on initialised/post-selected qubits has non-trivial effect
    tab.apply_V(Qubit(1), ChoiMixTableau::TableauSegment::Input);
    tab.apply_V(Qubit(2), ChoiMixTableau::TableauSegment::Output);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::X),
                              SpPauliStabiliser(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    // Affecting the input segment should give the same effect as for
    // UnitaryRevTableau (since lhs is transposed, +Y is flipped to -Y, and
    // phase is returned on rhs)
    REQUIRE(
        tab.get_row(2) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(1), Pauli::Y), SpPauliStabiliser({}, 2)});
    // Affecting the output segment should give the same effect as for
    // UnitaryTableau
    REQUIRE(
        tab.get_row(3) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(2), Pauli::Y, 2)});
    // Check V on identity
    tab.apply_V(Qubit(0), ChoiMixTableau::TableauSegment::Output);
    REQUIRE(
        tab.get_row(0) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::X),
                              SpPauliStabiliser(Qubit(0), Pauli::X)});
    // e^{-i X pi/4} Z C Z = (-iXZ) e^{-i X pi/4} C Z = -Y V C Z
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Y, 2)});
    REQUIRE(
        tab.get_row(2) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(1), Pauli::Y), SpPauliStabiliser({}, 2)});
    REQUIRE(
        tab.get_row(3) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(2), Pauli::Y, 2)});
    // Applying a V at the input end adds up to a net X
    tab.apply_V(Qubit(0), ChoiMixTableau::TableauSegment::Input);
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::X),
                              SpPauliStabiliser(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z, 2)});
    REQUIRE(
        tab.get_row(2) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(1), Pauli::Y), SpPauliStabiliser({}, 2)});
    REQUIRE(
        tab.get_row(3) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(2), Pauli::Y, 2)});
  }
  GIVEN("Applying CX gates") {
    ChoiMixTableau tab(4);
    tab.post_select(Qubit(2), ChoiMixTableau::TableauSegment::Output);
    tab.post_select(Qubit(3), ChoiMixTableau::TableauSegment::Input);
    // CX with control on initialised/post-selected qubits does nothing
    ChoiMixTableau orig = tab;
    tab.apply_CX(Qubit(2), Qubit(0), ChoiMixTableau::TableauSegment::Input);
    tab.apply_CX(Qubit(3), Qubit(1), ChoiMixTableau::TableauSegment::Output);
    tab.gaussian_form();
    REQUIRE(tab == orig);
    // Check CX on identity
    tab.apply_CX(Qubit(0), Qubit(1));
    REQUIRE(
        tab.get_row(0) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(0), Pauli::X),
            SpPauliStabiliser({{Qubit(0), Pauli::X}, {Qubit(1), Pauli::X}})});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::X),
                              SpPauliStabiliser(Qubit(1), Pauli::X)});
    REQUIRE(
        tab.get_row(3) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(1), Pauli::Z),
            SpPauliStabiliser({{Qubit(0), Pauli::Z}, {Qubit(1), Pauli::Z}})});
    REQUIRE(
        tab.get_row(4) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(2), Pauli::Z), {}});
    REQUIRE(
        tab.get_row(5) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(3), Pauli::Z)});
    // CX on input cancels back to original
    tab.apply_CX(Qubit(0), Qubit(1), ChoiMixTableau::TableauSegment::Input);
    tab.gaussian_form();
    REQUIRE(tab == orig);
    // CX with target on initialised/post-selected qubit still entangles
    tab.apply_CX(Qubit(0), Qubit(2), ChoiMixTableau::TableauSegment::Input);
    tab.apply_CX(Qubit(1), Qubit(3));
    tab.gaussian_form();
    REQUIRE(
        tab.get_row(0) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser({{Qubit(0), Pauli::X}, {Qubit(2), Pauli::X}}),
            SpPauliStabiliser(Qubit(0), Pauli::X)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(0), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(2) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(1), Pauli::X),
            SpPauliStabiliser({{Qubit(1), Pauli::X}, {Qubit(3), Pauli::X}})});
    REQUIRE(
        tab.get_row(3) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Z),
                              SpPauliStabiliser(Qubit(1), Pauli::Z)});
    REQUIRE(
        tab.get_row(4) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(2), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z)});
    REQUIRE(
        tab.get_row(5) ==
        ChoiMixTableau::row_tensor_t{
            {},
            SpPauliStabiliser({{Qubit(1), Pauli::Z}, {Qubit(3), Pauli::Z}})});
  }
  GIVEN("A full circuit") {
    Circuit circ = get_test_circ();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    ChoiMixTableau rev_tab = get_tableau_with_gates_applied_at_front();
    tab.gaussian_form();
    rev_tab.gaussian_form();
    REQUIRE(tab == rev_tab);
  }
  GIVEN("A PI/2 rotation at end") {
    Circuit circ = get_test_circ();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    SpPauliStabiliser pauli{
        {{Qubit(0), Pauli::X}, {Qubit(1), Pauli::Y}, {Qubit(2), Pauli::Z}}};
    tab.apply_pauli(pauli, 3);
    tab.gaussian_form();

    add_ops_list_two_to_circuit(circ, OpType::Sdg);
    ChoiMixTableau correct_tab = circuit_to_cm_tableau(circ);
    correct_tab.gaussian_form();
    REQUIRE(tab == correct_tab);
  }
  GIVEN("A PI/2 rotation at front") {
    ChoiMixTableau tab = get_tableau_with_gates_applied_at_front();
    SpPauliStabiliser pauli({Pauli::X, Pauli::Y, Pauli::Z});
    tab.apply_pauli(pauli, 1, ChoiMixTableau::TableauSegment::Input);
    tab.gaussian_form();

    Circuit circ(3);
    add_ops_list_two_to_circuit(circ);
    add_ops_list_one_to_circuit(circ);
    ChoiMixTableau correct_tab = circuit_to_cm_tableau(circ);
    correct_tab.gaussian_form();
    REQUIRE(tab == correct_tab);
  }
  GIVEN("Combining two non-unitary circuits via tableau compose") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    circ.qubit_discard(Qubit(2));
    ChoiMixTableau first = circuit_to_cm_tableau(circ);

    Circuit circ1(3);
    add_ops_list_two_to_circuit(circ1);
    circ1.qubit_create(Qubit(2));
    ChoiMixTableau second = circuit_to_cm_tableau(circ1);
    ChoiMixTableau correct = circuit_to_cm_tableau(circ >> circ1);
    ChoiMixTableau result = ChoiMixTableau::compose(first, second);
    result.canonical_column_order();
    result.gaussian_form();
    correct.canonical_column_order();
    correct.gaussian_form();
    REQUIRE(result == correct);
  }
  GIVEN("Testing more gates") {
    ChoiMixTableau tab(3);
    tab.apply_gate(
        OpType::Y, {Qubit(0)}, ChoiMixTableau::TableauSegment::Input);
    tab.apply_gate(
        OpType::noop, {Qubit(0)}, ChoiMixTableau::TableauSegment::Input);
    tab.apply_gate(
        OpType::BRIDGE, {Qubit(0), Qubit(1), Qubit(2)},
        ChoiMixTableau::TableauSegment::Input);
    tab.apply_gate(
        OpType::SWAP, {Qubit(0), Qubit(1)},
        ChoiMixTableau::TableauSegment::Input);
    tab.apply_gate(
        OpType::Reset, {Qubit(0)}, ChoiMixTableau::TableauSegment::Input);

    tab.canonical_column_order();
    tab.gaussian_form();
    REQUIRE(tab.get_n_rows() == 5);
    REQUIRE(
        tab.get_row(0) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(1), Pauli::X),
            SpPauliStabiliser(
                {{Qubit(0), Pauli::X}, {Qubit(2), Pauli::X}}, 2)});
    REQUIRE(
        tab.get_row(1) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Z),
                              SpPauliStabiliser(Qubit(0), Pauli::Z, 2)});
    REQUIRE(
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(2), Pauli::X),
                              SpPauliStabiliser(Qubit(2), Pauli::X)});
    REQUIRE(
        tab.get_row(3) ==
        ChoiMixTableau::row_tensor_t{
            SpPauliStabiliser(Qubit(2), Pauli::Z),
            SpPauliStabiliser(
                {{Qubit(0), Pauli::Z}, {Qubit(2), Pauli::Z}}, 2)});
    REQUIRE(
        tab.get_row(4) == ChoiMixTableau::row_tensor_t{
                              {}, SpPauliStabiliser(Qubit(1), Pauli::Z)});
  }
  GIVEN("Combining post-selections and discarding") {
    ChoiMixTableau tab(5);
    // Post-selecting an initialised qubit succeeds deterministically
    tab.post_select(Qubit(1), ChoiMixTableau::TableauSegment::Input);
    tab.post_select(Qubit(1), ChoiMixTableau::TableauSegment::Output);
    // Post-selecting a mixed qubit succeeds probabilistically
    tab.discard_qubit(Qubit(2), ChoiMixTableau::TableauSegment::Input);
    tab.post_select(Qubit(2), ChoiMixTableau::TableauSegment::Output);
    // Discarding an initialised qubit
    tab.discard_qubit(Qubit(3), ChoiMixTableau::TableauSegment::Input);
    tab.post_select(Qubit(3), ChoiMixTableau::TableauSegment::Output);
    // Discarding a mixed qubit
    tab.discard_qubit(Qubit(4), ChoiMixTableau::TableauSegment::Input);
    tab.discard_qubit(Qubit(4), ChoiMixTableau::TableauSegment::Output);
    REQUIRE(tab == ChoiMixTableau(1));
    // Test that impossible post-selection fails
    tab.post_select(Qubit(0), ChoiMixTableau::TableauSegment::Input);
    tab.apply_gate(OpType::X, {Qubit(0)});
    REQUIRE_THROWS(
        tab.post_select(Qubit(0), ChoiMixTableau::TableauSegment::Output));
  }
}

SCENARIO("Error handling in ChoiMixTableau generation") {
  GIVEN("Exceptions in ChoiMixTableau constructors") {
    MatrixXb xmat = MatrixXb::Zero(3, 3);
    VectorXb ph = VectorXb::Zero(3);
    // Different size components
    REQUIRE_THROWS_AS(
        ChoiMixTableau(xmat, MatrixXb::Zero(2, 3), ph), std::invalid_argument);
    // Rows not independent
    MatrixXb zmat(3, 3);
    zmat << true, true, false, true, false, true, false, true, true;
    REQUIRE_THROWS_AS(ChoiMixTableau(xmat, zmat, ph), std::invalid_argument);
    // Rows don't commute
    zmat(2, 2) = false;
    xmat(0, 0) = true;
    REQUIRE_THROWS_AS(ChoiMixTableau(xmat, zmat, ph), std::invalid_argument);
  }
  GIVEN("Add a non-clifford gate at end") {
    ChoiMixTableau tab(2);
    REQUIRE_THROWS_AS(tab.apply_gate(OpType::T, {Qubit(0)}), BadOpType);
  }
  GIVEN("Add a non-clifford gate at front") {
    ChoiMixTableau tab(2);
    REQUIRE_THROWS_AS(
        tab.apply_gate(
            OpType::Tdg, {Qubit(0)}, ChoiMixTableau::TableauSegment::Input),
        BadOpType);
  }
  GIVEN("Tableau from a non-Clifford circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CH, {1, 0});
    REQUIRE_THROWS_AS(circuit_to_cm_tableau(circ), BadOpType);
  }
}

SCENARIO("Synthesis of circuits from ChoiMixTableaus") {
  GIVEN("An identity circuit") {
    Circuit circ(3);
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    REQUIRE(res == circ);
  }
  GIVEN("Just some Pauli gates for phase tests") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::Z, {3});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    REQUIRE(res == circ);
  }
  GIVEN("Iterate through single-qubit Cliffords with all entanglements") {
    for (unsigned i = 0; i < 27; ++i) {
      Circuit circ(7);
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      circ.add_op<unsigned>(OpType::CY, {0, 2});
      circ.add_op<unsigned>(OpType::CZ, {0, 3});
      circ.add_op<unsigned>(OpType::H, {0});
      circ.add_op<unsigned>(OpType::CX, {0, 4});
      circ.add_op<unsigned>(OpType::CY, {0, 5});
      circ.add_op<unsigned>(OpType::CZ, {0, 6});
      circ.add_op<unsigned>(OpType::H, {0});
      if (i % 3 == 1) circ.add_op<unsigned>(OpType::S, {0});
      if (i % 3 == 2) circ.add_op<unsigned>(OpType::Sdg, {0});
      if ((i / 3) % 3 == 1) circ.add_op<unsigned>(OpType::V, {0});
      if ((i / 3) % 3 == 2) circ.add_op<unsigned>(OpType::Vdg, {0});
      if ((i / 9) % 3 == 1) circ.add_op<unsigned>(OpType::S, {0});
      if ((i / 9) % 3 == 2) circ.add_op<unsigned>(OpType::Sdg, {0});
      ChoiMixTableau tab = circuit_to_cm_tableau(circ);
      Circuit res = cm_tableau_to_circuit(tab).first;
      ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
      REQUIRE(res_tab == tab);
      REQUIRE(test_unitary_comparison(circ, res, true));
    }
  }
  GIVEN("A unitary circuit") {
    Circuit circ = get_test_circ();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(res_tab == tab);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("Check unitary equivalence by calculating matrix") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::ZZMax, {0, 1});
    circ.add_op<unsigned>(OpType::ECR, {2, 3});
    circ.add_op<unsigned>(OpType::ISWAPMax, {0, 3});
    circ.add_op<unsigned>(OpType::SX, {1});
    circ.add_op<unsigned>(OpType::SXdg, {2});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A Clifford state") {
    Circuit circ = get_test_circ();
    circ.qubit_create_all();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    tab.canonical_column_order();
    tab.gaussian_form();
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
  }
  GIVEN("A total diagonalisation circuit") {
    Circuit circ = get_test_circ();
    for (unsigned i = 0; i < circ.n_qubits(); ++i) {
      circ.add_op<unsigned>(OpType::Collapse, {i});
    }
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(res_tab == tab);
  }
  GIVEN("A partial diagonalisation circuit") {
    Circuit circ = get_test_circ();
    for (unsigned i = 0; i < circ.n_qubits(); ++i) {
      circ.add_op<unsigned>(OpType::Collapse, {i});
    }
    circ.qubit_discard(Qubit(0));
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    std::pair<Circuit, unit_map_t> res = cm_tableau_to_circuit(tab);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    qubit_map_t perm;
    for (const std::pair<const UnitID, UnitID>& p : res.second) {
      perm.insert({Qubit(p.second), Qubit(p.first)});
    }
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    tab.canonical_column_order();
    tab.gaussian_form();
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
  }
  GIVEN("Another circuit for extra test coverage in row reductions") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::Collapse, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::CZ, {0, 3});
    circ.add_op<unsigned>(OpType::Collapse, {1});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CY, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 3});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.qubit_discard(Qubit(2));
    circ.qubit_discard(Qubit(3));
    circ.add_op<unsigned>(OpType::Collapse, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_circuit(tab).first;
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    tab.canonical_column_order();
    tab.gaussian_form();
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
  }
  GIVEN("An isometry") {
    Circuit circ(5);
    circ.qubit_create(Qubit(1));
    circ.qubit_create(Qubit(2));
    circ.qubit_create(Qubit(3));
    circ.add_op<unsigned>(OpType::Collapse, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 1});
    circ.add_op<unsigned>(OpType::CX, {4, 2});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    std::pair<Circuit, unit_map_t> res = cm_tableau_to_circuit(tab);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    qubit_map_t perm;
    for (const std::pair<const UnitID, UnitID>& p : res.second) {
      perm.insert({Qubit(p.second), Qubit(p.first)});
    }
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    tab.canonical_column_order();
    tab.gaussian_form();
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
  }
  GIVEN("Extra coverage for isometries") {
    Circuit circ(5);
    circ.qubit_create(Qubit(1));
    circ.qubit_create(Qubit(2));
    circ.qubit_create(Qubit(3));
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::Collapse, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 1});
    circ.add_op<unsigned>(OpType::CX, {4, 2});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    std::pair<Circuit, unit_map_t> res = cm_tableau_to_circuit(tab);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    qubit_map_t perm;
    for (const std::pair<const UnitID, UnitID>& p : res.second) {
      perm.insert({Qubit(p.second), Qubit(p.first)});
    }
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    tab.canonical_column_order();
    tab.gaussian_form();
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
  }
}

}  // namespace test_ChoiMixTableau
}  // namespace tket
