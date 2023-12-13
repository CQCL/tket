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
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
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
static qubit_map_t inv_perm(const qubit_map_t& perm) {
  qubit_map_t inv;
  for (const std::pair<const Qubit, Qubit>& qp : perm)
    inv.insert({qp.second, qp.first});
  return inv;
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
            SpPauliStabiliser(Qubit(0), Pauli::Y)});
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
    // UnitaryRevTableau
    REQUIRE(
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Y), {}});
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
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Y), {}});
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
        tab.get_row(2) == ChoiMixTableau::row_tensor_t{
                              SpPauliStabiliser(Qubit(1), Pauli::Y), {}});
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
    Circuit res = cm_tableau_to_exact_circuit(tab).first;
    REQUIRE(res == circ);
    res = cm_tableau_to_unitary_extension_circuit(tab).first;
    REQUIRE(res == circ);
  }
  GIVEN("Just some Pauli gates for phase tests") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::Z, {3});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_exact_circuit(tab).first;
    REQUIRE(res == circ);
    res = cm_tableau_to_unitary_extension_circuit(tab).first;
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
      Circuit res = cm_tableau_to_exact_circuit(tab).first;
      Circuit res_uni = cm_tableau_to_unitary_extension_circuit(tab).first;
      REQUIRE(res == res_uni);
      ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
      REQUIRE(res_tab == tab);
      REQUIRE(test_unitary_comparison(circ, res, true));
    }
  }
  GIVEN("A unitary circuit") {
    Circuit circ = get_test_circ();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    std::pair<Circuit, qubit_map_t> res = cm_tableau_to_exact_circuit(tab);
    res.first.permute_boundary_output(inv_perm(res.second));
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(tab);
    res_uni.first.permute_boundary_output(inv_perm(res_uni.second));
    REQUIRE(res.first == res_uni.first);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    REQUIRE(res_tab == tab);
    REQUIRE(test_unitary_comparison(circ, res.first, true));
  }
  GIVEN("Check unitary equivalence by calculating matrix") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::ZZMax, {0, 1});
    circ.add_op<unsigned>(OpType::ECR, {2, 3});
    circ.add_op<unsigned>(OpType::ISWAPMax, {0, 3});
    circ.add_op<unsigned>(OpType::SX, {1});
    circ.add_op<unsigned>(OpType::SXdg, {2});
    circ.add_op<unsigned>(OpType::CY, {1, 3});
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    std::pair<Circuit, qubit_map_t> res = cm_tableau_to_exact_circuit(tab);
    res.first.permute_boundary_output(inv_perm(res.second));
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(tab);
    res_uni.first.permute_boundary_output(inv_perm(res_uni.second));
    REQUIRE(res.first == res_uni.first);
    REQUIRE(test_unitary_comparison(circ, res.first, true));
    THEN("Build the tableau manually for apply_gate coverage on inputs") {
      ChoiMixTableau rev_tab(4);
      rev_tab.apply_gate(
          OpType::CY, {Qubit(1), Qubit(3)},
          ChoiMixTableau::TableauSegment::Input);
      rev_tab.apply_gate(
          OpType::SXdg, {Qubit(2)}, ChoiMixTableau::TableauSegment::Input);
      rev_tab.apply_gate(
          OpType::SX, {Qubit(1)}, ChoiMixTableau::TableauSegment::Input);
      rev_tab.apply_gate(
          OpType::ISWAPMax, {Qubit(0), Qubit(3)},
          ChoiMixTableau::TableauSegment::Input);
      rev_tab.apply_gate(
          OpType::ECR, {Qubit(2), Qubit(3)},
          ChoiMixTableau::TableauSegment::Input);
      rev_tab.apply_gate(
          OpType::ZZMax, {Qubit(0), Qubit(1)},
          ChoiMixTableau::TableauSegment::Input);
      rev_tab.canonical_column_order();
      rev_tab.gaussian_form();
      REQUIRE(tab == rev_tab);
    }
  }
  GIVEN("A Clifford state") {
    Circuit circ = get_test_circ();
    circ.qubit_create_all();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_exact_circuit(tab).first;
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(res_tab == tab);
    Circuit res_uni =
        cm_tableau_to_unitary_extension_circuit(tab, circ.all_qubits()).first;
    REQUIRE(test_statevector_comparison(res, res_uni, true));
  }
  GIVEN("A partial Clifford state (tests mixed initialisations)") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    circ.add_op<unsigned>(OpType::Collapse, {1});
    circ.qubit_create_all();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_exact_circuit(tab).first;
    CHECK(res.created_qubits().size() == 3);
    CHECK(res.discarded_qubits().size() == 0);
    CHECK(res.count_gates(OpType::Collapse) == 1);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(res_tab == tab);
    Circuit res_uni =
        cm_tableau_to_unitary_extension_circuit(tab, circ.all_qubits()).first;
    Eigen::VectorXcd res_sv = tket_sim::get_statevector(res_uni);
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      Eigen::MatrixXcd outmat = rrow.second.to_sparse_matrix(3);
      CHECK((outmat * res_sv).isApprox(res_sv));
    }
  }
  GIVEN("A total diagonalisation circuit") {
    Circuit circ = get_test_circ();
    for (unsigned i = 0; i < circ.n_qubits(); ++i) {
      circ.add_op<unsigned>(OpType::Collapse, {i});
    }
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    Circuit res = cm_tableau_to_exact_circuit(tab).first;
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(res_tab == tab);
    // Test unitary synthesis by statevector of dagger
    Circuit as_state = get_test_circ().dagger();
    Circuit res_uni_dag =
        cm_tableau_to_unitary_extension_circuit(tab).first.dagger();
    REQUIRE(test_statevector_comparison(as_state, res_uni_dag, true));
  }
  GIVEN("A partial diagonalisation circuit") {
    Circuit circ = get_test_circ();
    for (unsigned i = 0; i < circ.n_qubits(); ++i) {
      circ.add_op<unsigned>(OpType::Collapse, {i});
    }
    circ.qubit_discard(Qubit(0));
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    std::pair<Circuit, qubit_map_t> res = cm_tableau_to_exact_circuit(tab);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    qubit_map_t perm;
    for (const std::pair<const Qubit, Qubit>& p : res.second) {
      perm.insert({p.second, p.first});
    }
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
    Circuit res_uni_dag =
        cm_tableau_to_unitary_extension_circuit(tab).first.dagger();
    Eigen::VectorXcd as_state = tket_sim::get_statevector(res_uni_dag);
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      CmplxSpMat rmat = rrow.first.to_sparse_matrix(3);
      if (rrow.second.is_real_negative()) rmat *= -1.;
      Eigen::MatrixXcd rmatd = rmat;
      CHECK((rmat * as_state).isApprox(as_state));
    }
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
    std::pair<Circuit, qubit_map_t> res = cm_tableau_to_exact_circuit(tab);
    res.first.permute_boundary_output(inv_perm(res.second));
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    REQUIRE(res_tab == tab);
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(tab);
    res_uni.first.permute_boundary_output(inv_perm(res_uni.second));
    res_tab = circuit_to_cm_tableau(res_uni.first);
    res_tab.tab_.row_mult(0, 1);
    Eigen::MatrixXcd res_u = tket_sim::get_unitary(res_uni.first);
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      CmplxSpMat inmat = rrow.first.to_sparse_matrix(5);
      Eigen::MatrixXcd inmatd = inmat;
      CmplxSpMat outmat = rrow.second.to_sparse_matrix(5);
      Eigen::MatrixXcd outmatd = outmat;
      CHECK((outmatd * res_u * inmatd).isApprox(res_u));
    }
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
    std::pair<Circuit, qubit_map_t> res = cm_tableau_to_exact_circuit(tab);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    qubit_map_t perm;
    for (const std::pair<const Qubit, Qubit>& p : res.second) {
      perm.insert({p.second, p.first});
    }
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(
            tab, {Qubit(1), Qubit(2), Qubit(3)});
    Eigen::MatrixXcd res_u = tket_sim::get_unitary(res_uni.first);
    Eigen::MatrixXcd init_proj = Eigen::MatrixXcd::Zero(32, 32);
    init_proj.block(0, 0, 2, 2) = Eigen::MatrixXcd::Identity(2, 2);
    init_proj.block(16, 16, 2, 2) = Eigen::MatrixXcd::Identity(2, 2);
    Eigen::MatrixXcd res_iso = res_u * init_proj;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      CmplxSpMat inmat = rrow.first.to_sparse_matrix(5);
      Eigen::MatrixXcd inmatd = inmat;
      QubitPauliMap outstr;
      for (const std::pair<const Qubit, Pauli>& qp : rrow.second.string)
        outstr.insert({res_uni.second.at(qp.first), qp.second});
      CmplxSpMat outmat = SpPauliString(outstr).to_sparse_matrix(5);
      Eigen::MatrixXcd outmatd = outmat;
      if (rrow.second.is_real_negative()) outmatd *= -1.;
      CHECK((outmatd * res_iso * inmatd).isApprox(res_iso));
    }
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
    std::pair<Circuit, qubit_map_t> res = cm_tableau_to_exact_circuit(tab);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res.first);
    qubit_map_t perm;
    for (const std::pair<const Qubit, Qubit>& p : res.second) {
      perm.insert({p.second, p.first});
    }
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(
            tab, {Qubit(1), Qubit(2), Qubit(3)});
    Eigen::MatrixXcd res_u = tket_sim::get_unitary(res_uni.first);
    Eigen::MatrixXcd init_proj = Eigen::MatrixXcd::Zero(32, 32);
    init_proj.block(0, 0, 2, 2) = Eigen::MatrixXcd::Identity(2, 2);
    init_proj.block(16, 16, 2, 2) = Eigen::MatrixXcd::Identity(2, 2);
    Eigen::MatrixXcd res_iso = res_u * init_proj;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      CmplxSpMat inmat = rrow.first.to_sparse_matrix(5);
      Eigen::MatrixXcd inmatd = inmat;
      QubitPauliMap outstr;
      for (const std::pair<const Qubit, Pauli>& qp : rrow.second.string)
        outstr.insert({res_uni.second.at(qp.first), qp.second});
      CmplxSpMat outmat = SpPauliString(outstr).to_sparse_matrix(5);
      Eigen::MatrixXcd outmatd = outmat;
      if (rrow.second.is_real_negative()) outmatd *= -1.;
      CHECK((outmatd * res_iso * inmatd).isApprox(res_iso));
    }
  }
  GIVEN("Synthesising a tableau requiring post-selection") {
    Circuit circ = get_test_circ();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    tab.post_select(Qubit(0), ChoiMixTableau::TableauSegment::Output);
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(tab, {}, {Qubit(0)});
    Eigen::MatrixXcd res_u = tket_sim::get_unitary(res_uni.first);
    // q[0] was removed from the tableau by postselection so need to infer
    // position in res_uni.second from the other qubits
    SpPauliString zzz({Pauli::Z, Pauli::Z, Pauli::Z});
    zzz.set(res_uni.second.at(Qubit(1)), Pauli::I);
    zzz.set(res_uni.second.at(Qubit(2)), Pauli::I);
    Eigen::MatrixXcd z0 = zzz.to_sparse_matrix(3);
    Eigen::MatrixXcd res_proj = 0.5 * (res_u + (z0 * res_u));
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      CmplxSpMat inmat = rrow.first.to_sparse_matrix(3);
      Eigen::MatrixXcd inmatd = inmat;
      QubitPauliMap outstr;
      for (const std::pair<const Qubit, Pauli>& qp : rrow.second.string)
        outstr.insert({res_uni.second.at(qp.first), qp.second});
      CmplxSpMat outmat = SpPauliString(outstr).to_sparse_matrix(3);
      Eigen::MatrixXcd outmatd = outmat;
      if (rrow.second.is_real_negative()) outmatd *= -1.;
      CHECK((outmatd * res_proj * inmatd).isApprox(res_proj));
    }
  }
  GIVEN("Synthesising a tableau with all post-selections") {
    Circuit circ = get_test_circ();
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    tab.post_select(Qubit(0), ChoiMixTableau::TableauSegment::Output);
    tab.post_select(Qubit(1), ChoiMixTableau::TableauSegment::Output);
    tab.post_select(Qubit(2), ChoiMixTableau::TableauSegment::Output);
    Circuit res = cm_tableau_to_unitary_extension_circuit(
                      tab, {}, {Qubit(0), Qubit(1), Qubit(2)})
                      .first.dagger();
    Eigen::VectorXcd res_sv = tket_sim::get_statevector(res);
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      Eigen::MatrixXcd inmat = rrow.first.to_sparse_matrix(3);
      if (rrow.second.is_real_negative()) inmat *= -1.;
      CHECK((inmat * res_sv).isApprox(res_sv));
    }
  }
  GIVEN("Initialisations, collapses, discards and post-selections") {
    Circuit circ(5);
    circ.qubit_create(Qubit(1));
    circ.qubit_create(Qubit(2));
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
    circ.qubit_discard(Qubit(0));
    ChoiMixTableau tab = circuit_to_cm_tableau(circ);
    tab.post_select(Qubit(3), ChoiMixTableau::TableauSegment::Output);
    tab.canonical_column_order();
    tab.gaussian_form();
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(tab, {Qubit(1)}, {Qubit(0)});
    // First rebuild tableau by initialising, post-selecting, etc.
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res_uni.first);
    qubit_map_t perm;
    for (const std::pair<const Qubit, Qubit>& p : res_uni.second)
      perm.insert({p.second, p.first});
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    // Post-select/initialise
    res_tab.post_select(Qubit(1), ChoiMixTableau::TableauSegment::Input);
    res_tab.post_select(Qubit(0), ChoiMixTableau::TableauSegment::Output);
    // Collapsing q[4] in X basis as per circ
    res_tab.apply_gate(
        OpType::H, {Qubit(4)}, ChoiMixTableau::TableauSegment::Output);
    res_tab.collapse_qubit(Qubit(4), ChoiMixTableau::TableauSegment::Output);
    res_tab.apply_gate(
        OpType::H, {Qubit(4)}, ChoiMixTableau::TableauSegment::Output);
    // Discarding q[0] also removes Z row for q[0], so recreate this by
    // XCollapse at input
    res_tab.apply_gate(
        OpType::H, {Qubit(0)}, ChoiMixTableau::TableauSegment::Input);
    res_tab.collapse_qubit(Qubit(0), ChoiMixTableau::TableauSegment::Input);
    res_tab.apply_gate(
        OpType::H, {Qubit(0)}, ChoiMixTableau::TableauSegment::Input);
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    REQUIRE(res_tab == tab);

    Eigen::MatrixXcd res_u = tket_sim::get_unitary(res_uni.first);
    qubit_vector_t res_qbs = res_uni.first.all_qubits();
    // q[1] has no input terms, so initialise it
    SpPauliString z1({Qubit(1), Pauli::Z});
    Eigen::MatrixXcd z1u = z1.to_sparse_matrix(res_qbs);
    res_u = 0.5 * (res_u + (res_u * z1u));
    // q[0] has no output terms, so postselect it
    SpPauliString z0({res_uni.second.at(Qubit(0)), Pauli::Z});
    Eigen::MatrixXcd z0u = z0.to_sparse_matrix(res_qbs);
    res_u = 0.5 * (res_u + (z0u * res_u));

    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      Eigen::MatrixXcd inmat = rrow.first.to_sparse_matrix(res_qbs);
      QubitPauliMap outstr;
      for (const std::pair<const Qubit, Pauli>& qp : rrow.second.string)
        outstr.insert({res_uni.second.at(qp.first), qp.second});
      Eigen::MatrixXcd outmat = SpPauliString(outstr).to_sparse_matrix(res_qbs);
      if (rrow.second.is_real_negative()) outmat *= -1.;
      CHECK((outmat * res_u * inmat).isApprox(res_u));
    }
  }
  GIVEN(
      "A custom tableau with overlapping initialised and post-selected "
      "qubits") {
    std::list<ChoiMixTableau::row_tensor_t> rows{
        {SpPauliStabiliser({Pauli::Z, Pauli::X, Pauli::I}), {}},
        {SpPauliStabiliser({Pauli::X, Pauli::Y, Pauli::Z}), {}},
        {{}, SpPauliStabiliser({Pauli::X, Pauli::X, Pauli::I})},
        {{}, SpPauliStabiliser({Pauli::I, Pauli::X, Pauli::X})},
        {SpPauliStabiliser({Pauli::I, Pauli::I, Pauli::Z}),
         SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::Z})},
        {SpPauliStabiliser({Pauli::Z, Pauli::I, Pauli::X}),
         SpPauliStabiliser({Pauli::I, Pauli::I, Pauli::X})},
    };
    ChoiMixTableau tab(rows);
    REQUIRE_THROWS(cm_tableau_to_unitary_extension_circuit(tab));
    std::pair<Circuit, qubit_map_t> res_uni =
        cm_tableau_to_unitary_extension_circuit(
            tab, {Qubit(3), Qubit(4)}, {Qubit(3), Qubit(4)});

    ChoiMixTableau res_tab = circuit_to_cm_tableau(res_uni.first);
    qubit_map_t perm;
    for (const std::pair<const Qubit, Qubit>& p : res_uni.second)
      perm.insert({p.second, p.first});
    res_tab.rename_qubits(perm, ChoiMixTableau::TableauSegment::Output);
    res_tab.post_select(Qubit(3), ChoiMixTableau::TableauSegment::Input);
    res_tab.post_select(Qubit(4), ChoiMixTableau::TableauSegment::Input);
    res_tab.post_select(Qubit(3), ChoiMixTableau::TableauSegment::Output);
    res_tab.post_select(Qubit(4), ChoiMixTableau::TableauSegment::Output);
    res_tab.canonical_column_order();
    res_tab.gaussian_form();
    tab.canonical_column_order();
    tab.gaussian_form();
    REQUIRE(res_tab == tab);

    Eigen::MatrixXcd res_u = tket_sim::get_unitary(res_uni.first);
    qubit_vector_t res_qbs = res_uni.first.all_qubits();
    // initialise q[3] and q[4]
    SpPauliString z3i{Qubit(3), Pauli::Z};
    Eigen::MatrixXcd z3iu = z3i.to_sparse_matrix(res_qbs);
    res_u = 0.5 * (res_u + (res_u * z3iu));
    SpPauliString z4i{Qubit(4), Pauli::Z};
    Eigen::MatrixXcd z4iu = z4i.to_sparse_matrix(res_qbs);
    res_u = 0.5 * (res_u + (res_u * z4iu));
    // post-select q[3] and q[4]
    SpPauliString z3o{res_uni.second.at(Qubit(3)), Pauli::Z};
    Eigen::MatrixXcd z3ou = z3o.to_sparse_matrix(res_qbs);
    res_u = 0.5 * (res_u + (z3ou * res_u));
    SpPauliString z4o{res_uni.second.at(Qubit(4)), Pauli::Z};
    Eigen::MatrixXcd z4ou = z4o.to_sparse_matrix(res_qbs);
    res_u = 0.5 * (res_u + (z4ou * res_u));

    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      ChoiMixTableau::row_tensor_t rrow = tab.get_row(r);
      Eigen::MatrixXcd inmat = rrow.first.to_sparse_matrix(res_qbs);
      QubitPauliMap outstr;
      for (const std::pair<const Qubit, Pauli>& qp : rrow.second.string)
        outstr.insert({res_uni.second.at(qp.first), qp.second});
      Eigen::MatrixXcd outmat = SpPauliString(outstr).to_sparse_matrix(res_qbs);
      if (rrow.second.is_real_negative()) outmat *= -1.;
      CHECK((outmat * res_u * inmat).isApprox(res_u));
    }
  }
}

SCENARIO("Conversions to and from UnitaryTableau and UnitaryRevTableau") {
  GIVEN("A round trip UnitaryTableau -> ChoiMixTableau -> UnitaryTableau") {
    Circuit circ = get_test_circ();
    UnitaryTableau utab = circuit_to_unitary_tableau(circ);
    ChoiMixTableau cmtab = unitary_tableau_to_cm_tableau(utab);
    UnitaryTableau utab2 = cm_tableau_to_unitary_tableau(cmtab);
    REQUIRE(utab == utab2);
  }
  GIVEN(
      "A round trip UnitaryRevTableau -> ChoiMixTableau -> UnitaryRevTableau") {
    Circuit circ = get_test_circ();
    UnitaryRevTableau utab = circuit_to_unitary_rev_tableau(circ);
    ChoiMixTableau cmtab = unitary_rev_tableau_to_cm_tableau(utab);
    UnitaryRevTableau utab2 = cm_tableau_to_unitary_rev_tableau(cmtab);
    REQUIRE(utab == utab2);
  }
  GIVEN("A non-unitary ChoiMixTableau") {
    Circuit circ = get_test_circ();
    circ.qubit_discard(Qubit(1));
    circ.qubit_create(Qubit(2));
    ChoiMixTableau cmtab = circuit_to_cm_tableau(circ);
    REQUIRE_THROWS(cm_tableau_to_unitary_tableau(cmtab));
    REQUIRE_THROWS(cm_tableau_to_unitary_rev_tableau(cmtab));
  }
}

}  // namespace test_ChoiMixTableau
}  // namespace tket
