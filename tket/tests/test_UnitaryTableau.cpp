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
#include <sstream>

#include "Clifford/UnitaryTableau.hpp"
#include "Converters/Converters.hpp"
#include "Converters/UnitaryTableauBox.hpp"
#include "testutil.hpp"

namespace tket::test_UnitaryTableau {

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
static UnitaryTableau get_tableau_with_gates_applied_at_front() {
  UnitaryTableau tab(3);
  tab.apply_gate_at_front(OpType::CX, {Qubit(1), Qubit(0)});
  tab.apply_gate_at_front(OpType::Vdg, {Qubit(1)});
  tab.apply_gate_at_front(OpType::CX, {Qubit(1), Qubit(2)});
  tab.apply_gate_at_front(OpType::CX, {Qubit(0), Qubit(1)});
  tab.apply_gate_at_front(OpType::S, {Qubit(1)});
  tab.apply_gate_at_front(OpType::CX, {Qubit(0), Qubit(1)});
  return tab;
}

SCENARIO("Correct creation of UnitaryTableau") {
  GIVEN("An identity circuit") {
    UnitaryTableau tab(3);
    REQUIRE(tab.get_zrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::Z, 1.));
    REQUIRE(tab.get_zrow(Qubit(1)) == QubitPauliTensor(Qubit(1), Pauli::Z, 1.));
    REQUIRE(tab.get_zrow(Qubit(2)) == QubitPauliTensor(Qubit(2), Pauli::Z, 1.));
    REQUIRE(tab.get_xrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::X, 1.));
    REQUIRE(tab.get_xrow(Qubit(1)) == QubitPauliTensor(Qubit(1), Pauli::X, 1.));
    REQUIRE(tab.get_xrow(Qubit(2)) == QubitPauliTensor(Qubit(2), Pauli::X, 1.));
  }
  GIVEN("A single S gate") {
    UnitaryTableau tab0(3);
    UnitaryTableau tab1(3);
    tab0.apply_S_at_front(Qubit(0));
    tab1.apply_S_at_end(Qubit(0));
    REQUIRE(
        tab0.get_zrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::Z, 1.));
    REQUIRE(
        tab0.get_xrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::Y, -1.));
    REQUIRE(tab0 == tab1);
  }
  GIVEN("A single V gate") {
    UnitaryTableau tab0(3);
    UnitaryTableau tab1(3);
    tab0.apply_V_at_front(Qubit(0));
    tab1.apply_V_at_end(Qubit(0));
    REQUIRE(
        tab0.get_zrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::Y, 1.));
    REQUIRE(
        tab0.get_xrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::X, 1.));
    REQUIRE(tab0 == tab1);
  }
  GIVEN("A single H gate") {
    UnitaryTableau tab0(3);
    UnitaryTableau tab1(3);
    tab0.apply_gate_at_front(OpType::H, {Qubit(0)});
    tab1.apply_gate_at_end(OpType::H, {Qubit(0)});
    REQUIRE(
        tab0.get_zrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::X, 1.));
    REQUIRE(
        tab0.get_xrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::Z, 1.));
    REQUIRE(tab0 == tab1);
  }
  GIVEN("A single CX gate") {
    UnitaryTableau tab0(3);
    UnitaryTableau tab1(3);
    tab0.apply_CX_at_front(Qubit(0), Qubit(1));
    tab1.apply_CX_at_end(Qubit(0), Qubit(1));
    REQUIRE(
        tab0.get_zrow(Qubit(0)) == QubitPauliTensor(Qubit(0), Pauli::Z, 1.));
    REQUIRE(
        tab0.get_xrow(Qubit(1)) == QubitPauliTensor(Qubit(1), Pauli::X, 1.));
    REQUIRE(
        tab0.get_zrow(Qubit(1)) ==
        QubitPauliTensor(std::list<Pauli>({Pauli::Z, Pauli::Z})));
    REQUIRE(
        tab0.get_xrow(Qubit(0)) ==
        QubitPauliTensor(std::list<Pauli>({Pauli::X, Pauli::X})));
    REQUIRE(tab0 == tab1);
  }
  GIVEN("A Clifford circuit") {
    Circuit circ = get_test_circ();
    UnitaryTableau tab = circuit_to_unitary_tableau(circ);
    UnitaryTableau rev_tab = get_tableau_with_gates_applied_at_front();
    REQUIRE(tab == rev_tab);
  }
  GIVEN("A PI/2 rotation") {
    Circuit circ = get_test_circ();
    UnitaryTableau tab = circuit_to_unitary_tableau(circ);
    QubitPauliTensor pauli =
        QubitPauliTensor(std::list<Pauli>{Pauli::X, Pauli::Y, Pauli::Z});
    tab.apply_pauli_at_end(pauli, 3);

    add_ops_list_two_to_circuit(circ, OpType::Sdg);
    UnitaryTableau correct_tab = circuit_to_unitary_tableau(circ);
    REQUIRE(tab == correct_tab);
  }
  GIVEN("A PI/2 rotation at front") {
    UnitaryTableau tab = get_tableau_with_gates_applied_at_front();
    QubitPauliTensor pauli =
        QubitPauliTensor(std::list<Pauli>({Pauli::X, Pauli::Y, Pauli::Z}));
    tab.apply_pauli_at_front(pauli, 1);

    Circuit circ(3);
    add_ops_list_two_to_circuit(circ);
    add_ops_list_one_to_circuit(circ);
    UnitaryTableau correct_tab = circuit_to_unitary_tableau(circ);
    REQUIRE(tab == correct_tab);
  }
  GIVEN("Combining two circuits via tableau compose") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    UnitaryTableau first = circuit_to_unitary_tableau(circ);

    Circuit circ1(3);
    add_ops_list_two_to_circuit(circ1);
    UnitaryTableau second = circuit_to_unitary_tableau(circ1);
    UnitaryTableau correct = circuit_to_unitary_tableau(circ >> circ1);
    UnitaryTableau result = UnitaryTableau::compose(first, second);
    REQUIRE(result == correct);
  }
}

SCENARIO("Error handling in UnitaryTableau generation") {
  GIVEN("Add a non-clifford gate at end") {
    UnitaryTableau tab(2);
    REQUIRE_THROWS_AS(tab.apply_gate_at_end(OpType::T, {Qubit(0)}), NotValid);
  }
  GIVEN("Add a non-clifford gate at front") {
    UnitaryTableau tab(2);
    REQUIRE_THROWS_AS(
        tab.apply_gate_at_front(OpType::Tdg, {Qubit(0)}), NotValid);
  }
  GIVEN("Tableau from a non-Clifford circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CH, {1, 0});
    REQUIRE_THROWS_AS(circuit_to_unitary_tableau(circ), NotValid);
  }
}

SCENARIO("Synthesis of circuits from UnitaryTableau") {
  GIVEN("Transform a circuit in a round trip") {
    Circuit circ(3);
    add_ops_list_one_to_circuit(circ);
    UnitaryTableau tab = circuit_to_unitary_tableau(circ);
    Circuit res = unitary_tableau_to_circuit(tab);
    UnitaryTableau res_tab = circuit_to_unitary_tableau(res);
    REQUIRE(res_tab == tab);
  }
}

SCENARIO("UnitaryTableauBoxes in Circuits") {
  GIVEN("A tableau from another circuit") {
    Circuit inner(3);
    add_ops_list_one_to_circuit(inner);
    UnitaryTableau tab = circuit_to_unitary_tableau(inner);
    Circuit circ(4);
    Op_ptr box = std::make_shared<const UnitaryTableauBox>(tab);
    circ.add_op<unsigned>(OpType::CZ, uvec({1, 2}));
    circ.add_op<unsigned>(box, uvec({0, 1, 3}));
    circ.add_op<unsigned>(OpType::CX, uvec({0, 2}));
    Circuit correct(4);
    correct.add_op<unsigned>(OpType::CZ, uvec({1, 2}));
    correct.add_op<unsigned>(OpType::SWAP, uvec({2, 3}));
    add_ops_list_one_to_circuit(correct);
    correct.add_op<unsigned>(OpType::SWAP, uvec({2, 3}));
    correct.add_op<unsigned>(OpType::CX, uvec({0, 2}));
    REQUIRE(test_unitary_comparison(circ, correct, true));
  }
}

SCENARIO("Unitary inversions") {
  GIVEN("Some unitary tableau") {
    Circuit inner(3);
    add_ops_list_one_to_circuit(inner);
    UnitaryTableau tab = circuit_to_unitary_tableau(inner);
    Op_ptr box = std::make_shared<const UnitaryTableauBox>(tab);
    WHEN("Dagger") {
      Op_ptr box_dagger = box->dagger();
      Circuit circ(3);
      circ.add_op<unsigned>(box_dagger, uvec{0, 1, 2});
      REQUIRE(test_unitary_comparison(circ, inner.dagger(), true));
    }
    WHEN("Transpose") {
      Op_ptr box_transpose = box->transpose();
      Circuit circ(3);
      circ.add_op<unsigned>(box_transpose, uvec{0, 1, 2});
      REQUIRE(test_unitary_comparison(circ, inner.transpose(), true));
    }
  }
}

SCENARIO("Compare SymplecticTableau and UnitaryTableau") {
  GIVEN("The same sequence of gates, compare string representations") {
    SymplecticTableau stab(PauliStabiliserList{
        {{Pauli::Z, Pauli::I, Pauli::I}, true},
        {{Pauli::I, Pauli::Z, Pauli::I}, true},
        {{Pauli::I, Pauli::I, Pauli::Z}, true}});
    // Paulis cancel with subsequent gadget
    stab.apply_gate(OpType::X, uvec{0});
    stab.apply_gate(OpType::Y, uvec{1});
    stab.apply_gate(OpType::Z, uvec{2});
    stab.apply_pauli_gadget({{Pauli::X, Pauli::Y, Pauli::Z}, true}, 2);
    // CY and CZ combine to Sdg(0), CX(0, 1)
    stab.apply_gate(OpType::CY, uvec{0, 1});
    stab.apply_gate(OpType::CZ, uvec{0, 1});
    // SWAP that will remain
    stab.apply_gate(OpType::SWAP, uvec{1, 2});
    // BRIDGE cancels CX from CY+CZ
    stab.apply_gate(OpType::BRIDGE, uvec{0, 1, 2});
    std::stringstream stabstr;
    stabstr << stab;
    CHECK(stabstr.str() == "0 0 0 1 0 0 0\n0 0 0 0 0 1 0\n0 0 0 0 1 0 0\n");
    UnitaryTableau utab(3);
    // Same sequence, but appended to the front instead of the end
    utab.apply_gate_at_front(OpType::BRIDGE, {Qubit(0), Qubit(1), Qubit(2)});
    utab.apply_gate_at_front(OpType::SWAP, {Qubit(1), Qubit(2)});
    utab.apply_gate_at_front(OpType::CZ, {Qubit(0), Qubit(1)});
    utab.apply_gate_at_front(OpType::CY, {Qubit(0), Qubit(1)});
    utab.apply_pauli_at_front(
        QubitPauliTensor{{Pauli::X, Pauli::Y, Pauli::Z}}, 2);
    utab.apply_gate_at_front(OpType::X, {Qubit(0)});
    utab.apply_gate_at_front(OpType::Y, {Qubit(1)});
    utab.apply_gate_at_front(OpType::Z, {Qubit(2)});
    std::stringstream utabstr;
    utabstr << utab;
    CHECK(
        utabstr.str() ==
        "X@q[0]\t->\t1 0 0   1 0 0   0\n"
        "X@q[1]\t->\t0 0 1   0 0 0   0\n"
        "X@q[2]\t->\t0 1 0   0 0 0   0\n"
        "--\n"
        "Z@q[0]\t->\t0 0 0   1 0 0   0\n"
        "Z@q[1]\t->\t0 0 0   0 0 1   0\n"
        "Z@q[2]\t->\t0 0 0   0 1 0   0\n");
  }
}

SCENARIO("Tableau serialisation") {
  GIVEN("A circuit containing a tableau") {
    MatrixXb xx(3, 3);
    MatrixXb xz(3, 3);
    VectorXb xph(3);
    MatrixXb zx(3, 3);
    MatrixXb zz(3, 3);
    VectorXb zph(3);
    xx << 1, 1, 0, 0, 1, 0, 0, 0, 1;
    xz << 0, 0, 0, 0, 0, 0, 0, 0, 1;
    xph << 0, 0, 1;
    zx << 0, 0, 0, 0, 1, 0, 0, 0, 0;
    zz << 1, 0, 0, 1, 1, 0, 0, 0, 1;
    zph << 1, 0, 1;
    Op_ptr box =
        std::make_shared<const UnitaryTableauBox>(xx, xz, xph, zx, zz, zph);
    Circuit circ(3);
    circ.add_op<unsigned>(box, uvec{0, 1, 2});

    nlohmann::json j_circ = circ;
    Circuit circ2 = j_circ.get<Circuit>();
    REQUIRE(circ2 == circ);
  }
}

}  // namespace tket::test_UnitaryTableau
