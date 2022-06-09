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

#include "ArchAwareSynth/SteinerForest.hpp"
#include "testutil.hpp"
namespace tket {
SCENARIO("Synthesise a CNOT-only steiner Forest") {
  GIVEN("Empty circuit") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    Circuit circ(4);
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("One CNOT") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("Small circuit") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("Medium circuit") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(2), Node(4)}});
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
  }
}
SCENARIO("Synthesis a Rz-only Steiner Forest") {
  GIVEN("A single Rz") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("Three Rzs") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
  }
}

SCENARIO("Build a steiner Forest") {
  GIVEN("construction") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    REQUIRE(sf.tree_count == 1);
  }
  GIVEN("simple test") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 0);
    REQUIRE(sf.synth_circuit.n_vertices() == 7);
  }
  GIVEN("complex") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(2), Node(4)},
         {Node(2), Node(5)}});

    Circuit circ(6);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 2);
    REQUIRE(sf.synth_circuit.n_vertices() == 13);
  }
  GIVEN("complex 2") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)},
         {Node(7), Node(8)}});
    Circuit circ(9);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 8});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {8});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {5});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {7, 8});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 3);
    REQUIRE(sf.synth_circuit.n_vertices() == 18);
    REQUIRE(sf.synth_circuit.depth() == 0);
  }
  GIVEN("add_row_globally 1") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 1);
    REQUIRE(sf.synth_circuit.n_vertices() == 8);
    REQUIRE(sf.synth_circuit.depth() == 0);

    sf.add_row_globally(0, 1);
    sf.add_row_globally(2, 3);
    sf.add_row_globally(3, 2);
    sf.add_row_globally(1, 0);

    REQUIRE(sf.tree_count == 1);
    REQUIRE(sf.synth_circuit.n_vertices() == 12);
    REQUIRE(sf.synth_circuit.depth() == 2);
  }
  GIVEN("add_row_globally 2") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 1);

    REQUIRE(sf.synth_circuit.n_vertices() == 8);
    REQUIRE(sf.synth_circuit.depth() == 0);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);

    REQUIRE(sf.tree_count == 1);
    REQUIRE(sf.synth_circuit.n_vertices() == 17);
    REQUIRE(sf.synth_circuit.depth() == 9);
  }
  GIVEN("add_row_globally 3") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 1);

    REQUIRE(sf.synth_circuit.n_vertices() == 8);
    REQUIRE(sf.synth_circuit.depth() == 0);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);
    sf.add_row_globally(0, 1);
    sf.add_row_globally(1, 0);

    REQUIRE(sf.tree_count == 1);
    REQUIRE(sf.synth_circuit.n_vertices() == 17);
    REQUIRE(sf.synth_circuit.depth() == 9);
  }
  GIVEN("add_operation_list 1") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 1);

    REQUIRE(sf.synth_circuit.n_vertices() == 8);
    REQUIRE(sf.synth_circuit.depth() == 0);

    aas::OperationList oplist = {
        std::pair(1, 0), std::pair(0, 1), std::pair(1, 0), std::pair(0, 1),
        std::pair(1, 0), std::pair(0, 1), std::pair(1, 0), std::pair(0, 1)};

    sf.add_operation_list(oplist);

    REQUIRE(sf.tree_count == 1);
    REQUIRE(sf.synth_circuit.n_vertices() == 16);
    REQUIRE(sf.synth_circuit.depth() == 8);
  }
  GIVEN("recursive_operation_search") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 1);

    REQUIRE(sf.synth_circuit.n_vertices() == 8);
    REQUIRE(sf.synth_circuit.depth() == 0);

    aas::PathHandler pathhand(archi);

    aas::OperationList oplist = sf.operations_available_at_index(pathhand, 3);
    aas::OperationList oplist2 = {std::pair(1, 0), std::pair(2, 3)};
    REQUIRE(oplist == oplist2);

    aas::CostedOperations cosop =
        recursive_operation_search(pathhand, sf, 2, oplist);

    aas::CostedOperations expectedResult = std::pair(2, oplist2);
    REQUIRE(cosop == expectedResult);
  }
  GIVEN("operations_available_at_index") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);

    REQUIRE(sf.tree_count == 1);

    REQUIRE(sf.synth_circuit.n_vertices() == 8);
    REQUIRE(sf.synth_circuit.depth() == 0);

    aas::PathHandler pathhand(archi);

    aas::OperationList oplist = sf.operations_available_at_index(pathhand, 3);
    aas::OperationList oplist2 = {std::pair(1, 0), std::pair(2, 3)};

    REQUIRE(oplist == oplist2);
  }
}
SCENARIO("check error in steiner Forest") {
  GIVEN("lookahead 0 - phase_poly_synthesis_int") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    REQUIRE_THROWS_AS(
        aas::phase_poly_synthesis_int(archi, ppbox, 0, aas::CNotSynthType::Rec),
        std::logic_error);
  }
  GIVEN("lookahead 0 - best_operations_lookahead") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    aas::PathHandler ph = aas::PathHandler(archi);
    REQUIRE_THROWS_AS(
        aas::best_operations_lookahead(ph, sf, 0), std::logic_error);
  }
  GIVEN("empty steinerforest") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    Circuit circ(2);
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    aas::PathHandler ph = aas::PathHandler(archi);
    REQUIRE_THROWS_AS(
        aas::best_operations_lookahead(ph, sf, 1), std::logic_error);
  }
}
SCENARIO("Synthesise a phase polynomial for a given architecture", "[.long]") {
  GIVEN("phase_poly_synthesis 1") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 2);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 8);
    REQUIRE(result.count_gates(OpType::CX) == 7);
  }
  GIVEN("phase_poly_synthesis 1 - swap") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    PhasePolyBox ppbox(circ);
    Circuit result =
        aas::phase_poly_synthesis(archi, ppbox, 2, aas::CNotSynthType::SWAP);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 8);
    REQUIRE(result.count_gates(OpType::CX) == 7);
  }
  GIVEN("phase_poly_synthesis 1 - hampath") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    PhasePolyBox ppbox(circ);
    Circuit result =
        phase_poly_synthesis(archi, ppbox, 2, aas::CNotSynthType::HamPath);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 6);
    REQUIRE(result.count_gates(OpType::CX) == 5);
  }
  GIVEN("phase_poly_synthesis 1 - rec") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    PhasePolyBox ppbox(circ);
    Circuit result =
        aas::phase_poly_synthesis(archi, ppbox, 2, aas::CNotSynthType::Rec);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 8);
    REQUIRE(result.count_gates(OpType::CX) == 7);
  }
  GIVEN("phase_poly_synthesis 2") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 2);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 11);
    REQUIRE(result.count_gates(OpType::CX) == 10);
  }
  GIVEN("phase_poly_synthesis 3") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::Rz, 0.11, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.12, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.13, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.14, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.6, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 12);
    REQUIRE(result.count_gates(OpType::CX) == 6);
  }
  GIVEN("phase_poly_synthesis 4") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PhasePolyBox ppbox(circ);
    aas::SteinerForest sf = aas::SteinerForest(archi, ppbox);
    REQUIRE(sf.tree_count == 0);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 1);
    REQUIRE(result.count_gates(OpType::CX) == 0);
  }
  GIVEN("phase_poly_synthesis 5") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(8);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {7});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 2);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 20);
    REQUIRE(result.count_gates(OpType::CX) == 19);
  }
  GIVEN("phase_poly_synthesis 6") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(8);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {5});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {6});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {7});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 2);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 14);
    REQUIRE(result.count_gates(OpType::CX) == 7);
  }
  GIVEN("phase_poly_synthesis 7") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)},
         {Node(7), Node(8)}});
    Circuit circ(9);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 8});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {8});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {5});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {7, 8});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 19);
    REQUIRE(result.count_gates(OpType::CX) == 16);
  }
  GIVEN("phase_poly_synthesis 8") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 7);
    REQUIRE(result.count_gates(OpType::CX) == 6);
  }
  GIVEN("phase_poly_synthesis 9") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 5);
    REQUIRE(result.count_gates(OpType::CX) == 4);
  }
  GIVEN("phase_poly_synthesis 10") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 2);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 10);
    REQUIRE(result.count_gates(OpType::CX) == 8);
  }
  GIVEN("phase_poly_synthesis 11") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 2});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 7);
    REQUIRE(result.count_gates(OpType::CX) == 6);
  }
  GIVEN("phase_poly_synthesis 12") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(3)}, {Node(3), Node(2)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 12);
    REQUIRE(result.count_gates(OpType::CX) == 10);
  }
  GIVEN("phase_poly_synthesis 13") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(8);
    circ.add_op<unsigned>(OpType::CX, {5, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
    circ.add_op<unsigned>(OpType::CX, {2, 7});
    circ.add_op<unsigned>(OpType::Rz, 0.6, {7});
    circ.add_op<unsigned>(OpType::CX, {7, 6});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {6});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 55);
    REQUIRE(result.count_gates(OpType::CX) == 48);
  }
  GIVEN("phase_poly_synthesis 14") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(2), Node(3)},
         {Node(3), Node(4)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 2);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 23);
    REQUIRE(result.count_gates(OpType::CX) == 21);
  }
  GIVEN("phase_poly_synthesis 15") {
    const Architecture archi({{Node(0), Node(1)}, {Node(1), Node(2)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 2});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 7);
    REQUIRE(result.count_gates(OpType::CX) == 6);
  }
  GIVEN("phase_poly_synthesis 16") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(3)}, {Node(3), Node(2)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 12);
    REQUIRE(result.count_gates(OpType::CX) == 10);
  }
  GIVEN("phase_poly_synthesis 17") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(0)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {4});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 11);
    REQUIRE(result.count_gates(OpType::CX) == 10);
  }
  GIVEN("phase_poly_synthesis 18") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(5), Node(7)},
         {Node(6), Node(7)},
         {Node(7), Node(8)},
         {Node(8), Node(9)},
         {Node(9), Node(0)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(10);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});

    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 5});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {6});

    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 9});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {9});

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 40);
    REQUIRE(result.count_gates(OpType::CX) == 36);
  }
  GIVEN("phase_poly_synthesis 19") {
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});

    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 8);
    REQUIRE(result.count_gates(OpType::CX) == 7);
  }
  GIVEN("phase_poly_synthesis 20") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(0)}});

    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 8);
    REQUIRE(result.count_gates(OpType::CX) == 7);
  }
  GIVEN("phase_poly_synthesis 21") {
    // without labelled hamilton path
    const Architecture archi(
        {{Node(10), Node(12)},
         {Node(9), Node(6)},
         {Node(12), Node(6)},
         {Node(9), Node(10)}});

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);

    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 16);
    REQUIRE(result.count_gates(OpType::CX) == 15);
  }
  GIVEN("phase_poly_synthesis 22") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(6)},
         {Node(2), Node(3)},
         {Node(2), Node(5)},
         {Node(2), Node(7)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(8);
    circ.add_op<unsigned>(OpType::CX, {5, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
    circ.add_op<unsigned>(OpType::CX, {2, 7});
    circ.add_op<unsigned>(OpType::Rz, 0.6, {7});
    circ.add_op<unsigned>(OpType::CX, {7, 6});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {6});

    PhasePolyBox ppbox(circ);
    Circuit result =
        aas::phase_poly_synthesis(archi, ppbox, 1, aas::CNotSynthType::Rec);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 88);
    REQUIRE(result.count_gates(OpType::CX) == 81);
  }
  GIVEN("phase_poly_synthesis 23") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(9)}});

    Circuit circ(10);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {4});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 77);
    REQUIRE(result.count_gates(OpType::CX) == 76);
  }
  GIVEN("phase_poly_synthesis 24") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(9)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(10);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {4});

    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {6, 1});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 19);
    REQUIRE(result.count_gates(OpType::CX) == 18);
  }
  GIVEN("phase_poly_synthesis 25") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)}, {Node(0), Node(2)}, {Node(0), Node(3)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {3});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 21);
    REQUIRE(result.count_gates(OpType::CX) == 20);
  }
  GIVEN("phase_poly_synthesis 26") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(6), Node(7)},
         {Node(7), Node(8)},
         {Node(0), Node(3)},
         {Node(3), Node(6)},
         {Node(1), Node(4)},
         {Node(4), Node(7)},
         {Node(2), Node(5)},
         {Node(5), Node(8)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(9);
    circ.add_op<unsigned>(OpType::CX, {5, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::CX, {3, 0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
    circ.add_op<unsigned>(OpType::CX, {2, 7});
    circ.add_op<unsigned>(OpType::Rz, 0.6, {7});
    circ.add_op<unsigned>(OpType::CX, {7, 6});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {6});
    circ.add_op<unsigned>(OpType::CX, {6, 8});
    circ.add_op<unsigned>(OpType::Rz, 0.7, {8});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 154);
    REQUIRE(result.count_gates(OpType::CX) == 146);
  }
  GIVEN("phase_poly_synthesis 27") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(9)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(10);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {4});

    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {6, 1});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 19);
    REQUIRE(result.count_gates(OpType::CX) == 18);
  }
  GIVEN("phase_poly_synthesis 27") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(15);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 10});
    circ.add_op<unsigned>(OpType::CX, {10, 11});
    circ.add_op<unsigned>(OpType::CX, {11, 12});
    circ.add_op<unsigned>(OpType::CX, {12, 13});
    circ.add_op<unsigned>(OpType::CX, {13, 14});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {14});

    PhasePolyBox ppbox(circ);
    Circuit result = aas::phase_poly_synthesis(archi, ppbox, 1);
    // the following line is not executed for performance reasons
    // REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 329);
    REQUIRE(result.count_gates(OpType::CX) == 328);
  }
  GIVEN("phase_poly_synthesis 27 - swap") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(15);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 10});
    circ.add_op<unsigned>(OpType::CX, {10, 11});
    circ.add_op<unsigned>(OpType::CX, {11, 12});
    circ.add_op<unsigned>(OpType::CX, {12, 13});
    circ.add_op<unsigned>(OpType::CX, {13, 14});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {14});

    PhasePolyBox ppbox(circ);
    Circuit result =
        aas::phase_poly_synthesis(archi, ppbox, 1, aas::CNotSynthType::SWAP);
    // the following line is not executed for performance reasons
    // REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 1156);
    REQUIRE(result.count_gates(OpType::CX) == 1155);
  }
  GIVEN("phase_poly_synthesis 27 - ham path") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(15);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 10});
    circ.add_op<unsigned>(OpType::CX, {10, 11});
    circ.add_op<unsigned>(OpType::CX, {11, 12});
    circ.add_op<unsigned>(OpType::CX, {12, 13});
    circ.add_op<unsigned>(OpType::CX, {13, 14});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {14});

    PhasePolyBox ppbox(circ);
    REQUIRE_THROWS_AS(
        phase_poly_synthesis(archi, ppbox, 1, aas::CNotSynthType::HamPath),
        std::logic_error);
  }
  GIVEN("phase_poly_synthesis 27 - loop") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(15);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {4, 10});
    circ.add_op<unsigned>(OpType::CX, {10, 11});
    circ.add_op<unsigned>(OpType::CX, {11, 12});
    circ.add_op<unsigned>(OpType::CX, {12, 13});
    circ.add_op<unsigned>(OpType::CX, {13, 14});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {14});

    PhasePolyBox ppbox(circ);
    Circuit result =
        aas::phase_poly_synthesis(archi, ppbox, 1, aas::CNotSynthType::Rec);
    // the following line is not executed for performance reasons
    // REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 329);
    REQUIRE(result.count_gates(OpType::CX) == 328);
  }
  GIVEN("phase_poly_synthesis 24 - swap") {
    // No Hamiltonian path
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(9)}});
    aas::PathHandler pathhand(archi);

    Circuit circ(10);
    circ.add_op<unsigned>(OpType::CX, {6, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {7, 4});

    circ.add_op<unsigned>(OpType::Rz, 0.7, {4});

    circ.add_op<unsigned>(OpType::CX, {7, 4});
    circ.add_op<unsigned>(OpType::CX, {9, 7});
    circ.add_op<unsigned>(OpType::CX, {3, 9});
    circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::CX, {2, 5});
    circ.add_op<unsigned>(OpType::CX, {8, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 8});
    circ.add_op<unsigned>(OpType::CX, {6, 1});

    PhasePolyBox ppbox(circ);
    Circuit result =
        aas::phase_poly_synthesis(archi, ppbox, 1, aas::CNotSynthType::SWAP);
    REQUIRE(test_unitary_comparison(circ, result));
    REQUIRE(result.n_gates() == 19);
    REQUIRE(result.count_gates(OpType::CX) == 18);
  }
}
}  // namespace tket
