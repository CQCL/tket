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
#include <fstream>

#include "../testutil.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/PauliGraphRefactor/Converters.hpp"
#include "tket/PauliGraphRefactor/PauliGraph.hpp"

namespace tket {
namespace test_PauliGraph3 {

using namespace pg;

bool comp_seqs(
    const std::list<PGOp_ptr>& seq1, const std::list<PGOp_ptr>& seq2) {
  if (seq1.size() != seq2.size()) return false;
  std::list<PGOp_ptr>::const_iterator it2 = seq2.begin();
  for (const PGOp_ptr& op1 : seq1) {
    if (*op1 != **it2) return false;
    ++it2;
  }
  return true;
}

SCENARIO("Correct creation of refactored PauliGraphs") {
  GIVEN("A Clifford circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Vdg, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 1qb circuit") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with no interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rx, 1.6, {1});
    circ.add_op<unsigned>(OpType::Rz, 1.3, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with some anti-commuting interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {1});
    circ.add_op<unsigned>(OpType::XXPhase, 1.1, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with some commuting interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.1, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit a Clifford-angled ZZPhase") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 0.5, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 1qb circuit with stuff to merge") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A 2qb circuit with stuff to merge") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.1, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.8, {0});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.6, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A circuit with Cliffords and non-Cliffords") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {0});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.8, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A dense example") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {1});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {1});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {3});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A more interesting example (tof_3)") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::T, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Tdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::Tdg, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::Tdg, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::T, {2});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::T, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Tdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("A circuit with a PauliExpBox") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::ZZPhase, 0.2, {0, 1});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    PauliExpBox peb(SymPauliTensor({Pauli::Y, Pauli::X}, 0.333));
    circ.add_box(peb, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("Teleportation") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_measure(0, 0);
    circ.add_measure(1, 1);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{2}, {1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Z, {}, uvec{2}, {0}, 1);
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    std::list<PGOp_ptr> sequence = pg.pgop_sequence();
    std::list<PGOp_ptr> correct_sequence{
        std::make_shared<PGInputTableau>(ChoiMixTableau(3)),
        std::make_shared<PGMeasure>(
            SpPauliStabiliser({Pauli::Z, Pauli::X, Pauli::I}), Bit(1)),
        std::make_shared<PGMeasure>(
            SpPauliStabiliser({Pauli::X, Pauli::Z, Pauli::X}), Bit(0)),
        std::make_shared<PGConditional>(
            std::make_shared<PGCliffordRot>(
                SpPauliStabiliser({Pauli::I, Pauli::I, Pauli::X}), 2),
            bit_vector_t{Bit(1)}, 1),
        std::make_shared<PGConditional>(
            std::make_shared<PGCliffordRot>(
                SpPauliStabiliser({Pauli::I, Pauli::X, Pauli::Z}), 2),
            bit_vector_t{Bit(0)}, 1),
        std::make_shared<PGOutputTableau>(ChoiMixTableau({
            {SpPauliStabiliser({Pauli::X, Pauli::Z, Pauli::X}),
             SpPauliStabiliser(Qubit(0), Pauli::Z)},
            {SpPauliStabiliser({Pauli::Z, Pauli::I, Pauli::I}),
             SpPauliStabiliser(Qubit(0), Pauli::X)},
            {SpPauliStabiliser({Pauli::Z, Pauli::X, Pauli::I}),
             SpPauliStabiliser(Qubit(1), Pauli::Z)},
            {SpPauliStabiliser({Pauli::I, Pauli::Z, Pauli::X}),
             SpPauliStabiliser(Qubit(1), Pauli::X)},
            {SpPauliStabiliser({Pauli::I, Pauli::X, Pauli::Z}),
             SpPauliStabiliser(Qubit(2), Pauli::Z)},
            {SpPauliStabiliser({Pauli::I, Pauli::I, Pauli::X}),
             SpPauliStabiliser(Qubit(2), Pauli::X)},
        }))};
    REQUIRE_NOTHROW(pg.verify());
    CHECK(comp_seqs(sequence, correct_sequence));
    THEN("Print diagram to file") {
      std::ofstream dot_file("pauligraph.dot");
      pg.to_graphviz(dot_file);
      dot_file.close();
      remove("pauligraph.dot");
    }
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    PauliGraph res_pg = circuit_to_pauli_graph3(res);
    std::list<PGOp_ptr> res_sequence = res_pg.pgop_sequence();
    CHECK(comp_seqs(res_sequence, correct_sequence));
  }
  GIVEN("A conjugated Reset and Collapse") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::Reset, {1});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::Collapse, {2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(res.count_gates(OpType::Reset) == 1);
    ChoiMixTableau circ_tab = circuit_to_cm_tableau(circ);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(circ_tab == res_tab);
  }
  GIVEN("A conjugated Box") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::Sycamore, {1, 2});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(res.count_gates(OpType::Sycamore) == 1);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
  GIVEN("Some stabiliser assertions") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 1.5, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliStabiliser pauli1 = {{Pauli::X, Pauli::X}, 0};
    PauliStabiliser pauli2 = {{Pauli::Z, Pauli::Z}, 0};
    PauliStabiliser pauli3 = {{Pauli::Y, Pauli::Y}, 2};
    PauliStabiliserVec stabilisers = {pauli1, pauli2, pauli3};
    StabiliserAssertionBox box(stabilisers);
    circ.add_assertion(box, {Qubit(0), Qubit(2)}, Qubit(1));
    circ.add_assertion(box, {Qubit(0), Qubit(2)}, Qubit(1));
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    std::list<PGOp_ptr> sequence = pg.pgop_sequence();
    SpPauliStabiliser anc_z(Qubit(1), Pauli::Z);
    SpPauliStabiliser anc_x(DensePauliMap{Pauli::X, Pauli::X});
    std::list<PGOp_ptr> correct_sequence{
        std::make_shared<PGInputTableau>(ChoiMixTableau(3)),
        std::make_shared<PGRotation>(
            SpPauliStabiliser(Qubit(0), Pauli::Z), 1.5),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::X, Pauli::I, Pauli::X}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name(), 0)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::Z}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name(), 1)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Y, Pauli::Z, Pauli::Y}, 2), anc_z, anc_x,
            Bit(c_debug_one_prefix() + "_" + c_debug_default_name(), 0)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::X, Pauli::I, Pauli::X}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)",
                0)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::Z}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)",
                1)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Y, Pauli::Z, Pauli::Y}, 2), anc_z, anc_x,
            Bit(c_debug_one_prefix() + "_" + c_debug_default_name() + "(1)",
                0)),
        std::make_shared<PGOutputTableau>(ChoiMixTableau({
            {SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::I}),
             SpPauliStabiliser(Qubit(0), Pauli::Z)},
            {SpPauliStabiliser(Qubit(0), Pauli::X),
             SpPauliStabiliser(Qubit(0), Pauli::X)},
            {SpPauliStabiliser(Qubit(1), Pauli::Z),
             SpPauliStabiliser(Qubit(1), Pauli::Z)},
            {SpPauliStabiliser({Pauli::X, Pauli::X, Pauli::I}),
             SpPauliStabiliser(Qubit(1), Pauli::X)},
            {SpPauliStabiliser(Qubit(2), Pauli::Z),
             SpPauliStabiliser(Qubit(2), Pauli::Z)},
            {SpPauliStabiliser(Qubit(2), Pauli::X),
             SpPauliStabiliser(Qubit(2), Pauli::X)},
        }))};
    CHECK(comp_seqs(sequence, correct_sequence));
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(res.count_gates(OpType::StabiliserAssertionBox) == 6);
    PauliGraph res_pg = circuit_to_pauli_graph3(res);
    std::list<PGOp_ptr> res_sequence = res_pg.pgop_sequence();
    CHECK(comp_seqs(res_sequence, correct_sequence));
  }
  GIVEN("Don't collect cliffords") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Y, {0});
    circ.add_op<unsigned>(OpType::Sdg, {1});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CY, {2, 0});
    circ.add_op<unsigned>(OpType::PhaseGadget, 0.198, {0, 1, 2});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::Vdg, {2});
    circ.add_op<unsigned>(OpType::CZ, {1, 2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {0, 2});
    circ.add_op<unsigned>(OpType::YYPhase, 1.387, {0, 1});
    circ.add_op<unsigned>(OpType::TK1, {0.98, 0.2, 1.87}, {1});
    circ.add_op<unsigned>(OpType::TK2, {1.34, 0.23, 1.42}, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ, false);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(test_unitary_comparison(circ, res, true));
  }
}

}  // namespace test_PauliGraph3
}  // namespace tket
