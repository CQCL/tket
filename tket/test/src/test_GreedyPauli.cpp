// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testutil.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Gate/SymTable.hpp"
#include "tket/Ops/ClassicalOps.hpp"
#include "tket/Predicates/PassGenerators.hpp"
#include "tket/Transformations/GreedyPauliOptimisation.hpp"
#include "tket/Utils/Expression.hpp"

namespace tket {
namespace test_GreedyPauliSimp {

SCENARIO("Clifford synthesis") {
  GIVEN("Empty circuit") {
    Circuit circ(3);
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("1Q Simple Clifford") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Sdg, {0});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("2Q Simple Clifford") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Y, {0});
    circ.add_op<unsigned>(OpType::Vdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("3Q Simple Clifford") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Y, {0});
    circ.add_op<unsigned>(OpType::Sdg, {2});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {0, 2});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("5Q Simple Clifford") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CZ, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::X, {3});
    circ.add_op<unsigned>(OpType::CZ, {0, 4});
    circ.add_op<unsigned>(OpType::CY, {0, 1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::Y, {4});
    circ.add_op<unsigned>(OpType::CY, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 0});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("Clifford with swaps") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CZ, {1, 3});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::SWAP, {3, 1});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
}
SCENARIO("Complete synthesis") {
  GIVEN("1Q Simple Circuit") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Sdg, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {0});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("Symbolic Circuit") {
    Circuit circ(2);
    auto a = SymTable::fresh_symbol("a");
    auto b = SymTable::fresh_symbol("b");
    auto ea = Expr(a);
    auto eb = Expr(b);
    circ.add_op<unsigned>(OpType::Sdg, {0});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::Ry, eb, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rx, ea, {0});
    Circuit d(circ);
    symbol_map_t symbol_map;
    symbol_map[a] = Expr(0.5);
    symbol_map[b] = Expr(0.7);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    circ.symbol_substitution(symbol_map);
    d.symbol_substitution(symbol_map);
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("4Q PauliExp Circuit") {
    Circuit circ(4);
    circ.add_box(
        PauliExpBox(SymPauliTensor({Pauli::X, Pauli::X}, 0.3)), {0, 1});
    circ.add_box(
        PauliExpBox(SymPauliTensor({Pauli::Z, Pauli::Y}, -0.1)), {2, 3});
    circ.add_box(
        PauliExpPairBox(
            SymPauliTensor({Pauli::X, Pauli::Z}, 1.0),
            SymPauliTensor({Pauli::Z, Pauli::X}, 0.4)),
        {0, 2});
    circ.add_box(
        PauliExpCommutingSetBox({
            {{Pauli::I, Pauli::Y, Pauli::I}, -0.1},
            {{Pauli::X, Pauli::Y, Pauli::Z}, -1.2},
            {{Pauli::X, Pauli::Y, Pauli::Z}, 0.5},
        }),
        {1, 2, 3});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CZ, {1, 3});
    Circuit d = Transforms::GreedyPauliSimp::greedy_pauli_graph_synthesis(circ);
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("Arbitrary Circuit") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.15, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {4});
    circ.add_op<unsigned>(OpType::CZ, {1, 4});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::ZZPhase, 0.7, {3, 2});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::SWAP, {3, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("Circuit with trivial Pauli exps") {
    Circuit circ(4);
    circ.add_box(PauliExpBox(SymPauliTensor({Pauli::X, Pauli::X}, 2)), {0, 1});
    circ.add_box(
        PauliExpPairBox(
            SymPauliTensor({Pauli::I, Pauli::I}, 1.2),
            SymPauliTensor({Pauli::Z, Pauli::X}, -2)),
        {0, 2});
    circ.add_box(
        PauliExpCommutingSetBox({
            {{Pauli::I, Pauli::Y, Pauli::I}, 0},
            {{Pauli::X, Pauli::Y, Pauli::Z}, 0},
            {{Pauli::I, Pauli::I, Pauli::I}, 0.5},
        }),
        {1, 2, 3});
    Circuit d = Transforms::GreedyPauliSimp::greedy_pauli_graph_synthesis(circ);
    REQUIRE(test_unitary_comparison(circ, d, true));
    REQUIRE(d.n_gates() == 0);
  }
  GIVEN("Circuit with non-default UnitIDs") {
    Circuit circ;
    register_t reg_a = circ.add_q_register("a", 2);
    register_t reg_b = circ.add_q_register("b", 2);
    circ.add_op<UnitID>(OpType::CX, {reg_a[0], reg_b[1]});
    circ.add_op<UnitID>(OpType::SWAP, {reg_b[0], reg_a[1]});
    circ.add_op<UnitID>(OpType::Rz, 0.3, {reg_a[1]});
    circ.add_op<UnitID>(OpType::CX, {reg_a[1], reg_b[1]});
    circ.add_op<UnitID>(OpType::Ry, 0.2, {reg_b[1]});
    circ.add_op<UnitID>(OpType::H, {reg_b[1]});
    circ.add_op<UnitID>(OpType::Rz, 0.3, {reg_a[0]});
    circ.add_op<UnitID>(OpType::CY, {reg_a[0], reg_a[1]});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  GIVEN("Circuit with measurements") {
    Circuit circ(4, 4);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.15, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    // For circuit d, we add measurement after synthesis
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(d));
    REQUIRE(d.has_implicit_wireswaps());
    for (unsigned i = 0; i < 4; i++) {
      d.add_op<UnitID>(OpType::Measure, {Qubit(i), Bit(i)});
    }
    // For circuit g, we add measurement before synthesis
    Circuit g(circ);
    for (unsigned i = 0; i < 4; i++) {
      g.add_op<UnitID>(OpType::Measure, {Qubit(i), Bit(i)});
    }
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(g));
    REQUIRE(d == g);
  }
  GIVEN("Circuit with conditional gates") {
    Circuit circ(2, 2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_conditional_gate<unsigned>(OpType::Rz, {0.5}, {0}, {0}, 0);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit d(2, 2);
    Op_ptr cond = std::make_shared<Conditional>(
        std::make_shared<PauliExpBox>(
            SymPauliTensor({Pauli::Z, Pauli::I}, 0.5)),
        1, 0);
    d.add_conditional_gate<unsigned>(OpType::Sdg, {}, {0}, {0}, 0);
    d.add_conditional_gate<unsigned>(OpType::Z, {}, {0}, {0}, 0);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with conditional gates 2") {
    Circuit circ(2, 1);
    circ.add_box(
        PauliExpBox(SymPauliTensor({Pauli::X, Pauli::Z}, 0.12)), {0, 1});
    Op_ptr cond = std::make_shared<Conditional>(
        std::make_shared<PauliExpBox>(
            SymPauliTensor({Pauli::Z, Pauli::Z}, 0.5)),
        1, 0);
    circ.add_op<unsigned>(cond, {0, 0, 1});
    // two boxes anti-commute hence simultaneous diagonalisation
    Circuit d(2, 1);
    d.add_op<unsigned>(OpType::CY, {1, 0});
    d.add_op<unsigned>(OpType::Rx, 0.12, {0});
    d.add_conditional_gate<unsigned>(OpType::Sdg, {}, {0}, {0}, 0);
    d.add_conditional_gate<unsigned>(OpType::Z, {}, {0}, {0}, 0);
    d.add_op<unsigned>(OpType::CY, {1, 0});
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with conditional gates and measures") {
    Circuit circ(2, 2);
    Op_ptr cond1 = std::make_shared<Conditional>(
        std::make_shared<PauliExpBox>(
            SymPauliTensor({Pauli::Z, Pauli::X}, 0.5)),
        1, 0);
    Op_ptr cond2 = std::make_shared<Conditional>(
        std::make_shared<PauliExpBox>(
            SymPauliTensor({Pauli::Z, Pauli::Y}, 0.12)),
        1, 0);
    circ.add_op<unsigned>(cond1, {0, 0, 1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    // can commute to the front
    circ.add_op<unsigned>(cond2, {1, 0, 1});
    Circuit d(2, 2);
    // ZX 0.5
    d.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    d.add_conditional_gate<unsigned>(OpType::Sdg, {}, {0}, {0}, 0);
    d.add_conditional_gate<unsigned>(OpType::V, {}, {1}, {0}, 0);
    d.add_conditional_gate<unsigned>(OpType::Z, {}, {0}, {0}, 0);
    d.add_op<unsigned>(OpType::Measure, {0, 0});
    // ZY 0.12
    // compute
    d.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {1}, 0);
    d.add_conditional_gate<unsigned>(OpType::CY, {}, {0, 1}, {1}, 0);
    d.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {1}, 0);
    // action
    d.add_conditional_gate<unsigned>(OpType::Rz, {0.12}, {0}, {1}, 0);
    // uncompute
    d.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {1}, 0);
    d.add_conditional_gate<unsigned>(OpType::CY, {}, {0, 1}, {1}, 0);
    d.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {1}, 0);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }

  GIVEN("Conditionals merging") {
    Circuit circ(2, 2);
    Op_ptr cond1 = std::make_shared<Conditional>(
        std::make_shared<PauliExpBox>(
            SymPauliTensor({Pauli::Z, Pauli::X}, 0.25)),
        1, 0);
    Op_ptr cond2 = std::make_shared<Conditional>(
        std::make_shared<PauliExpBox>(
            SymPauliTensor({Pauli::Z, Pauli::X}, -0.25)),
        1, 0);
    circ.add_op<unsigned>(cond1, {0, 0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(cond2, {0, 0, 1});
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_op<unsigned>(OpType::Rz, -0.3, {0});
    // should all be canceled
    Circuit d(2, 2);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with classical gates") {
    Circuit circ(1, 4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(ClassicalX(), {1});
    circ.add_op<unsigned>(ClassicalCX(), {0, 1});
    circ.add_op<unsigned>(AndWithOp(), {2, 3});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with WASMs") {
    std::string wasm_file = "string/with/path/to/wasm/file";
    std::string wasm_func = "stringNameOfWASMFunc";
    std::vector<unsigned> uv = {2, 1};
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(6, 1, uv, uv, wasm_func, wasm_file);
    Circuit circ(1, 7);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<UnitID>(
        wop_ptr,
        {Bit(0), Bit(1), Bit(2), Bit(3), Bit(4), Bit(5), WasmState(0)});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with mid-circuit measurements") {
    Circuit circ(2, 2);
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_op<unsigned>(OpType::Tdg, {0});
    circ.add_op<unsigned>(OpType::Measure, {1, 1});
    Circuit d(2, 2);
    d.add_op<unsigned>(OpType::Measure, {0, 0});
    d.add_op<unsigned>(OpType::Measure, {1, 1});
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with resets") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Reset, {0});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
  GIVEN("Circuit with measures, classicals, and resets") {
    Circuit circ(3, 1);
    circ.add_box(
        PauliExpBox(SymPauliTensor({Pauli::X, Pauli::Z, Pauli::Z}, 0.3)),
        {0, 1, 2});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_op<unsigned>(ClassicalX(), {0});
    circ.add_op<unsigned>(OpType::Reset, {1});
    Circuit d(3, 1);
    d.add_op<unsigned>(OpType::CZ, {0, 2});
    d.add_op<unsigned>(OpType::CZ, {0, 1});
    d.add_op<unsigned>(OpType::Rx, 0.3, {0});
    d.add_op<unsigned>(OpType::Measure, {0, 0});
    d.add_op<unsigned>(ClassicalX(), {0});
    d.add_op<unsigned>(OpType::CZ, {1, 0});
    d.add_op<unsigned>(OpType::CZ, {0, 2});
    d.add_op<unsigned>(OpType::Reset, {1});
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(circ == d);
  }
}
SCENARIO("Test GreedyPauliSimp for individual gates") {
  Circuit circ(1);
  circ.add_op<unsigned>(OpType::Z, {0});
  std::vector<Op_ptr> ops_0q = {
      get_op_ptr(OpType::Phase, 0.25),
  };
  std::vector<Op_ptr> ops_1q = {
      get_op_ptr(OpType::noop),
      get_op_ptr(OpType::Z),
      get_op_ptr(OpType::X),
      get_op_ptr(OpType::Y),
      get_op_ptr(OpType::S),
      get_op_ptr(OpType::V),
      get_op_ptr(OpType::Sdg),
      get_op_ptr(OpType::Vdg),
      get_op_ptr(OpType::H),
      get_op_ptr(OpType::Rz, 0.25),
      get_op_ptr(OpType::Rz, 0.5),
      get_op_ptr(OpType::Rx, 1),
      get_op_ptr(OpType::Rx, 0.15),
      get_op_ptr(OpType::Ry, 0.25),
      get_op_ptr(OpType::Ry, -0.5),
      get_op_ptr(OpType::PhasedX, {0.15, 0.2}),
      get_op_ptr(OpType::PhasedX, {0.5, -0.5}),
      get_op_ptr(OpType::PhasedX, {0.2, 1}),
      get_op_ptr(OpType::T),
      get_op_ptr(OpType::Tdg),
  };
  std::vector<Op_ptr> ops_2q = {
      get_op_ptr(OpType::SWAP),
      get_op_ptr(OpType::CX),
      get_op_ptr(OpType::CY),
      get_op_ptr(OpType::CZ),
      get_op_ptr(OpType::ZZMax),
      get_op_ptr(OpType::ZZPhase, 0.25),
      get_op_ptr(OpType::ZZPhase, 0.5),
      get_op_ptr(OpType::PhaseGadget, 0.5, 2),
      get_op_ptr(OpType::XXPhase, 0.25),
      get_op_ptr(OpType::XXPhase, 0.5),
      get_op_ptr(OpType::YYPhase, 0.25),
      get_op_ptr(OpType::YYPhase, 1),
  };
  for (Op_ptr op : ops_0q) {
    Circuit circ(1);
    circ.add_op<unsigned>(op, {});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  for (Op_ptr op : ops_1q) {
    Circuit circ(1);
    circ.add_op<unsigned>(op, {0});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
  for (Op_ptr op : ops_2q) {
    Circuit circ(2);
    circ.add_op<unsigned>(op, {0, 1});
    Circuit d(circ);
    REQUIRE(Transforms::greedy_pauli_optimisation().apply(circ));
    REQUIRE(test_unitary_comparison(circ, d, true));
  }
}
SCENARIO("Test GreedyPauliSimp pass construction") {
  // test pass construction
  GIVEN("A circuit") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Rz, 0.5, {1});
    CompilationUnit cu(c);
    CHECK(gen_greedy_pauli_simp(0.3, 0.5)->apply(cu));
    REQUIRE(test_unitary_comparison(c, cu.get_circ_ref(), true));
  }
}
}  // namespace test_GreedyPauliSimp
}  // namespace tket
