// Copyright 2019-2021 Cambridge Quantum Computing
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

#include <catch2/catch.hpp>
#include <numeric>

#include "Circuit/Boxes.hpp"
#include "CircuitsForTesting.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_Rebase {

SCENARIO("Building rebases with rebase_factory") {
  GIVEN("A circuit with all gates in basis set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_op<unsigned>(OpType::Rx, 0.4, {0});
    Circuit copy = c;
    Circuit blank(2);
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::S, OpType::V, OpType::Rx};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, blanker);
    REQUIRE(!t.apply(c));
    REQUIRE(copy == c);
  }
  GIVEN("Rebasing a CZ circuit to a CX gate set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::CZ, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::S, OpType::V, OpType::H};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, blanker);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::CZ) == 0);
    REQUIRE(c.count_gates(OpType::CX) == 1);
    const auto s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a CX circuit to a CZ gate set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit cx(2);
    cx.add_op<unsigned>(OpType::H, {0});
    cx.add_op<unsigned>(OpType::CZ, {0, 1});
    cx.add_op<unsigned>(OpType::H, {0});
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CRz};
    OpTypeSet singleqs = {OpType::S, OpType::V, OpType::H};
    Transform t = Transform::rebase_factory(multiqs, cx, singleqs, blanker);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::CZ) == 1);
    REQUIRE(c.count_gates(OpType::CX) == 0);
    const StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a CY circuit to a CZ gate set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::X, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit cx(2);
    cx.add_op<unsigned>(OpType::H, {0});
    cx.add_op<unsigned>(OpType::CZ, {0, 1});
    cx.add_op<unsigned>(OpType::H, {0});
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CZ};
    OpTypeSet singleqs = {OpType::S, OpType::X, OpType::H, OpType::Sdg};
    Transform t = Transform::rebase_factory(multiqs, cx, singleqs, blanker);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::CZ) == 1);
    REQUIRE(c.count_gates(OpType::CX) == 0);
    REQUIRE(c.count_gates(OpType::CY) == 0);
    const StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a Controlled rotation circuit to a CX gate set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::CRx, 0.5, {0, 1});
    c.add_op<unsigned>(OpType::CRy, 0.5, {0, 1});
    c.add_op<unsigned>(OpType::CRz, 0.5, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::S, OpType::V, OpType::H};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, blanker);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::CZ) == 0);
    REQUIRE(c.count_gates(OpType::CX) == 6);
    const StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a CV and CVdg ciruit to a CX gate set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::CV, {0, 1});
    c.add_op<unsigned>(OpType::CVdg, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::S, OpType::V, OpType::H};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, blanker);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::CX) == 4);
    const StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a CSX and CSXdg ciruit to a CX gate set") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::CSX, {0, 1});
    c.add_op<unsigned>(OpType::CSXdg, {0, 1});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::SXdg, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::S, OpType::V, OpType::H};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, blanker);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::CX) == 4);
    StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a Rx/T sequence to TK1") {
    Circuit c(1);
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::Rx, 0.34, {0});
    c.add_op<unsigned>(OpType::T, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto tk1_map = [](const Expr& theta, const Expr& phi, const Expr& lambda) {
      Circuit u(1);
      std::vector<Expr> params = {theta, phi, lambda};
      u.add_op<unsigned>(OpType::TK1, params, {0});
      return u;
    };
    auto blanker = [](const Expr&, const Expr&, const Expr&) {
      return Circuit(1);
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::TK1};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, tk1_map);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::T) == 0);
    REQUIRE(c.count_gates(OpType::Rx) == 0);
    REQUIRE(c.count_gates(OpType::TK1) == 3);
    StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a U3 sequence to Rz/Rx") {
    Circuit c(1);
    std::vector<Expr> params0 = {0.19, 1.23, 0.58};
    c.add_op<unsigned>(OpType::U3, params0, {0});
    std::vector<Expr> params1 = {1.76, 1.05, 0.24};
    c.add_op<unsigned>(OpType::U3, params1, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto rzrx_map = [](const Expr& alpha, const Expr& beta, const Expr& gamma) {
      Circuit u(1);
      u.add_op<unsigned>(OpType::Rz, gamma, {0});
      u.add_op<unsigned>(OpType::Rx, beta, {0});
      u.add_op<unsigned>(OpType::Rz, alpha, {0});
      return u;
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::Rz, OpType::Rx};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, rzrx_map);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::U3) == 0);
    REQUIRE(c.count_gates(OpType::Rx) == 2);
    REQUIRE(c.count_gates(OpType::Rz) == 4);
    StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("Rebasing a Rx/T sequence to Rz/Rx") {
    Circuit c(1);
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::Rx, 0.34, {0});
    c.add_op<unsigned>(OpType::T, {0});
    const auto s0 = tket_sim::get_statevector(c);
    Circuit blank(2);
    auto rzrx_map = [](const Expr& alpha, const Expr& beta, const Expr& gamma) {
      Circuit u(1);
      u.add_op<unsigned>(OpType::Rz, gamma, {0});
      u.add_op<unsigned>(OpType::Rx, beta, {0});
      u.add_op<unsigned>(OpType::Rz, alpha, {0});
      Transform::remove_redundancies().apply(u);
      return u;
    };
    OpTypeSet multiqs = {OpType::CX};
    OpTypeSet singleqs = {OpType::Rz, OpType::Rx};
    Transform t = Transform::rebase_factory(multiqs, blank, singleqs, rzrx_map);
    REQUIRE(t.apply(c));
    REQUIRE(c.count_gates(OpType::T) == 0);
    REQUIRE(c.count_gates(OpType::U3) == 0);
    REQUIRE(c.count_gates(OpType::Rx) == 1);
    REQUIRE(c.count_gates(OpType::Rz) == 2);
    StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("A circuit to be rebased to ProjectQ gateset") {
    Circuit c(2);
    std::vector<Expr> params = {0.5, 0., 1.};
    c.add_op<unsigned>(OpType::U3, params, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    params = {1., 0., 1.};
    c.add_op<unsigned>(OpType::U3, params, {1});
    const auto s0 = tket_sim::get_statevector(c);
    REQUIRE(Transform::rebase_projectq().apply(c));
    REQUIRE(c.count_gates(OpType::U3) == 0);
    StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("A circuit to be rebased to OQC gateset") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::Rx, 0.34, {0});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::Ry, 0.34, {0});
    c.add_op<unsigned>(OpType::H, {0});
    const auto s0 = tket_sim::get_statevector(c);
    REQUIRE(Transform::rebase_OQC().apply(c));
    REQUIRE(
        c.count_gates(OpType::ECR) + c.count_gates(OpType::Rz) +
            c.count_gates(OpType::SX) ==
        c.n_gates());
    StateVector s1 = tket_sim::get_statevector(c);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("A UCCSD example") {
    auto circ = CircuitsForTesting::get().uccsd;
    const StateVector s0 = tket_sim::get_statevector(circ);
    Transform::rebase_tket().apply(circ);
    REQUIRE(circ.count_gates(OpType::Rz) == 0);
    REQUIRE(circ.count_gates(OpType::Rx) == 0);
    const StateVector s1 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
    Transform::decompose_ZX().apply(circ);
    REQUIRE(circ.count_gates(OpType::TK1) == 0);
    const StateVector s2 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s2));
    Transform::decompose_cliffords_std().apply(circ);
    REQUIRE(circ.count_gates(OpType::Rz) == 2);
    REQUIRE(circ.count_gates(OpType::Rx) == 0);
    REQUIRE(circ.count_gates(OpType::TK1) == 0);
    const StateVector s3 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s3));
  }
  GIVEN("A circuit with conditional gates") {
    Circuit circ(2, 1);
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_conditional_gate<unsigned>(OpType::H, {}, {1}, {0}, 1);
    Transform::rebase_tket().apply(circ);
    Circuit correct(2, 1);
    correct.add_op<unsigned>(OpType::TK1, {0, 0, 0.25}, {0});
    correct.add_conditional_gate<unsigned>(
        OpType::TK1, {0.5, 0.5, 0.5}, {1}, {0}, 1);
    correct.add_phase(0.625);
    REQUIRE(circ == correct);
  }
}

SCENARIO("Decompose all boxes") {
  GIVEN("A quantum-only CircBox") {
    Circuit u(2);
    u.add_op<unsigned>(OpType::Ry, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    CircBox ubox(u);
    Circuit v(2);
    v.add_box(ubox, {0, 1});
    bool success = Transform::decomp_boxes().apply(v);
    REQUIRE(success);
    REQUIRE(u == v);
  }
  GIVEN("A mixed CircBox") {
    Circuit u(2, 1);
    u.add_op<unsigned>(OpType::Ry, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    u.add_measure(0, 0);
    CircBox ubox(u);
    Circuit v(2, 1);
    v.add_box(ubox, {/*qubits*/ 0, 1, /*bits*/ 0});
    bool success = Transform::decomp_boxes().apply(v);
    REQUIRE(success);
    REQUIRE(u == v);
  }
  GIVEN("A custom gate") {
    Circuit u(2);
    Sym a = SymEngine::symbol("a");
    Expr a_expr(a);
    u.add_op<unsigned>(OpType::Ry, a_expr - 0.3, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(u.is_symbolic());
    composite_def_ptr_t def = CompositeGateDef::define_gate("g", u, {a});
    Circuit v(2);
    v.add_box(CustomGate(def, {0.5}), {0, 1});
    REQUIRE(!v.is_symbolic());
    symbol_map_t smap = {{a, 0.5}};
    u.symbol_substitution(smap);
    REQUIRE(!u.is_symbolic());
    bool success = Transform::decomp_boxes().apply(v);
    REQUIRE(success);
    REQUIRE(u == v);
  }
  GIVEN("A conditional box") {
    Circuit u(2, 1);
    u.add_op<unsigned>(OpType::Ry, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    u.add_measure(0, 0);
    CircBox ubox(u);
    Conditional cond(std::make_shared<CircBox>(ubox), 2, 1);
    CHECK(cond.n_qubits() == 2);
    Circuit v(2, 3);
    v.add_op<UnitID>(
        std::make_shared<Conditional>(cond),
        {Bit(0), Bit(1), Qubit(0), Qubit(1), Bit(2)});
    bool success = Transform::decomp_boxes().apply(v);
    REQUIRE(success);
    Circuit compare(2, 3);
    compare.add_conditional_gate<unsigned>(OpType::Ry, {-0.75}, {0}, {0, 1}, 1);
    compare.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 1);
    compare.add_conditional_gate<unsigned>(
        OpType::Measure, {}, {0, 2}, {0, 1}, 1);
    REQUIRE(v == compare);
  }
  GIVEN("A conditional box using an existing bit") {
    Circuit u(2, 1);
    u.add_op<unsigned>(OpType::Ry, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    u.add_conditional_gate<unsigned>(OpType::X, {}, {1}, {0}, 1);
    CircBox ubox(u);
    Conditional cond(std::make_shared<CircBox>(ubox), 2, 1);
    Circuit v(2, 2);
    v.add_op<UnitID>(
        std::make_shared<Conditional>(cond),
        {Bit(0), Bit(1), Qubit(0), Qubit(1), Bit(0)});
    bool success = Transform::decomp_boxes().apply(v);
    REQUIRE(success);
    REQUIRE_NOTHROW(v.get_commands());
  }
}

SCENARIO("Check each Clifford case for tk1_to_rzh") {
  GIVEN("Each case") {
    struct RzHTestCase {
      Expr alpha;
      Expr beta;
      Expr gamma;
      unsigned expected_gates;
    };
    const std::vector<RzHTestCase> cases = {
        {0.234, 0., 0.953, 1},    {0.234, 0.5, 0.953, 3},
        {0.234, 1., 0.953, 4},    {0.234, 1.5, 0.953, 3},
        {0.234, 2., 0.953, 1},    {0.234, 2.5, 0.953, 3},
        {0.234, 3., 0.953, 4},    {0.234, 3.5, 0.953, 3},
        {0.234, 0.354, 0.953, 5},
    };
    for (const RzHTestCase& test : cases) {
      Circuit correct(1);
      correct.add_op<unsigned>(
          OpType::TK1, {test.alpha, test.beta, test.gamma}, {0});
      Circuit result = Transform::tk1_to_rzh(test.alpha, test.beta, test.gamma);
      REQUIRE(result.n_gates() == test.expected_gates);
      REQUIRE(test_unitary_comparison(correct, result));
    }
  }
}

SCENARIO("Check cases for tk1_to_rzsx") {
  GIVEN("Each case") {
    struct RzSXTestCase {
      Expr alpha;
      Expr beta;
      Expr gamma;
      unsigned expected_gates;
    };
    const std::vector<RzSXTestCase> cases = {
        {1.5, 0., 2.53, 1},     {1.5, 2., 2.53, 1},       {1.5, 4., 2.53, 1},

        {1.5, 0., 2.5, 0},      {1.5, 2., 2.5, 0},        {1.5, 4., 2.5, 0},

        {2., 1., 0., 2},        {2, -1., 0., 2},          {2., 3., 0., 2},
        {4., 1., 0., 2},        {4., -1., 0., 2},         {4., 3., 0., 2},
        {2., 1., 8., 2},        {2., -1., 8., 2},         {2., 3., 8., 2},

        {3., 1., 0.5, 4},       {3., -1., 0.5, 4},        {7., 3., 0.5, 4},
        {3.5, 1., 0.5, 4},      {0.3, -1., 3.7, 4},       {0.3, 3., 1, 4},

        {0.5, 1.3, 0.5, 3},     {2.5, 1.3, 0.5, 3},       {0.5, 1.3, 2.5, 3},
        {2.5, 1.3, 2.5, 3},     {-1.5, 1.3, -1.5, 3},

        {0., 0.5, 0., 1},       {0., 0.5, 2., 1},         {0., 2.5, 0., 1},
        {0., 2.5, 2., 1},       {2., 0.5, 0., 1},         {2., 0.5, 2., 1},
        {2., 2.5, 0., 1},       {2., 2.5, 2., 1},

        {7.55, 1.3, 1.55, 5},   {3.55, 1.3, 0.5, 5},      {8.53, 1.3, 2.25, 5},
        {0., -0.5, 0, 5},       {0., 1.5, 0, 5},          {0., 3.5, 0, 5},
        {2., -0.5, 0, 5},       {2., 1.5, 0, 5},          {2., 3.5, 0, 5},
        {0.234, 3.5, 0.953, 5}, {0.234, 0.354, 0.953, 5},
    };
    const Expr a_expr(SymEngine::symbol("a"));
    const Expr b_expr(SymEngine::symbol("b"));
    const Expr c_expr(SymEngine::symbol("c"));
    const std::vector<RzSXTestCase> symbolic_cases = {
        {a_expr, 2, c_expr, 1},
        {a_expr, 1, c_expr, 4},
        {a_expr, 1, a_expr, 2},
        {2.5, b_expr, 0.5, 3},
        {a_expr, b_expr, c_expr, 5}};

    for (const RzSXTestCase& test : cases) {
      Circuit correct(1);
      correct.add_op<unsigned>(
          OpType::TK1, {test.alpha, test.beta, test.gamma}, {0});
      Circuit result =
          Transform::tk1_to_rzsx(test.alpha, test.beta, test.gamma);
      REQUIRE(result.n_gates() == test.expected_gates);
      REQUIRE(test_unitary_comparison(correct, result));
    }

    for (const RzSXTestCase& test : symbolic_cases) {
      Circuit correct(1);
      correct.add_op<unsigned>(
          OpType::TK1, {test.alpha, test.beta, test.gamma}, {0});
      Circuit result =
          Transform::tk1_to_rzsx(test.alpha, test.beta, test.gamma);
      REQUIRE(result.n_gates() == test.expected_gates);
    }
  }
}

}  // namespace test_Rebase
}  // namespace tket
