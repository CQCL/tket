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
#include <set>

#include "Circuit/CircUtils.hpp"
#include "CircuitsForTesting.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/CliffordOptimisation.hpp"
#include "Transformations/CliffordReductionPass.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/PauliStrings.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_Clifford {

SCENARIO("Test decomposition into Clifford gates", "[transform]") {
  GIVEN(
      "STD FORM: A TK1 instance for each set of parameters fitting "
      "multiples of pi/2") {
    for (int alpha = 0; alpha < 4; alpha++) {
      for (int beta = 0; beta < 4; beta++) {
        for (int gamma = 0; gamma < 4; gamma++) {
          Circuit circ(1);
          std::vector<Expr> params({alpha * 0.5, beta * 0.5, gamma * 0.5});
          circ.add_op<unsigned>(OpType::TK1, params, {0});
          Eigen::Matrix2cd m_before = get_matrix_from_circ(circ);
          REQUIRE(Transforms::decompose_cliffords_std().apply(circ));
          Transforms::decompose_single_qubits_TK1().apply(circ);
          Eigen::Matrix2cd m_after = get_matrix_from_circ(circ);
          REQUIRE(m_before.isApprox(m_after));
        }
      }
    }
  }

  GIVEN("STD FORM: An incompatible circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::U1, 1e-6, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {1});
    bool success = Transforms::decompose_cliffords_std().apply(circ);
    REQUIRE(!success);
  }

  GIVEN("STD FORM: Negative parameters") {
    Circuit circ(1);
    std::vector<Expr> params = {0.5, -0.5, 0.5};
    circ.add_op<unsigned>(OpType::U3, params, {0});
    bool success = Transforms::decompose_cliffords_std().apply(circ);
    REQUIRE(success);
    VertexVec vertices = circ.vertices_in_order();
    REQUIRE(circ.get_OpType_from_Vertex(vertices[1]) == OpType::V);
  }
}

SCENARIO("Check that singleq_clifford_sweep reduces to standard forms") {
  GIVEN("A circuit in the standard form") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::TK1, {0, 0, 0.31}, {1});
    Circuit circ2(circ);
    Transforms::singleq_clifford_sweep().apply(circ);
    REQUIRE(circ2 == circ);
  }
  GIVEN("Some U3s with only pi/2 angles") {
    Circuit circ(2);
    std::vector<Expr> p = {0.5, 1., 0.};
    circ.add_op<unsigned>(OpType::U3, p, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    p = {0.5, 0., 0.5};
    circ.add_op<unsigned>(OpType::U3, p, {0});
    REQUIRE(Transforms::singleq_clifford_sweep().apply(circ));
    Circuit correct(2);
    correct.add_op<unsigned>(OpType::Z, {0});
    correct.add_op<unsigned>(OpType::X, {0});
    correct.add_op<unsigned>(OpType::S, {0});
    correct.add_op<unsigned>(OpType::V, {0});
    correct.add_op<unsigned>(OpType::S, {0});
    correct.add_op<unsigned>(OpType::CX, {0, 1});
    correct.add_op<unsigned>(OpType::V, {0});
    correct.add_op<unsigned>(OpType::S, {0});
    REQUIRE(circ == correct);
  }
  GIVEN("Some Zs, Xs, and rotations on qubit 0 to commute/copy") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    REQUIRE(Transforms::singleq_clifford_sweep().apply(circ));
    Circuit correct(2);
    correct.add_op<unsigned>(OpType::Z, {0});
    correct.add_op<unsigned>(OpType::X, {0});
    correct.add_op<unsigned>(OpType::X, {1});
    correct.add_op<unsigned>(OpType::S, {0});
    correct.add_op<unsigned>(OpType::CX, {0, 1});
    correct.add_op<unsigned>(OpType::V, {0});
    REQUIRE(circ == correct);
  }
  GIVEN("Some Zs, Xs, and rotations on qubit 1 to commute/copy") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    REQUIRE(Transforms::singleq_clifford_sweep().apply(circ));
    Circuit correct(2);
    correct.add_op<unsigned>(OpType::Z, {0});
    correct.add_op<unsigned>(OpType::Z, {1});
    correct.add_op<unsigned>(OpType::X, {1});
    correct.add_op<unsigned>(OpType::V, {1});
    correct.add_op<unsigned>(OpType::CX, {0, 1});
    correct.add_op<unsigned>(OpType::S, {1});
    REQUIRE(circ == correct);
  }
  GIVEN("Mixtures of copying and commuting") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::S, {0});
    REQUIRE(Transforms::singleq_clifford_sweep().apply(circ));
    Circuit correct(2);
    correct.add_op<unsigned>(OpType::S, {0});
    correct.add_op<unsigned>(OpType::Z, {1});
    correct.add_op<unsigned>(OpType::X, {1});
    correct.add_op<unsigned>(OpType::V, {1});
    correct.add_op<unsigned>(OpType::CX, {0, 1});
    correct.add_op<unsigned>(OpType::S, {1});
    REQUIRE(circ == correct);
  }
}

SCENARIO("Rewriting Clifford subcircuits", "[transform]") {
  Circuit circ(2);
  for (int nn = 0; nn < 4; ++nn) {
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
  }
  REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
}

SCENARIO("valid_insertion_point returns space-like interaction points") {
  Circuit circ(4);
  Vertex cx1 = circ.add_op<unsigned>(OpType::CX, {0, 1});
  Vertex cx2 = circ.add_op<unsigned>(OpType::CX, {2, 3});
  Vertex cx3 = circ.add_op<unsigned>(OpType::CX, {1, 2});
  Vertex cx4 = circ.add_op<unsigned>(OpType::CX, {1, 2});
  CliffordReductionPassTester clifford_pass(circ);

  auto e1 = circ.get_nth_out_edge(cx1, 1);       // edge cx1 - cx3
  auto e2 = circ.get_nth_out_edge(cx2, 0);       // edge cx2 - cx3
  auto e_final = circ.get_nth_out_edge(cx3, 0);  // edge  cx3 - cx4

  auto to_ip = [cx1](Edge e) {
    // we only care about the edge e
    return InteractionPoint{e, cx1, Pauli::I, false};
  };
  std::list<InteractionPoint> seq0{to_ip(e1), to_ip(e_final)};
  std::list<InteractionPoint> seq1{to_ip(e2)};

  auto ips = clifford_pass.valid_insertion_point(seq0, seq1);

  REQUIRE(ips);
  REQUIRE((ips->first).e == e1);
  REQUIRE((ips->second).e == e2);
}

SCENARIO("ham3tc.qasm file was breaking for canonical clifford transform") {
  Circuit circ(5);
  circ.add_op<unsigned>(OpType::H, {1});
  circ.add_op<unsigned>(OpType::H, {4});
  circ.add_op<unsigned>(OpType::T, {1});
  circ.add_op<unsigned>(OpType::S, {1});
  circ.add_op<unsigned>(OpType::Z, {1});
  circ.add_op<unsigned>(OpType::CX, {2, 1});
  circ.add_op<unsigned>(OpType::CX, {3, 1});
  circ.add_op<unsigned>(OpType::T, {1});
  circ.add_op<unsigned>(OpType::CX, {2, 1});
  circ.add_op<unsigned>(OpType::T, {2});
  circ.add_op<unsigned>(OpType::T, {1});
  circ.add_op<unsigned>(OpType::S, {1});
  circ.add_op<unsigned>(OpType::CX, {3, 1});
  circ.add_op<unsigned>(OpType::CX, {3, 2});
  circ.add_op<unsigned>(OpType::T, {2});
  circ.add_op<unsigned>(OpType::S, {2});
  circ.add_op<unsigned>(OpType::CX, {3, 2});
  circ.add_op<unsigned>(OpType::CX, {2, 1});
  circ.add_op<unsigned>(OpType::T, {1});
  circ.add_op<unsigned>(OpType::S, {1});
  circ.add_op<unsigned>(OpType::CX, {2, 1});
  circ.add_op<unsigned>(OpType::T, {3});
  circ.add_op<unsigned>(OpType::S, {4});
  circ.add_op<unsigned>(OpType::Z, {4});
  circ.add_op<unsigned>(OpType::CX, {4, 1});
  circ.add_op<unsigned>(OpType::S, {1});
  add_2qb_gates(
      circ, OpType::CX,
      {{4, 1}, {1, 4}, {4, 1}, {1, 4}, {3, 2}, {2, 3}, {1, 3}, {3, 2}});
  circ.add_op<unsigned>(OpType::H, {4});
  circ.add_op<unsigned>(OpType::Collapse, {4});
  WHEN("Hyper Clifford Squash") {
    REQUIRE(Transforms::canonical_hyper_clifford_squash().apply(circ));
  }
  WHEN("Clifford Simp") { REQUIRE(Transforms::clifford_simp().apply(circ)); }
}

SCENARIO("Test multiq clifford replacements") {
  GIVEN("Replacement number 1") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    Circuit replacement1(2);
    replacement1.add_op<unsigned>(OpType::Z, {1});
    replacement1.add_op<unsigned>(OpType::S, {0});
    replacement1.add_op<unsigned>(OpType::S, {1});
    replacement1.add_op<unsigned>(OpType::CX, {0, 1});
    replacement1.add_op<unsigned>(OpType::V, {0});
    replacement1.add_op<unsigned>(OpType::S, {0});
    replacement1.add_op<unsigned>(OpType::S, {1});
    REQUIRE(replacement1 == circ);
  }
  GIVEN("Replacement number 2") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    Circuit replacement2(2);
    replacement2.add_op<unsigned>(OpType::X, {1});
    replacement2.add_op<unsigned>(OpType::V, {0});
    replacement2.add_op<unsigned>(OpType::V, {1});
    replacement2.add_op<unsigned>(OpType::CX, {1, 0});
    replacement2.add_op<unsigned>(OpType::S, {0});
    replacement2.add_op<unsigned>(OpType::V, {0});
    replacement2.add_op<unsigned>(OpType::V, {1});
    replacement2.add_phase(0.75);
    REQUIRE(circ == replacement2);
  }
  GIVEN("Replacement number 3") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::S) == 2);
    REQUIRE(circ.count_gates(OpType::V) == 2);
    REQUIRE(circ.n_vertices() == 8);
  }
  GIVEN("Replacement number 5") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    Circuit replacement5(2);
    replacement5.add_op<unsigned>(OpType::S, {0});
    replacement5.add_op<unsigned>(OpType::Z, {1});
    replacement5.add_op<unsigned>(OpType::S, {1});
    replacement5.add_op<unsigned>(OpType::V, {1});
    replacement5.add_op<unsigned>(OpType::S, {1});
    replacement5.add_op<unsigned>(OpType::CX, {0, 1});
    replacement5.add_op<unsigned>(OpType::S, {1});
    replacement5.add_op<unsigned>(OpType::V, {1});
    REQUIRE(circ == replacement5);
  }
  GIVEN("Replacement number 6") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::U1, 0.2, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::U1, 0.4, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::V) == 4);
    REQUIRE(circ.count_gates(OpType::X) == 1);
    REQUIRE(circ.count_gates(OpType::S) == 2);
    REQUIRE(circ.count_gates(OpType::U1) == 2);
  }
  GIVEN("Replacement number 7") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});

    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));

    Circuit comp(2);
    comp.add_op<unsigned>(OpType::V, {0});
    comp.add_op<unsigned>(OpType::S, {0});
    comp.add_op<unsigned>(OpType::CX, {1, 0});
    comp.add_op<unsigned>(OpType::S, {0});
    comp.add_op<unsigned>(OpType::S, {0});
    comp.add_op<unsigned>(OpType::V, {0});
    comp.add_op<unsigned>(OpType::X, {0});
    comp.add_op<unsigned>(OpType::V, {0});
    comp.add_op<unsigned>(OpType::S, {1});
    comp.add_op<unsigned>(OpType::V, {1});
    comp.add_op<unsigned>(OpType::S, {1});
    comp.add_op<unsigned>(OpType::Z, {1});
    comp.add_op<unsigned>(OpType::V, {1});
    comp.add_op<unsigned>(OpType::S, {1});
    comp.add_op<unsigned>(OpType::V, {1});
    comp.add_op<unsigned>(OpType::CX, {1, 0});
    comp.add_op<unsigned>(OpType::V, {0});
    comp.add_op<unsigned>(OpType::Z, {0});
    comp.add_op<unsigned>(OpType::S, {0});
    comp.add_op<unsigned>(OpType::S, {1});
    comp.add_op<unsigned>(OpType::V, {1});
    comp.add_op<unsigned>(OpType::S, {1});
    comp.add_op<unsigned>(OpType::CX, {0, 1});
    comp.add_op<unsigned>(OpType::V, {0});
    comp.add_op<unsigned>(OpType::S, {0});
    comp.add_op<unsigned>(OpType::S, {1});
    REQUIRE(circ.count_gates(OpType::Z) == comp.count_gates(OpType::Z));
    REQUIRE(circ.count_gates(OpType::S) == comp.count_gates(OpType::S));
    REQUIRE(circ.count_gates(OpType::CX) == comp.count_gates(OpType::CX));
    REQUIRE(circ.count_gates(OpType::X) == comp.count_gates(OpType::X));
  }
  GIVEN("Test that replacements will not break causal ordering") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {2, 0});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    REQUIRE(!Transforms::multiq_clifford_replacement(true).apply(circ));
    REQUIRE_NOTHROW(circ.depth_by_type(OpType::CX));
  }
}

SCENARIO("Test clifford reduction") {
  GIVEN("Replacement number 1") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Replacement number 2") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Replacement number 3") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 0);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Replacement number 5") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Replacement number 6") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::U1, 0.2, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::U1, 0.4, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 0);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Replacement number 7") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Test that replacements will not break causal ordering") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    add_2qb_gates(circ, OpType::CX, {{2, 0}, {3, 2}, {2, 1}, {1, 0}});
    REQUIRE_FALSE(Transforms::clifford_reduction(true).apply(circ));
    REQUIRE_NOTHROW(circ.depth_by_type(OpType::CX));
  }
  GIVEN("Circuit with a selection of Clifford gates") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CY, {1, 0});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::ZZMax, {0, 1});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction(true).apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::CY) == 0);
    REQUIRE(circ.count_gates(OpType::CZ) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Circuit with non-Clifford gates") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {2, 1});
    circ.add_op<unsigned>(OpType::ZZPhase, 0.3, {0, 1});
    circ.add_op<unsigned>(OpType::Rx, 0.1, {2});
    circ.add_op<unsigned>(OpType::CH, {2, 1});
    circ.add_op<unsigned>(OpType::CnRy, 0.2, {1, 0});
    circ.add_op<unsigned>(OpType::CY, {2, 0});
    Circuit copy(circ);
    REQUIRE(Transforms::clifford_reduction(true).apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::CY) == 0);
    REQUIRE(circ.count_gates(OpType::CZ) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    Transforms::rebase_tket().apply(circ);
    Transforms::rebase_tket().apply(copy);
    REQUIRE(test_unitary_comparison(circ, copy));
  }
  GIVEN("Circuit with no possible reductions from this method") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {0, 1}, {1, 2}});
    circ.add_op<unsigned>(OpType::Rx, 0.2, {1});
    add_2qb_gates(circ, OpType::ZZMax, {{0, 1}, {1, 2}, {2, 0}});
    circ.add_op<unsigned>(OpType::Ry, 0.1, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    REQUIRE_FALSE(Transforms::clifford_reduction(true).apply(circ));
  }
}

SCENARIO("Test clifford replacements that allow for SWAPs") {
  GIVEN("allow_swaps 1") {
    Circuit circ(2);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}});
    Circuit original = circ;
    REQUIRE_FALSE(Transforms::clifford_reduction(false).apply(circ));
    REQUIRE(Transforms::clifford_reduction(true).apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(original, circ));
  }
  GIVEN("allow_swaps 2") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit original = circ;
    REQUIRE_FALSE(Transforms::clifford_reduction(false).apply(circ));
    REQUIRE(Transforms::clifford_reduction(true).apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 1);
    REQUIRE(test_unitary_comparison(original, circ));
  }
  GIVEN("Test them both") {
    Circuit circ(4);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {0, 1}});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {2});
    add_2qb_gates(circ, OpType::CX, {{0, 2}, {1, 3}, {3, 1}, {1, 3}});
    Circuit original = circ;
    REQUIRE(Transforms::clifford_reduction(true).apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 1);
    REQUIRE(circ.count_gates(OpType::ZZMax) == 0);
    REQUIRE(test_unitary_comparison(original, circ));
  }
}

SCENARIO("Test Clifford matching plays well with commuting gates") {
  GIVEN("A commuting section at start on first qubit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN("A commuting section at end on first qubit (matching CX direction)") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN("A commuting section at start on second qubit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN(
      "A commuting section at end on second qubit (matching CX "
      "direction)") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN("A commuting section at end on first qubit (opposite CX direction)") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::CX, {2, 0});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN(
      "A commuting section at end on second qubit (opposite CX "
      "direction)") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN("A mixture of all commuting regions") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {0, 2}, {2, 1}});
    circ.add_op<unsigned>(OpType::X, {2});
    add_2qb_gates(circ, OpType::CX, {{2, 0}, {1, 2}, {1, 0}});
    REQUIRE(Transforms::multiq_clifford_replacement(true).apply(circ));

    Circuit correct(3);
    add_2qb_gates(correct, OpType::CX, {{0, 2}, {2, 1}, {1, 0}});
    correct.add_op<unsigned>(OpType::X, {2});
    add_2qb_gates(correct, OpType::CX, {{2, 1}, {0, 2}});
    REQUIRE(
        circ.circuit_equality(correct, {Circuit::Check::ImplicitPermutation}));
  }
  GIVEN("A Hadamard in the non-useful decomposition") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(Transforms::multiq_clifford_replacement().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 3);
  }
}

SCENARIO("Testing full clifford_simp") {
  GIVEN("A UCCSD example") {
    auto circ = CircuitsForTesting::get().uccsd;
    const StateVector s0 = tket_sim::get_statevector(circ);
    Transforms::optimise_via_PhaseGadget(CXConfigType::Tree).apply(circ);
    Transforms::clifford_simp().apply(circ);
    circ.assert_valid();
    REQUIRE(circ.count_gates(OpType::CX) == 8);
    const StateVector s1 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
  GIVEN("A TK2 equivalent to 1CX") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::TK2, {0.5, 0, 0}, {0, 1});
    REQUIRE(Transforms::clifford_simp().apply(circ));
    REQUIRE(circ.count_gates(OpType::CX) == 1);
  }
}

}  // namespace test_Clifford
}  // namespace tket
