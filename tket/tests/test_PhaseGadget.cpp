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
#include <utility>

#include "Circuit/CircUtils.hpp"
#include "CircuitsForTesting.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/CliffordOptimisation.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/PhaseOptimisation.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/EigenConfig.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_PhaseGadget {

SCENARIO("Convert into PhaseGadgets", "[transform]") {
  GIVEN("A circuit which should be converted") {
    Circuit circ(2);
    Vertex v1 = circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 1e-4, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    bool success = Transforms::decompose_PhaseGadgets().apply(circ);
    REQUIRE(success);
    REQUIRE(circ.n_vertices() == 5);
    REQUIRE(circ.get_OpType_from_Vertex(v1) == OpType::PhaseGadget);
  }
  GIVEN("A circuit which shouldn't be converted") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 1e-4, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    bool success = Transforms::decompose_PhaseGadgets().apply(circ);
    REQUIRE(!success);
  }
}

SCENARIO("Smash CXs using PhaseGadgets", "[transform]") {
  GIVEN("A circuit which can be converted then smashed") {
    WHEN("The circuit is converted into PhaseGadgets") {
      Circuit circ(3);
      circ.add_op<unsigned>(OpType::CX, {2, 0});
      Vertex v1 = circ.add_op<unsigned>(OpType::CX, {0, 1});
      circ.add_op<unsigned>(OpType::Rz, 1e-4, {1});
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      circ.add_op<unsigned>(OpType::CX, {2, 0});
      bool success1 = Transforms::decompose_PhaseGadgets().apply(circ);
      REQUIRE(success1);
      THEN("The circuit is smashed") {
        bool success2 = Transforms::smash_CX_PhaseGadgets().apply(circ);
        REQUIRE(success2);
        REQUIRE(test_equiv_val(
            (circ.get_Op_ptr_from_Vertex(v1))->get_params()[0], 1e-4));
        REQUIRE(circ.n_in_edges(v1) == 3);
      }
    }
  }
  GIVEN("A bigger circuit which can be smashed") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{1, 0}, {2, 0}, {4, 0}, {3, 0}});
    circ.add_op<unsigned>(OpType::Rz, 1e-3, {0});
    add_2qb_gates(circ, OpType::CX, {{3, 0}, {4, 0}, {2, 0}, {1, 0}});
    REQUIRE(verify_n_qubits_for_ops(circ));
    REQUIRE(Transforms::decompose_PhaseGadgets().apply(circ));
    REQUIRE(verify_n_qubits_for_ops(circ));
    REQUIRE(Transforms::smash_CX_PhaseGadgets().apply(circ));
    REQUIRE(verify_n_qubits_for_ops(circ));
    REQUIRE(circ.n_vertices() == 11);
    REQUIRE(circ.n_edges() == 10);
  }
  GIVEN("A circuit which cannot be smashed") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{1, 0}, {2, 0}, {4, 0}, {3, 0}});
    circ.add_op<unsigned>(OpType::Rz, 1e-3, {0});
    add_2qb_gates(circ, OpType::CX, {{3, 0}, {1, 0}, {2, 0}, {4, 0}});
    REQUIRE(verify_n_qubits_for_ops(circ));
    REQUIRE(Transforms::decompose_PhaseGadgets().apply(circ));
    REQUIRE(verify_n_qubits_for_ops(circ));
    REQUIRE(!Transforms::smash_CX_PhaseGadgets().apply(circ));
  }
}

SCENARIO("Aligning ports on PhaseGadgets", "[transform]") {
  GIVEN("A sequence of PhaseGadgets, trying to align all ports") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::PhaseGadget, 0.5, {0, 1, 2, 3});
    add_1qb_gates(circ, OpType::X, {1, 2});
    circ.add_op<unsigned>(OpType::PhaseGadget, 0.25, {3, 2, 1, 0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::PhaseGadget, 0.75, {0, 1, 2});
    Transforms::align_PhaseGadgets().apply(circ);
    VertexVec vertices = circ.vertices_in_order();
    Vertex first_gadget = vertices[4];
    bool matching_ports = true;
    for (int i = 0; i < 4; i++) {
      Edge e_i = circ.get_nth_out_edge(first_gadget, i);
      Vertex next = boost::target(e_i, circ.dag);
      if (circ.get_OpType_from_Vertex(next) == OpType::PhaseGadget) {
        matching_ports &=
            (circ.dag[e_i].ports.first == circ.dag[e_i].ports.second);
      } else {
        Edge to_next_gadget = circ.get_next_edge(next, e_i);
        matching_ports &=
            (circ.dag[e_i].ports.first ==
             circ.dag[to_next_gadget].ports.second);
      }
    }
    REQUIRE(matching_ports);
  }
}

SCENARIO("Full optimise_via_PhaseGadget") {
  GIVEN("A UCCSD example") {
    auto circ = CircuitsForTesting::get().uccsd;
    const auto s0 = tket_sim::get_unitary(circ);
    Transforms::optimise_via_PhaseGadget(CXConfigType::Tree).apply(circ);
    REQUIRE(circ.count_gates(OpType::CX) == 12);
    REQUIRE(circ.depth() == 13);
    const auto s1 = tket_sim::get_unitary(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
}

SCENARIO("Constructing Pauli gadgets") {
  GIVEN("XY") {
    double t = 0.3;
    Eigen::Matrix4cd a;  // kronecker product of X and Y
    // clang-format off
        a << 0., 0., 0., -1.*i_,
             0., 0., 1.*i_, 0.,
             0., -1.*i_, 0., 0.,
             1.*i_, 0., 0., 0.;
    // clang-format on
    Eigen::Matrix4cd m = (+0.5 * PI * i_ * t * a).exp();
    Circuit circ = pauli_gadget({Pauli::X, Pauli::Y}, t);
    const Eigen::Matrix4cd u = tket_sim::get_unitary(circ);
    Eigen::Matrix4cd v = m * u;
    REQUIRE((v - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
}

SCENARIO("Identifying and synthesising Pauli gadgets") {
  const auto get_test_circ = []() {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    add_1qb_gates(circ, OpType::V, {2, 3});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {2, 3}});
    circ.add_op<unsigned>(OpType::Rz, 0.6, {3});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {1, 2}, {0, 1}});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    add_1qb_gates(circ, OpType::Vdg, {2, 3});
    return circ;
  };
  GIVEN("A single Pauli gadget") {
    auto circ = get_test_circ();
    Transforms::pairwise_pauli_gadgets().apply(circ);
    Transforms::singleq_clifford_sweep().apply(circ);
    Transforms::synthesise_tket().apply(circ);
    Circuit expected(4);
    expected.add_op<unsigned>(OpType::V, {0});
    expected.add_op<unsigned>(OpType::H, {1});
    add_1qb_gates(expected, OpType::V, {2, 3});
    expected.add_op<unsigned>(OpType::PhaseGadget, 0.6, {0, 1, 2, 3});
    expected.add_op<unsigned>(OpType::Vdg, {0});
    expected.add_op<unsigned>(OpType::H, {1});
    add_1qb_gates(expected, OpType::Vdg, {2, 3});
    Transforms::decompose_multi_qubits_CX().apply(expected);
    Transforms::singleq_clifford_sweep().apply(expected);
    Transforms::synthesise_tket().apply(expected);
    const auto m1 = tket_sim::get_unitary(circ);
    const auto m2 = tket_sim::get_unitary(expected);
    Eigen::MatrixXcd m = m1 * m2.conjugate().transpose();
    REQUIRE(
        (m - Eigen::MatrixXcd::Identity(16, 16)).cwiseAbs().sum() < ERR_EPS);

    // TODO:: Handle measurements once we have classical info in Commands
  }
  GIVEN("Pauli-gadget simplification of a symbolic circuit") {
    // --X--Rz(a)--X--  ==>  --U1(-a)--
    Circuit circ(1);
    Sym a = SymEngine::symbol("alpha");
    Expr alpha(a);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Rz, alpha, {0});
    circ.add_op<unsigned>(OpType::X, {0});
    REQUIRE(Transforms::pairwise_pauli_gadgets().apply(circ));
    REQUIRE(circ.n_gates() == 1);
    Vertex v = *circ.get_gates_of_type(OpType::TK1).begin();
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    std::vector<Expr> angles = op->get_params();
    REQUIRE(test_equiv_0(angles[0]));
    REQUIRE(test_equiv_0(angles[1]));
    REQUIRE(test_equiv_0(alpha + angles[2]));
  }
  GIVEN("Another symbolic circuit") {
    // --Rz(a)--Z--Ry(0.5)--Z--  ==>  --U3(-0.5, 0, a)--
    Circuit circ(1);
    Sym a = SymEngine::symbol("alpha");
    Expr alpha(a);
    circ.add_op<unsigned>(OpType::Rz, alpha, {0});
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.5, {0});
    circ.add_op<unsigned>(OpType::Z, {0});
    REQUIRE(Transforms::pairwise_pauli_gadgets().apply(circ));
    REQUIRE(circ.n_gates() == 1);
    Vertex v = *circ.get_gates_of_type(OpType::TK1).begin();
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    std::vector<Expr> angles = op->get_params();
    REQUIRE(test_equiv_0(angles[0] - 0.5));
    REQUIRE(test_equiv_0(angles[1] + 0.5));
    REQUIRE(test_equiv_0(angles[2] + 0.5 - alpha));
  }
  GIVEN("A pair of Pauli gadgets") {
    auto circ = get_test_circ();
    add_1qb_gates(circ, OpType::V, {0, 1, 2});
    circ.add_op<unsigned>(OpType::H, {3});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {2, 3}});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {1, 2}, {0, 1}});
    add_1qb_gates(circ, OpType::Vdg, {0, 1, 2});
    circ.add_op<unsigned>(OpType::H, {3});
    Transforms::pairwise_pauli_gadgets().apply(circ);
    REQUIRE(circ.count_gates(OpType::CX) == 6);
  }
  GIVEN("A sequence of 5 Pauli gadgets") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.6, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rx, 0.2, {1});

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rx, 1.25, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    Circuit copy(circ);
    Transforms::pairwise_pauli_gadgets().apply(circ);
    REQUIRE(test_statevector_comparison(circ, copy));
    REQUIRE(circ.count_gates(OpType::CX) == 2);
  }
  GIVEN("A UCCSD example") {
    auto circ = CircuitsForTesting::get().uccsd;
    const auto s0 = tket_sim::get_statevector(circ);
    Transforms::pairwise_pauli_gadgets().apply(circ);
    REQUIRE(circ.count_gates(OpType::CX) == 6);
    const auto s1 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s0, s1));
  }
}

SCENARIO("Decompose phase gadgets") {
  auto get_circ_decomposition = [](unsigned n_qubits, CXConfigType config) {
    std::vector<Pauli> nZ;
    for (unsigned i = 0; i < n_qubits; ++i) {
      nZ.push_back(Pauli::Z);
    }
    Circuit circ = pauli_gadget(nZ, 0.2);
    return std::make_pair(circ, phase_gadget(n_qubits, 0.2, config));
  };
  GIVEN("Using Star strategy") {
    CXConfigType config = CXConfigType::Star;
    auto [circ, decomp] = get_circ_decomposition(4, config);
    THEN("decomposition has 6x CX") {
      REQUIRE(decomp.count_gates(OpType::CX) == 6);
    }
    THEN("unitary is correct") {
      auto u1 = tket_sim::get_unitary(circ);
      auto u2 = tket_sim::get_unitary(decomp);
      REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    }
  }
  GIVEN("Using Tree strategy") {
    CXConfigType config = CXConfigType::Tree;
    auto [circ, decomp] = get_circ_decomposition(4, config);
    THEN("decomposition has 6x CX") {
      REQUIRE(decomp.count_gates(OpType::CX) == 6);
    }
    THEN("unitary is correct") {
      auto u1 = tket_sim::get_unitary(circ);
      auto u2 = tket_sim::get_unitary(decomp);
      REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    }
  }
  GIVEN("Using Snake strategy") {
    CXConfigType config = CXConfigType::Snake;
    auto [circ, decomp] = get_circ_decomposition(4, config);
    THEN("decomposition has 6x CX") {
      REQUIRE(decomp.count_gates(OpType::CX) == 6);
    }
    THEN("unitary is correct") {
      auto u1 = tket_sim::get_unitary(circ);
      auto u2 = tket_sim::get_unitary(decomp);
      REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    }
  }
  GIVEN("Using MultiQGate strategy, 4 qubits") {
    CXConfigType config = CXConfigType::MultiQGate;
    auto [circ, decomp] = get_circ_decomposition(4, config);
    THEN("decomposition has 2x XXPhase3") {
      REQUIRE(decomp.count_gates(OpType::XXPhase3) == 2);
    }
    THEN("decomposition has 2x CX") {
      REQUIRE(decomp.count_gates(OpType::CX) == 2);
    }
    THEN("unitary is correct") {
      auto u1 = tket_sim::get_unitary(circ);
      auto u2 = tket_sim::get_unitary(decomp);
      REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    }
  }
  GIVEN("Using MultiQGate strategy, 5 qubits") {
    CXConfigType config = CXConfigType::MultiQGate;
    auto [circ, decomp] = get_circ_decomposition(5, config);
    THEN("decomposition has 4x XXPhase3") {
      REQUIRE(decomp.count_gates(OpType::XXPhase3) == 4);
    }
    THEN("unitary is correct") {
      auto u1 = tket_sim::get_unitary(circ);
      auto u2 = tket_sim::get_unitary(decomp);
      REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    }
  }
}

}  // namespace test_PhaseGadget
}  // namespace tket
