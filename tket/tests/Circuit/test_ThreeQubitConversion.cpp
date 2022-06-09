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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "../Simulation/ComparisonFunctions.hpp"
#include "Circuit/Boxes.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/Command.hpp"
#include "Circuit/ThreeQubitConversion.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/ThreeQubitSquash.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Constants.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/UnitID.hpp"

namespace tket {
namespace test_ThreeQubitConversion {

static void check_three_qubit_synthesis(const Eigen::MatrixXcd &U) {
  static const std::set<OpType> expected_1q_gates = {
      OpType::TK1, OpType::H, OpType::Ry, OpType::Rz};
  Circuit c = three_qubit_synthesis(U);
  Transforms::decompose_TK2().apply(c);
  unsigned n_cx = 0;
  for (const Command &cmd : c) {
    OpType optype = cmd.get_op_ptr()->get_type();
    if (optype == OpType::CX) {
      n_cx++;
    } else {
      CHECK(expected_1q_gates.contains(optype));
    }
  }
  CHECK(n_cx <= 20);
  Eigen::MatrixXcd U1 = tket_sim::get_unitary(c);
  CHECK(tket_sim::compare_statevectors_or_unitaries(U, U1));
  Eigen::MatrixXcd U2 = get_3q_unitary(c);
  CHECK(tket_sim::compare_statevectors_or_unitaries(U, U2));
}

SCENARIO("Three-qubit circuits") {
  GIVEN("Round trip from a fixed 8x8 unitary") {
    // Randomly generated with scipy.stats.unitary_group.rvs.
    Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(8, 8);
    U(0, 0) = {0.13498182298658645, 0.07123133847729184};
    U(0, 1) = {-0.09573194703845724, -0.16102732488948038};
    U(0, 2) = {-0.12661200210472828, -0.46442090447446166};
    U(0, 3) = {-0.10016687393638907, 0.043966671778211466};
    U(0, 4) = {0.10455413844039041, -0.16217425591741186};
    U(0, 5) = {0.0633772123096529, -0.46717004599095074};
    U(0, 6) = {-0.3410867513707317, -0.37030410986723145};
    U(0, 7) = {-0.16778244522138677, 0.3959986334035144};
    U(1, 0) = {-0.23011745423147706, -0.003641966363325857};
    U(1, 1) = {0.03439963700002191, 0.061288252287784575};
    U(1, 2) = {0.015309416106435249, -0.4495585084322802};
    U(1, 3) = {0.3381657205950156, -0.1552501023149442};
    U(1, 4) = {-0.23941007652704477, -0.2573524278971172};
    U(1, 5) = {0.02545709649955111, 0.1830581970106428};
    U(1, 6) = {0.33852199408646483, 0.349686626423723};
    U(1, 7) = {-0.07832111105914577, 0.44786083899517665};
    U(2, 0) = {0.3147416482655379, -0.16247520773314159};
    U(2, 1) = {0.11224992323269103, 0.09442187933640739};
    U(2, 2) = {-0.06374407011546189, 0.09228924438291877};
    U(2, 3) = {0.1369446561864521, -0.3249687937032188};
    U(2, 4) = {-0.4922409153443705, -0.00965783139898696};
    U(2, 5) = {0.41360674134209213, -0.29229837481954035};
    U(2, 6) = {-0.2368018348299072, 0.24786480635214522};
    U(2, 7) = {-0.2350630029083468, -0.2107482928000711};
    U(3, 0) = {-0.14877249805087756, -0.2199731421948551};
    U(3, 1) = {0.04154897490233493, -0.251279632899556};
    U(3, 2) = {-0.3199118413634433, 0.16284854009182598};
    U(3, 3) = {0.2684142932688534, -0.028702106737609374};
    U(3, 4) = {-0.1547556915451631, -0.3280382237651188};
    U(3, 5) = {-0.3620989171490807, -0.10115751979391131};
    U(3, 6) = {0.25344952205706556, -0.3611206486270939};
    U(3, 7) = {-0.35302039564841275, -0.26589934096705037};
    U(4, 0) = {-0.2676622786385665, -0.41277317850222667};
    U(4, 1) = {-0.028377376885981128, -0.08427370381831906};
    U(4, 2) = {0.33827611224559684, -0.22103186906164543};
    U(4, 3) = {-0.42872033550794214, -0.1278382923865302};
    U(4, 4) = {-0.28539802330491526, -0.3406383162546306};
    U(4, 5) = {-0.007863981716221241, 0.03529699964705767};
    U(4, 6) = {-0.10847018386667288, -0.1233473498862589};
    U(4, 7) = {0.3447910329375026, -0.20489725275845655};
    U(5, 0) = {0.1875216586884449, -0.4977057328730517};
    U(5, 1) = {-0.4813140413884028, -0.33830562470128933};
    U(5, 2) = {-0.16158965857709123, -0.23235304921374156};
    U(5, 3) = {0.165439103189247, 0.0021286169941654096};
    U(5, 4) = {0.11061265636434567, 0.33680737520179843};
    U(5, 5) = {0.18844227595822957, 0.3134105313153881};
    U(5, 6) = {0.012041313388873987, -0.011704500833644174};
    U(5, 7) = {0.0028355715943290116, -0.061934598193077534};
    U(6, 0) = {0.3782021515664354, -0.06942953467359816};
    U(6, 1) = {-0.023349998945656567, 0.4265165406448352};
    U(6, 2) = {-0.2356282775967536, -0.16573704862370536};
    U(6, 3) = {-0.518189529099242, -0.07026032108380222};
    U(6, 4) = {0.00869795163707885, -0.1302526383364964};
    U(6, 5) = {-0.15910479099628008, 0.3448190658870644};
    U(6, 6) = {0.14166309775100366, -0.017596954959736386};
    U(6, 7) = {-0.36253361765167524, -0.013042170262335118};
    U(7, 0) = {0.11855668222327743, 0.20234520882323884};
    U(7, 1) = {-0.5017121940061847, -0.29220670756666767};
    U(7, 2) = {0.03271492594441479, 0.3115394934232476};
    U(7, 3) = {-0.2874418183011849, 0.26219877599860086};
    U(7, 4) = {-0.16506300521522488, -0.3127613010108169};
    U(7, 5) = {0.21010833148777708, -0.1493832562187949};
    U(7, 6) = {0.32103151296141164, 0.21632916190339818};
    U(7, 7) = {-0.04551679980729478, 0.12225020102655798};
    check_three_qubit_synthesis(U);
  }
  GIVEN("Round trip from a small circuit") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    auto u = tket_sim::get_unitary(c);
    check_three_qubit_synthesis(u);
  }
  GIVEN("Round trip from a larger circuit") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Y, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Y, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Y, {0});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Y, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Y, {0});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::Y, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Y, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::Y, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    auto u = tket_sim::get_unitary(c);
    check_three_qubit_synthesis(u);
  }
}

static void check_3q_unitary(const Circuit &c) {
  Eigen::MatrixXcd U = get_3q_unitary(c);
  Eigen::MatrixXcd U1 = tket_sim::get_unitary(c);
  CHECK(tket_sim::compare_statevectors_or_unitaries(U, U1));
}

SCENARIO("Unitary from circuits") {
  GIVEN("Empty circuit") {
    Circuit c(3);
    check_3q_unitary(c);
  }
  GIVEN("Simple circuit") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::S, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    check_3q_unitary(c);
  }
}

static bool check_3q_squash(const Circuit &c) {
  unsigned n_cx = c.count_gates(OpType::CX);
  Eigen::MatrixXcd U = tket_sim::get_unitary(c);
  Circuit c1 = c;
  bool success = Transforms::three_qubit_squash().apply(c1);
  unsigned n_cx1 = c1.count_gates(OpType::CX);
  if (success) {
    CHECK(n_cx1 < n_cx);
    Eigen::MatrixXcd U1 = tket_sim::get_unitary(c1);
    CHECK(tket_sim::compare_statevectors_or_unitaries(U, U1));
  } else {
    CHECK(c == c1);
  }
  return success;
}

SCENARIO("Three-qubit squash") {
  GIVEN("Empty circuit") {
    Circuit c(2);
    CHECK_FALSE(check_3q_squash(c));
  }
  GIVEN("1-qubit circuit, 1 gate") {
    Circuit c(1);
    c.add_op<unsigned>(OpType::H, {0});
    CHECK_FALSE(check_3q_squash(c));
  }
  GIVEN("1-qubit circuit, 2 gates") {
    Circuit c(1);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, 0.25, {0});
    CHECK_FALSE(check_3q_squash(c));
  }
  GIVEN("2-qubit circuit that cannot be squashed") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    CHECK_FALSE(check_3q_squash(c));
  }
  GIVEN("2-qubit circuit that can be squashed") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::Rz, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, 0.25, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::Rz, 0.25, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    CHECK(check_3q_squash(c));
  }
  GIVEN("3-qubit circuit that cannot be squashed") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Rz, 0.25, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Rz, 0.25, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    CHECK_FALSE(check_3q_squash(c));
  }
  GIVEN("Three-qubit circuit that can be squashed") {
    Circuit c(3);
    for (unsigned i = 0; i < 21; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    CHECK(check_3q_squash(c));
  }
  GIVEN("A complex circuit that can be squashed") {
    Circuit c2q(2);
    for (unsigned i = 0; i < 4; i++) {
      c2q.add_op<unsigned>(OpType::Rz, 0.25, {i % 2});
      c2q.add_op<unsigned>(OpType::CX, {i % 2, (i + 1) % 2});
    }
    CircBox c2qbox(c2q);
    Circuit c3q(3);
    for (unsigned i = 0; i < 21; i++) {
      c3q.add_op<unsigned>(OpType::Rz, 0.25, {i % 3});
      c3q.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
    }
    CircBox c3qbox(c3q);
    Circuit c(5);
    c.add_box(c2qbox, {1, 3});
    c.add_box(c3qbox, {3, 0, 2});
    c.add_box(c2qbox, {4, 2});
    c.add_box(c3qbox, {4, 3, 0});
    Transforms::decomp_boxes().apply(c);
    CHECK(check_3q_squash(c));
  }
  GIVEN("A circuit with measurements") {
    Circuit c(3, 3);
    for (unsigned i = 0; i < 22; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    for (unsigned q = 0; q < 3; q++) {
      c.add_op<unsigned>(OpType::Measure, {q, q});
    }
    CHECK(Transforms::three_qubit_squash().apply(c));
    CHECK(c.count_gates(OpType::CX) <= 20);
  }
  GIVEN("A circuit with classical control") {
    Circuit c(3, 1);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    for (unsigned i = 0; i < 11; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    c.add_conditional_gate<unsigned>(OpType::X, {}, {0}, {0}, 1);
    for (unsigned i = 11; i < 22; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    CHECK_FALSE(Transforms::three_qubit_squash().apply(c));
  }
  GIVEN("A circuit with a barrier") {
    Circuit c(3);
    for (unsigned i = 0; i < 11; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    c.add_barrier({0, 1, 2});
    for (unsigned i = 11; i < 22; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    CHECK_FALSE(Transforms::three_qubit_squash().apply(c));
  }
  GIVEN("A symbolic circuit") {
    Sym a = SymEngine::symbol("alpha");
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Ry, 2 * Expr(a), {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    CHECK(Transforms::three_qubit_squash().apply(c));
  }
}

SCENARIO("Special cases") {
  GIVEN("A 3-qubit circuit with no interaction between qb 0 and qbs 1,2") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::U3, {0.6, 0.7, 0.8}, {0});
    c.add_op<unsigned>(OpType::Rz, 0.1, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Rz, 0.2, {2});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    c.add_op<unsigned>(OpType::Rz, 0.3, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Rz, 0.4, {2});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    Eigen::MatrixXcd U = tket_sim::get_unitary(c);
    Circuit c1 = three_qubit_synthesis(U);
    REQUIRE(c1.count_gates(OpType::CX) <= 3);
    for (const Command &com : c1) {
      qubit_vector_t qbs = com.get_qubits();
      if (qbs.size() == 2) {
        CHECK(qbs[0] != Qubit(0));
        CHECK(qbs[1] != Qubit(0));
      }
    }
    Eigen::Matrix U1 = tket_sim::get_unitary(c1);
    CHECK(U.isApprox(U1));
  }
  GIVEN("A 3-qubit circuit with no interaction between qb 1 and qbs 0,2") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::U3, {0.6, 0.7, 0.8}, {1});
    c.add_op<unsigned>(OpType::Rz, 0.1, {0});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::Rz, 0.2, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::Rz, 0.3, {0});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::Rz, 0.4, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    Eigen::MatrixXcd U = tket_sim::get_unitary(c);
    Circuit c1 = three_qubit_synthesis(U);
    REQUIRE(c1.count_gates(OpType::CX) <= 3);
    for (const Command &com : c1) {
      qubit_vector_t qbs = com.get_qubits();
      if (qbs.size() == 2) {
        CHECK(qbs[0] != Qubit(1));
        CHECK(qbs[1] != Qubit(1));
      }
    }
    Eigen::Matrix U1 = tket_sim::get_unitary(c1);
    CHECK(U.isApprox(U1));
  }
  GIVEN("A 3-qubit circuit with no interaction between qb 0 and qbs 1,2") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::U3, {0.6, 0.7, 0.8}, {2});
    c.add_op<unsigned>(OpType::Rz, 0.1, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rz, 0.2, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Rz, 0.3, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rz, 0.4, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    Eigen::MatrixXcd U = tket_sim::get_unitary(c);
    Circuit c1 = three_qubit_synthesis(U);
    REQUIRE(c1.count_gates(OpType::CX) <= 3);
    for (const Command &com : c1) {
      qubit_vector_t qbs = com.get_qubits();
      if (qbs.size() == 2) {
        CHECK(qbs[0] != Qubit(2));
        CHECK(qbs[1] != Qubit(2));
      }
    }
    Eigen::Matrix U1 = tket_sim::get_unitary(c1);
    CHECK(U.isApprox(U1));
  }
}

}  // namespace test_ThreeQubitConversion
}  // namespace tket
