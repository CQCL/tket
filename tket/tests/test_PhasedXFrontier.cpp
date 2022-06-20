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
#include <optional>

#include "Circuit/Circuit.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/PhasedXFrontier.hpp"
#include "Utils/Expression.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_PhasedXFrontier {

SCENARIO("PhasedXFrontier: move frontier forward") {
  GIVEN("A small circuit") {
    Circuit circ(2);
    Vertex v1 = circ.add_op<unsigned>(OpType::PhasedX, {0.3, 0.5}, {0});
    Vertex v2 = circ.add_op<unsigned>(OpType::PhasedX, {0.2, 0.3}, {1});
    Vertex v3 = circ.add_op<unsigned>(OpType::PhasedX, {0.3, 0.5}, {0});
    Vertex cx = circ.add_op<unsigned>(OpType::CX, {0, 1});
    Vertex v4 = circ.add_op<unsigned>(OpType::PhasedX, {0.3, 0.5}, {0});

    Transforms::PhasedXFrontier frontier(circ);

    std::vector<std::optional<Vertex>> exp{v1, v2};
    REQUIRE(frontier.get_all_beta_vertices() == exp);

    frontier.squash_intervals();

    auto vs = frontier.get_all_beta_vertices();
    REQUIRE(vs.size() == 2);
    REQUIRE(vs[0].has_value());
    REQUIRE(vs[1] == v2);
    Op_ptr op = circ.get_Op_ptr_from_Vertex(*vs[0]);
    REQUIRE(op->get_type() == OpType::PhasedX);
    REQUIRE(equiv_expr(op->get_params()[0], 0.6, 4));

    // move forward
    REQUIRE(frontier.are_phasedx_left());
    frontier.next_multiqb(cx);
    REQUIRE(!frontier.are_phasedx_left());

    exp = {v4, std::nullopt};
    REQUIRE(frontier.get_all_beta_vertices() == exp);

    frontier.next_interval(0);
    exp = {std::nullopt, std::nullopt};
    REQUIRE(frontier.get_all_beta_vertices() == exp);
  }
  GIVEN("A larger circuit") {
    Circuit circ(4);
    //        ---------   --------
    // q0: ---|       |---|      |---C------------------
    //        |       |   |      |   |
    // q1: ---| NPhX  |---o- - - o---|------------------
    //        |       |   | NphX |   |   --------
    // q2: ---|       |---o- - - o---Z---|      |---C---
    //        ---------   |      |       | NPhX |   |
    // q3: ---------------|      |-------|      |---Z---
    //                    --------       --------
    Vertex v1 = circ.add_op<unsigned>(OpType::NPhasedX, {0.3, 0.5}, {0, 1, 2});
    Vertex v2 = circ.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {0, 3});
    Vertex cz1 = circ.add_op<unsigned>(OpType::CZ, {0, 2});
    Vertex v3 = circ.add_op<unsigned>(OpType::NPhasedX, {0.6, 1.3}, {2, 3});
    Vertex cz2 = circ.add_op<unsigned>(OpType::CZ, {2, 3});

    Transforms::PhasedXFrontier frontier(circ);

    WHEN("Getting the initial PhasedX frontier") {
      REQUIRE(frontier.get_all_betas() == std::vector<Expr>{0.3, 0.3, 0.3, 0.});
      auto vs = frontier.get_all_beta_vertices();
      // qubit 3 is shadowed
      std::vector<std::optional<Vertex>> exp = {v1, v1, v1, std::nullopt};
      THEN("There are as many phases as qubits") {
        REQUIRE(vs.size() == exp.size());
      }
      THEN("The vertices match expecations") { REQUIRE(vs == exp); }
    }
    WHEN("Before insertion") {
      REQUIRE(circ.count_gates(OpType::NPhasedX) == 3);
    }
    frontier.insert_2_phasedx();

    WHEN("After insertion") {
      REQUIRE(circ.count_gates(OpType::NPhasedX) == 4);
    }
    WHEN("Getting the current PhasedX frontier") {
      REQUIRE(frontier.get_all_betas() == std::vector<Expr>{0.2, 0., 0., 0.2});
      auto vs = frontier.get_all_beta_vertices();
      std::vector<std::optional<Vertex>> exp = {
          v2, std::nullopt, std::nullopt, v2};
      THEN("There are as many phases as qubits") {
        REQUIRE(vs.size() == exp.size());
      }
      THEN("The vertices match expecations") { REQUIRE(vs == exp); }
    }

    frontier.insert_2_phasedx();

    WHEN("Getting the PhasedX frontier with no remaining NPhasedX") {
      REQUIRE(frontier.get_all_betas() == std::vector<Expr>{0., 0., 0., 0.});
      auto vs = frontier.get_all_beta_vertices();
      std::vector<std::optional<Vertex>> exp = {
          std::nullopt, std::nullopt, std::nullopt, std::nullopt};
      THEN("There are as many phases as qubits") {
        REQUIRE(vs.size() == exp.size());
      }
      THEN("The vertices match expecations") { REQUIRE(vs == exp); }
    }

    frontier.next_multiqb(cz1);

    WHEN("Getting the second current PhasedX frontier") {
      REQUIRE(frontier.get_all_betas() == std::vector<Expr>{0., 0., 0.6, 0.6});
      auto vs = frontier.get_all_beta_vertices();
      std::vector<std::optional<Vertex>> exp = {
          std::nullopt, std::nullopt, v3, v3};
      THEN("There are as many phases as qubits") {
        REQUIRE(vs.size() == exp.size());
      }
      THEN("The vertices match expecations") { REQUIRE(vs == exp); }
    }
    frontier.insert_2_phasedx();

    WHEN("Getting the PhasedX frontier with no remaining NPhasedX") {
      REQUIRE(frontier.get_all_betas() == std::vector<Expr>{0., 0., 0., 0.});
      auto vs = frontier.get_all_beta_vertices();
      std::vector<std::optional<Vertex>> exp = {
          std::nullopt, std::nullopt, std::nullopt, std::nullopt};
      THEN("There are as many phases as qubits") {
        REQUIRE(vs.size() == exp.size());
      }
      THEN("The vertices match expecations") { REQUIRE(vs == exp); }
    }

    frontier.next_multiqb(cz2);
    WHEN("Getting the final PhasedX frontier") {
      REQUIRE(frontier.get_all_betas() == std::vector<Expr>{0., 0., 0., 0.});
      auto vs = frontier.get_all_beta_vertices();
      std::vector<std::optional<Vertex>> exp = {
          std::nullopt, std::nullopt, std::nullopt, std::nullopt};
      THEN("There are as many phases as qubits") {
        REQUIRE(vs.size() == exp.size());
      }
      THEN("The vertices match expecations") { REQUIRE(vs == exp); }
    }
  }
}

SCENARIO("PhasedXFrontier: replacing gates") {
  GIVEN("3x PhasedX gates in parallel") {
    Circuit c1(3);
    c1.add_op<unsigned>(OpType::PhasedX, {0.3, 0.5}, {0});
    c1.add_op<unsigned>(OpType::PhasedX, {0.2, 0.3}, {1});
    c1.add_op<unsigned>(OpType::PhasedX, {0.3, 0.5}, {2});

    WHEN("Replacing with 2x NPhasedX") {
      Circuit c2 = c1;
      Transforms::PhasedXFrontier frontier(c2);
      frontier.insert_2_phasedx();
      THEN("There are a total of 2x (global) NPhasedX") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        REQUIRE(c2.count_gates(OpType::PhasedX) == 0);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Replacing q1 with 1x NPhasedX") {
      Circuit c2 = c1;
      Transforms::PhasedXFrontier frontier(c2);
      frontier.insert_1_phasedx(0);
      THEN("There are the expected number of PhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 1);
        REQUIRE(c2.count_gates(OpType::PhasedX) == 1);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
      WHEN("Replacing q2 with 1x NPhasedX") {
        frontier.insert_1_phasedx(1);
        THEN("There are the expected number of PhasedX gates") {
          REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
          REQUIRE(c2.count_gates(OpType::PhasedX) == 2);
        }
        THEN("The resulting unitaries are equal") {
          auto u1 = tket_sim::get_unitary(c1);
          auto u2 = tket_sim::get_unitary(c2);
          REQUIRE(u1.isApprox(u2));
        }
      }
    }
  }
  GIVEN("2x 2-NPhasedX in parallel") {
    Circuit c1(4);
    c1.add_op<unsigned>(OpType::NPhasedX, {0.3, 0.5}, {0, 1});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {2, 3});

    WHEN("Replacing q1-q2 with 2x global NPhasedX") {
      Circuit c2 = c1;
      Transforms::PhasedXFrontier frontier(c2);
      frontier.insert_2_phasedx();
      THEN("There are a total of 2x (global) NPhasedX") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        REQUIRE(c2.count_gates(OpType::PhasedX) == 0);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Replacing q1-q2 with 1x NPhasedX") {
      Circuit c2 = c1;
      Transforms::PhasedXFrontier frontier(c2);
      frontier.insert_1_phasedx(0);
      THEN("There are the expected number of PhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        REQUIRE(c2.count_gates(OpType::PhasedX) == 0);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("2x 2-NPhasedX behind eachother") {
    Circuit c1(3);
    c1.add_op<unsigned>(OpType::NPhasedX, {0.3, 0.5}, {0, 1});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {1, 2});

    WHEN("Replacing first one with 1x global NPhasedX") {
      Circuit c2 = c1;
      Transforms::PhasedXFrontier frontier(c2);
      frontier.insert_1_phasedx(0);
      THEN("There are a total of 2x NPhasedX and 1x PhasedX") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        REQUIRE(c2.count_gates(OpType::PhasedX) == 1);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
      WHEN("Replacing second one with 2x NPhasedX") {
        frontier.insert_2_phasedx();
        frontier.insert_2_phasedx();
        THEN("There are 3x global PhasedX gates") {
          REQUIRE(c2.count_gates(OpType::NPhasedX) == 5);
          REQUIRE(c2.count_gates(OpType::PhasedX) == 0);
        }
        THEN("The resulting unitaries are equal") {
          auto u1 = tket_sim::get_unitary(c1);
          auto u2 = tket_sim::get_unitary(c2);
          REQUIRE(u1.isApprox(u2));
        }
      }
    }
  }
}

SCENARIO("PhasedXFrontier: Squashing PhasedX") {
  GIVEN("A simple circuit") {
    Circuit c1(1);
    c1.add_op<unsigned>(OpType::PhasedX, {0.4, 0.3}, {0});
    c1.add_op<unsigned>(OpType::PhasedX, {0.4, 0.3}, {0});
    Circuit c2 = c1;

    Transforms::PhasedXFrontier frontier(c2);
    WHEN("Squashing the interval") {
      frontier.squash_intervals();
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
}

}  // namespace test_PhasedXFrontier
}  // namespace tket
