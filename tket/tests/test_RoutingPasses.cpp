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

#include <catch2/catch.hpp>
#include <numeric>
#include <optional>

#include "Characterisation/DeviceCharacterisation.hpp"
#include "Circuit/Circuit.hpp"
#include "OpType/OpType.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/Predicates.hpp"
#include "Mapping/MappingManager.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/Verification.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/HelperFunctions.hpp"
#include "testutil.hpp"

namespace tket {

using Connection = Architecture::Connection;

SCENARIO("Test decompose_SWAP_to_CX pass", "[routing]") {
  Architecture arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
  GIVEN("A single SWAP gate. Finding if correct number of vertices added") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    int original_vertices = circ.n_vertices();
    reassign_boundary(circ);
    Transforms::decompose_SWAP_to_CX().apply(circ);
    int decompose_vertices = circ.n_vertices();
    REQUIRE(decompose_vertices - original_vertices == 2);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A single SWAP gate, finding if correct path is preserved.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    // check output boundary
    Vertex boundary_0 = circ.get_out(Qubit(0));
    Vertex boundary_1 = circ.get_out(Qubit(1));
    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(circ.get_out(Qubit(0)) == boundary_0);
    REQUIRE(circ.get_out(Qubit(1)) == boundary_1);
    // check output boundary is the same
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Transforms::decompose_SWAP_to_CX().apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[0], all[1]};
    REQUIRE(circ.get_commands()[2].get_args() == cor);
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture, opposite case.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Transforms::decompose_SWAP_to_CX().apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[1], all[0]};
    REQUIRE(circ.get_commands()[2].get_args() == cor);
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture, opposite SWAP.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Transforms::decompose_SWAP_to_CX().apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[0], all[1]};
    REQUIRE(circ.get_commands()[2].get_args() == cor);
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture, opposite case, opposite SWAP.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    Transforms::decompose_SWAP_to_CX().apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[1], all[0]};
    REQUIRE(circ.get_commands()[2].get_args() == cor);
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture, opposite SWAP, pre CX.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    Transforms::decompose_SWAP_to_CX().apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[0], all[1]};
    REQUIRE(circ.get_commands()[1].get_args() == cor);
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture, opposite case, opposite SWAP, pre CX.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    Transforms::decompose_SWAP_to_CX().apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[1], all[0]};
    REQUIRE(circ.get_commands()[1].get_args() == cor);
  }
  GIVEN(
      "A circuit that facilitates some CX annihilation for an undirected "
      "architecture, opposite case, opposite SWAP, pre CX, directed bool "
      "on.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    reassign_boundary(circ);
    Transforms::decompose_SWAP_to_CX(arc).apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[1], all[0]};
    REQUIRE(circ.get_commands()[1].get_args() == cor);
  }
  GIVEN("A circuit that with no CX gates, but with directed architecture.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    reassign_boundary(circ);
    Transforms::decompose_SWAP_to_CX(arc).apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[0], all[1]};
    REQUIRE(circ.get_commands()[0].get_args() == cor);
  }
  GIVEN(
      "A circuit that with no CX gates, but with directed architecture, "
      "opposite case.") {
    Architecture dummy_arc({{1, 0}});
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    reassign_boundary(circ);
    Transforms::decompose_SWAP_to_CX(dummy_arc).apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[1], all[0]};
    REQUIRE(circ.get_commands()[0].get_args() == cor);
  }
  // TEST CIRCUIT
  Circuit circ(10);
  int count = 0;
  for (unsigned x = 0; x < 10; ++x) {
    for (unsigned y = 0; y + 1 < x; ++y) {
      count += 2;
      if (x % 2) {
        add_2qb_gates(circ, OpType::SWAP, {{x, y}, {y + 1, y}});
      } else {
        add_2qb_gates(circ, OpType::SWAP, {{y, x}, {y, y + 1}});
      }
    }
  }

  GIVEN("A network of SWAP gates.") {
    int original_vertices = circ.n_vertices();
    std::vector<Vertex> original_boundary;
    for (unsigned i = 0; i < circ.n_qubits(); i++) {
      original_boundary.push_back(circ.get_out(Qubit(i)));
    }
    Transforms::decompose_SWAP_to_CX().apply(circ);
    int decompose_vertices = circ.n_vertices();
    for (unsigned i = 0; i < circ.n_qubits(); i++) {
      REQUIRE(original_boundary[i] == circ.get_out(Qubit(i)));
    }
    REQUIRE(decompose_vertices - original_vertices == 2 * count);
  }
  GIVEN("A routed network of SWAP gates.") {
    SquareGrid grid(2, 5);
    MappingManager mm(std::make_shared<Architecture>(grid));
    REQUIRE(mm.route_circuit(circ, {std::make_shared<LexiRouteRoutingMethod>()}));
    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, grid, false, true));
    GIVEN("Directed CX gates") {
      Transforms::decompose_SWAP_to_CX().apply(circ);
      Transforms::decompose_BRIDGE_to_CX().apply(circ);
      Transforms::decompose_CX_directed(grid).apply(circ);
      REQUIRE(respects_connectivity_constraints(circ, grid, true));
    }
  }
}

SCENARIO("Test redirect_CX_gates pass", "[routing]") {
  Architecture arc({{1, 0}, {1, 2}});
  GIVEN("A circuit that requires no redirection.") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{1, 0}, {1, 2}});
    reassign_boundary(circ);
    Transforms::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A circuit that requires redirection.") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {2, 1}});
    reassign_boundary(circ);
    Transforms::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A circuit that requires no redirection, with SWAP.") {
    Circuit circ(3);

    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 1});

    swap_v = circ.add_op<unsigned>(OpType::SWAP, {0, 2});
    swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {2, 1});
    reassign_boundary(circ);
    Transforms::decompose_SWAP_to_CX(arc).apply(circ);
    Transforms::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A circuit that requires redirection, with SWAP.") {
    Circuit circ(3);

    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {1, 0});

    swap_v = circ.add_op<unsigned>(OpType::SWAP, {0, 2});
    swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {1, 2});

    reassign_boundary(circ);
    Transforms::decompose_SWAP_to_CX(arc).apply(circ);
    Transforms::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A complicated circuit of CX gates, routed.") {
    Circuit circ(12);
    SquareGrid grid(3, 4);

    for (unsigned x = 0; x < 12; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    MappingManager mm(std::make_shared<Architecture>(grid));
    REQUIRE(mm.route_circuit(circ, {std::make_shared<LexiRouteRoutingMethod>()}));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    Transforms::decompose_SWAP_to_CX(arc).apply(circ);
    Transforms::decompose_CX_directed(grid).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, grid, true));
  }
}



SCENARIO("Routing preserves the number of qubits") {
  std::vector<std::pair<Node, Node>> cons;
  cons.push_back({Node("x", 1), Node("x", 0)});
  cons.push_back({Node("x", 2), Node("x", 1)});
  Architecture arc(
      std::vector<std::pair<Node, Node>>(cons.begin(), cons.end()));
  PassPtr pass = gen_default_mapping_pass(arc, false);
  Circuit c(3);
  c.add_op<unsigned>(OpType::CnX, {2, 1});
  CompilationUnit cu(c);
  bool applied = pass->apply(cu);
  const Circuit &c1 = cu.get_circ_ref();
  REQUIRE(c.n_qubits() == c1.n_qubits());
}

SCENARIO("Default mapping pass delays measurements") {
  std::vector<std::pair<Node, Node>> cons;
  cons.push_back({Node("x", 0), Node("x", 2)});
  cons.push_back({Node("x", 1), Node("x", 2)});
  cons.push_back({Node("x", 2), Node("x", 3)});
  cons.push_back({Node("x", 3), Node("x", 0)});
  Architecture arc(
      std::vector<std::pair<Node, Node>>(cons.begin(), cons.end()));
  PassPtr pass = gen_default_mapping_pass(arc, false);
  Circuit c(4, 4);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::CX, {1, 2});
  c.add_op<unsigned>(OpType::CX, {2, 3});
  c.add_op<unsigned>(OpType::CX, {3, 0});
  for (unsigned nn = 0; nn <= 3; ++nn) {
    c.add_measure(nn, nn);
  }
  Circuit c2(c);
  CompilationUnit cu(c);
  REQUIRE(pass->apply(cu));
  CompilationUnit cu2(c2);
  // delay_measures is default to true
  PassPtr pass2 = gen_default_mapping_pass(arc);
  REQUIRE(pass2->apply(cu2));
  PredicatePtr mid_meas_pred = std::make_shared<NoMidMeasurePredicate>();
  REQUIRE(!mid_meas_pred->verify(cu.get_circ_ref()));
  REQUIRE(mid_meas_pred->verify(cu2.get_circ_ref()));
}
SCENARIO(
    "Does copying decompose_SWAP_to_CX pass and applying it to a routed "
    "Circuit work correctly?") {
  GIVEN("A simple circuit and architecture.") {
    Circuit circ(5);
    add_2qb_gates(
        circ, OpType::CX,
        {{0, 3},
         {1, 4},
         {0, 1},
         {2, 0},
         {2, 1},
         {1, 0},
         {0, 4},
         {2, 1},
         {0, 3}});
    Architecture arc({{1, 0}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {4, 3}});
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(circ, {std::make_shared<LexiRouteRoutingMethod>()}));
    
    Transform T_1 = Transforms::decompose_SWAP_to_CX();
    T_1.apply(circ);
    REQUIRE(circ.count_gates(OpType::SWAP) == 0);
  }
}
}  // namespace tket
