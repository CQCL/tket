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

#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <optional>

#include "Simulation/ComparisonFunctions.hpp"
#include "testutil.hpp"
#include "tket/Characterisation/DeviceCharacterisation.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Mapping/LexiLabelling.hpp"
#include "tket/Mapping/LexiRoute.hpp"
#include "tket/Mapping/MappingManager.hpp"
#include "tket/Mapping/Verification.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/Predicates/CompilerPass.hpp"
#include "tket/Predicates/PassGenerators.hpp"
#include "tket/Predicates/Predicates.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/OptimisationPass.hpp"
#include "tket/Transformations/Rebase.hpp"
#include "tket/Transformations/Transform.hpp"
#include "tket/Utils/HelperFunctions.hpp"

namespace tket {

using Connection = Architecture::Connection;

SCENARIO("Test decompose_SWAP_to_CX pass", ) {
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
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
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
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    Transforms::decompose_SWAP_to_CX(arc).apply(circ);
    Transforms::decompose_CX_directed(grid).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, grid, true));
  }
}

SCENARIO("Routing preserves the number of qubits in given instance of CnX op") {
  std::vector<std::pair<Node, Node>> cons;
  cons.push_back({Node("x", 1), Node("x", 0)});
  cons.push_back({Node("x", 2), Node("x", 1)});
  Architecture arc(
      std::vector<std::pair<Node, Node>>(cons.begin(), cons.end()));
  PassPtr pass = gen_default_mapping_pass(arc, false);
  Circuit c(3);
  c.add_op<unsigned>(OpType::CnX, {2, 1});
  CompilationUnit cu(c);
  pass->apply(cu);
  const Circuit &c1 = cu.get_circ_ref();
  REQUIRE(c.n_qubits() == c1.n_qubits());
}
SCENARIO(
    "Default mapping pass deals with circuits with only single qubit gates.") {
  std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
  Architecture architecture(edges);

  Circuit circuit(1);
  circuit.add_op<unsigned>(OpType::S, {0});
  PassPtr pass = gen_default_mapping_pass(architecture, true);
  CompilationUnit cu(circuit);
  bool applied = pass->apply(cu);
  REQUIRE(applied);

  const unit_bimap_t &initial_map = cu.get_initial_map_ref();
  const unit_bimap_t &final_map = cu.get_final_map_ref();
  REQUIRE(initial_map.left.size() == 1);
  REQUIRE(final_map.left.size() == 1);
  REQUIRE(initial_map.left.find(Qubit(0))->second == Node(0));
  REQUIRE(final_map.left.find(Qubit(0))->second == Node(0));
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
    "Methods related to correct routing and decomposition of circuits with "
    "classical wires.") {
  GIVEN("A circuit with classical wires on CX gates.") {
    Architecture test_arc({{0, 1}, {1, 2}});
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {2, 1}, {0, 1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 2);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {2, 1}, {1, 0}, 3);
    circ.add_conditional_barrier({1, 2}, {1}, {0}, 1, "");
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 2}, {0, 1}, 0);
    MappingManager mm(std::make_shared<Architecture>(test_arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));

    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, test_arc, false, false));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, test_arc, false, false));
  }
  GIVEN(
      "A circuit that requires modification to satisfy architecture "
      "constraints.") {
    Architecture arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    Circuit circ(5, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {1, 3}, {1, 4}, {0, 1}});
    circ.add_conditional_barrier({0, 1, 2}, {}, {0}, 1, "");

    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));

    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, false));
    Command classical_com = circ.get_commands()[0];
    REQUIRE(classical_com.get_args()[0] == circ.all_bits()[0]);
  }
  GIVEN("A single Bridge gate with multiple classical wires, decomposed.") {
    Architecture arc({{0, 1}, {1, 2}});
    Circuit circ(3, 3);
    circ.add_conditional_gate<unsigned>(
        OpType::BRIDGE, {}, {0, 1, 2}, {0, 1, 2}, 1);
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
    for (Command com : circ.get_commands()) {
      REQUIRE(com.get_args()[0] == circ.all_bits()[0]);
      REQUIRE(com.get_args()[1] == circ.all_bits()[1]);
      REQUIRE(com.get_args()[2] == circ.all_bits()[2]);
    }
  }
  GIVEN("A directed architecture, a single CX gate that requires flipping.") {
    Architecture arc(std::vector<std::pair<unsigned, unsigned>>{{0, 1}});
    Circuit circ(2, 2);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {1, 0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {1, 0}, {0, 1}, 1);
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, false));
    REQUIRE(!respects_connectivity_constraints(circ, arc, true, false));
    Transforms::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true, false));
    std::vector<Command> all_coms = circ.get_commands();
    REQUIRE(all_coms[0].get_args()[0] == circ.all_bits()[1]);
    REQUIRE(all_coms[0].get_args()[1] == circ.all_bits()[0]);
    REQUIRE(all_coms[1].get_args()[0] == circ.all_bits()[0]);
    REQUIRE(all_coms[1].get_args()[1] == circ.all_bits()[1]);
  }
}

SCENARIO(
    "Methods related to correct routing and decomposition of circuits with "
    "classical wires - long test.",
    "[.long]") {
  GIVEN(
      "A large circuit, with a mixture of conditional CX and CZ with "
      "multiple classical wires, non conditional CX and CZ, and single "
      "qubit gates.") {
    SquareGrid arc(5, 10);
    Circuit circ(50, 10);
    for (unsigned i = 0; i < 48; i++) {
      circ.add_op<unsigned>(OpType::CX, {i, i + 1});
      circ.add_conditional_gate<unsigned>(
          OpType::CX, {}, {i + 2, i}, {0, 2, 3, 5}, 1);
      circ.add_conditional_gate<unsigned>(OpType::H, {}, {i}, {0, 7}, 1);
      circ.add_conditional_gate<unsigned>(
          OpType::CX, {}, {i + 2, i + 1}, {1, 2, 3, 5, 9}, 0);
      circ.add_conditional_gate<unsigned>(OpType::S, {}, {i + 1}, {1, 2, 7}, 1);
      circ.add_conditional_gate<unsigned>(
          OpType::CZ, {}, {i, i + 1}, {4, 6, 8, 7, 9}, 0);
      circ.add_conditional_gate<unsigned>(OpType::X, {}, {i + 2}, {0, 3}, 0);
    }
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));

    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
  }
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
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    Transform T_1 = Transforms::decompose_SWAP_to_CX();
    T_1.apply(circ);
    REQUIRE(circ.count_gates(OpType::SWAP) == 0);
  }
}

SCENARIO("Test various large, unplaced, circuits with conditional gates.") {
  GIVEN("10x10x0 qubit Square Grid Architecture, 45 qubits, 5 bits.") {
    SquareGrid architecture(10, 10);
    Circuit circ(45, 5);
    std::vector<std::pair<unsigned, unsigned>> edges_0;
    for (unsigned i = 4; i < 38; i++) {
      edges_0.push_back({i, i + 4});
      if (38 - i != i - 4) {
        edges_0.push_back({38 - i, i - 4});
      }
    }
    add_2qb_gates(circ, OpType::CX, edges_0);
    std::vector<unsigned> barrier_args(44);
    std::iota(barrier_args.begin(), barrier_args.end(), 1);
    circ.add_barrier(barrier_args);
    add_2qb_gates(circ, OpType::CZ, edges_0);
    circ.add_measure(Qubit(q_default_reg(), 10), Bit(c_default_reg(), 0));
    circ.add_measure(Qubit(q_default_reg(), 12), Bit(c_default_reg(), 1));
    circ.add_measure(Qubit(q_default_reg(), 20), Bit(c_default_reg(), 2));
    circ.add_measure(Qubit(q_default_reg(), 18), Bit(c_default_reg(), 3));
    circ.add_measure(Qubit(q_default_reg(), 32), Bit(c_default_reg(), 4));
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {5, 10}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {7, 15}, {1, 0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {44, 3}, {2}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {42, 43}, {3}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {12, 13}, {4}, 1);
    circ.add_barrier(barrier_args);
    add_2qb_gates(circ, OpType::ZZMax, edges_0);
    MappingManager mm(std::make_shared<Architecture>(architecture));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(respects_connectivity_constraints(circ, architecture, false, true));
    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, architecture, false, true));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(
        respects_connectivity_constraints(circ, architecture, false, false));
  }
  GIVEN("39 qubit Ring Architecture, 23 qubits, 10 bits.") {
    RingArch architecture(39);
    Circuit circ(23, 10);
    std::vector<std::pair<unsigned, unsigned>> edges_0;
    for (unsigned i = 2; i < 20; i++) {
      edges_0.push_back({i, i + 1});
      if (21 - i != i - 2) {
        edges_0.push_back({21 - i, i - 2});
      }
    }
    add_2qb_gates(circ, OpType::CX, edges_0);
    std::vector<unsigned> barrier_args(18);
    std::iota(barrier_args.begin(), barrier_args.end(), 1);
    circ.add_barrier(barrier_args);
    add_2qb_gates(circ, OpType::ZZMax, edges_0);
    circ.add_measure(Qubit(q_default_reg(), 9), Bit(c_default_reg(), 0));
    circ.add_measure(Qubit(q_default_reg(), 12), Bit(c_default_reg(), 1));
    circ.add_measure(Qubit(q_default_reg(), 10), Bit(c_default_reg(), 2));
    circ.add_measure(Qubit(q_default_reg(), 11), Bit(c_default_reg(), 3));
    circ.add_measure(Qubit(q_default_reg(), 4), Bit(c_default_reg(), 4));
    circ.add_measure(Qubit(q_default_reg(), 2), Bit(c_default_reg(), 5));
    circ.add_measure(Qubit(q_default_reg(), 1), Bit(c_default_reg(), 6));
    circ.add_measure(Qubit(q_default_reg(), 20), Bit(c_default_reg(), 7));
    circ.add_measure(Qubit(q_default_reg(), 19), Bit(c_default_reg(), 8));
    circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {5, 10}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {7, 15}, {1, 0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {14, 3}, {2, 1, 0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::H, {}, {2}, {3, 4, 8}, 1);
    circ.add_conditional_gate<unsigned>(OpType::S, {}, {1}, {4, 5, 6, 7}, 1);
    circ.add_barrier(barrier_args);
    add_2qb_gates(circ, OpType::ZZMax, edges_0);
    MappingManager mm(std::make_shared<Architecture>(architecture));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(respects_connectivity_constraints(circ, architecture, false, true));
    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, architecture, false, true));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(
        respects_connectivity_constraints(circ, architecture, false, false));
  }
  GIVEN("3x3x5 qubit Square Grid Architecture, 14 qubits, 4 bits.") {
    SquareGrid architecture(3, 3, 5);
    Circuit circ(14, 4);
    std::vector<std::pair<unsigned, unsigned>> edges_0;
    for (unsigned i = 2; i < 10; i++) {
      edges_0.push_back({i, i + 1});
      if (10 - i != i - 2) {
        edges_0.push_back({10 - i, i - 2});
      }
      edges_0.push_back({i + 1, i});
    }
    add_2qb_gates(circ, OpType::CX, edges_0);
    std::vector<unsigned> barrier_args(10);
    std::iota(barrier_args.begin(), barrier_args.end(), 1);
    circ.add_barrier(barrier_args);
    add_2qb_gates(circ, OpType::ZZMax, edges_0);
    circ.add_measure(Qubit(q_default_reg(), 9), Bit(c_default_reg(), 0));
    circ.add_measure(Qubit(q_default_reg(), 12), Bit(c_default_reg(), 1));
    circ.add_measure(Qubit(q_default_reg(), 10), Bit(c_default_reg(), 2));
    circ.add_measure(Qubit(q_default_reg(), 11), Bit(c_default_reg(), 3));
    circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {5, 10}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {7, 2}, {1, 0}, 1);
    circ.add_conditional_gate<unsigned>(
        OpType::CZ, {}, {1, 3}, {2, 1, 0, 3}, 1);
    circ.add_barrier(barrier_args);
    circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {9, 8}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {8, 7}, {1, 0}, 1);
    circ.add_barrier(barrier_args);
    circ.add_conditional_gate<unsigned>(
        OpType::CZ, {}, {7, 6}, {2, 1, 0, 3}, 1);
    circ.add_barrier(barrier_args);
    circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {11, 0}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {12, 3}, {1, 0}, 1);
    add_2qb_gates(circ, OpType::CX, edges_0);
    MappingManager mm(std::make_shared<Architecture>(architecture));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(respects_connectivity_constraints(circ, architecture, false, true));
    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, architecture, false, true));
    Transforms::decompose_BRIDGE_to_CX().apply(circ);
    REQUIRE(
        respects_connectivity_constraints(circ, architecture, false, false));
  }
}

}  // namespace tket
