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
#include <optional>

#include "Characterisation/DeviceCharacterisation.hpp"
#include "Circuit/Circuit.hpp"
#include "OpType/OpType.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/Predicates.hpp"
#include "Routing/Routing.hpp"
#include "Routing/Verification.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/HelperFunctions.hpp"
#include "testutil.hpp"

namespace tket {

using Connection = Architecture::Connection;

Interactions RoutingTester::get_interaction(const RoutingFrontier &sf) {
  return router->generate_interaction_frontier(sf);
}

// Wrappers of private methods for testing?
void RoutingTester::increment_distance(
    graphs::dist_vec &new_dist_vector, const Swap &pair, int increment) const {
  router->increment_distance(new_dist_vector, pair, increment);
}

graphs::dist_vec RoutingTester::generate_distance_vector(
    const Interactions &inter) const {
  return router->generate_distance_vector(inter);
}

graphs::dist_vec RoutingTester::update_distance_vector(
    const Swap &nodes, graphs::dist_vec new_dist_vector,
    const Interactions &inte) const {
  return router->update_distance_vector(nodes, new_dist_vector, inte);
}

const std::pair<unsigned, unsigned> RoutingTester::pair_dists(
    const Node &n1, const Node &p1, const Node &n2, const Node &p2) const {
  return router->pair_dists(n1, p1, n2, p2);
}

bool RoutingTester::swap_decreases(
    const Swap &nodes, const Interactions &inte) const {
  return router->swap_decreases(nodes, inte);
}

std::vector<Swap> RoutingTester::candidate_swaps(
    const std::vector<Connection> &trial_edges,
    const Interactions &inte) const {
  return router->candidate_swaps(trial_edges, inte);
}

std::vector<Swap> RoutingTester::cowtan_et_al_heuristic(
    std::vector<Swap> &candidate_swaps, const graphs::dist_vec &base_dists,
    const Interactions &interac) const {
  return router->cowtan_et_al_heuristic(candidate_swaps, base_dists, interac);
}

void RoutingTester::update_qmap(qubit_bimap_t &map, const Swap &swap) {
  router->update_qmap(map, swap);
}

std::vector<Swap> RoutingTester::path_to_swaps(
    const std::vector<Node> &path) const {
  return router->path_to_swaps(path);
}

qubit_bimap_t default_qubit_map(const Circuit &circ) {
  qubit_bimap_t qmap;
  unsigned node = 0;
  for (const Qubit &qb : circ.all_qubits()) {
    qmap.insert({qb, Node(node)});
    node++;
  }
  return qmap;
}
qubit_bimap_t RoutingTester::set_default_initial_map(
    std::optional<node_vector_t> canonical_node_order) {
  qubit_bimap_t qmap;
  unsigned node = 0;
  for (const Qubit &qb : router->circ_.all_qubits()) {
    if (canonical_node_order.has_value()) {
      qmap.insert({qb, canonical_node_order->at(node)});
    } else {
      qmap.insert({qb, Node(node)});
    }
    node++;
  }
  router->init_map = qmap;
  router->qmap = qmap;
  return qmap;
}

void RoutingTester::initialise_slicefrontier() {
  router->slice_frontier_.init();
}

void RoutingTester::add_distributed_cx(
    const Node &control_node, const Node &target_node,
    const Node &central_node) {
  router->add_distributed_cx(control_node, target_node, central_node);
}

std::pair<std::pair<bool, Node>, std::pair<bool, Node>>
RoutingTester::check_distributed_cx(const Swap &nodes) {
  return router->check_distributed_cx(nodes);
}

void RoutingTester::advance_frontier() { router->advance_frontier(); }

void RoutingTester::set_interaction() {
  router->interaction =
      router->generate_interaction_frontier(router->slice_frontier_);
}
void RoutingTester::set_qmap(qubit_bimap_t _qmap) { router->qmap = _qmap; }
void RoutingTester::set_config(const RoutingConfig &_config) {
  router->config_ = _config;
}
void RoutingTester::next_sf(RoutingFrontier &sf) { sf.next_slicefrontier(); }
Circuit *RoutingTester::get_circ() { return &(router->circ_); }

namespace test_Routing {

SCENARIO(
    "Test validity of circuit against architecture using "
    "respects_connectivity_constraints method.",
    "[routing]") {
  Architecture arc({{1, 0}, {1, 2}});

  GIVEN("A simple CX circuit and a line_placement map.") {
    tket::Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {0, 3}, {2, 4}, {1, 4}, {0, 4}});
    tket::Architecture test_arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    LinePlacement lp_obj(test_arc);
    // qubit_mapping_t lm = lp_obj.place_get_maps(circ)[0];
    lp_obj.place(circ);
    tket::Routing router(circ, test_arc);
    std::pair<tket::Circuit, bool> outcirc = router.solve();
    REQUIRE(outcirc.second == true);
    CHECK(respects_connectivity_constraints(outcirc.first, test_arc, false));
  }
  GIVEN("A failing case, undirected") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A working case, undirected") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A failing case, directed") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A working case, directed") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A failing case, undirected, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 2});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A working case, undirected, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 2});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 2});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A failing case, directed, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 0});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {1, 0});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A working case, directed, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 0});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
}

SCENARIO("Test decompose_SWAP_to_CX pass", "[routing]") {
  Architecture arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
  GIVEN("A single SWAP gate. Finding if correct number of vertices added") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    int original_vertices = circ.n_vertices();
    reassign_boundary(circ);
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
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
    Transform::decompose_SWAP_to_CX(arc).apply(circ);
    qubit_vector_t all = circ.all_qubits();
    unit_vector_t cor = {all[1], all[0]};
    REQUIRE(circ.get_commands()[1].get_args() == cor);
  }
  GIVEN("A circuit that with no CX gates, but with directed architecture.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {1, 0});
    reassign_boundary(circ);
    Transform::decompose_SWAP_to_CX(arc).apply(circ);
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
    Transform::decompose_SWAP_to_CX(dummy_arc).apply(circ);
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
    Transform::decompose_SWAP_to_CX().apply(circ);
    int decompose_vertices = circ.n_vertices();
    for (unsigned i = 0; i < circ.n_qubits(); i++) {
      REQUIRE(original_boundary[i] == circ.get_out(Qubit(i)));
    }
    REQUIRE(decompose_vertices - original_vertices == 2 * count);
  }
  GIVEN("A routed network of SWAP gates.") {
    SquareGrid grid(2, 5);
    Routing router(circ, grid);
    std::pair<Circuit, bool> output = router.solve();
    REQUIRE(output.second);
    circ = output.first;
    Transform::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, grid, false, true));
    GIVEN("Directed CX gates") {
      Transform::decompose_SWAP_to_CX().apply(output.first);
      Transform::decompose_BRIDGE_to_CX().apply(output.first);
      Transform::decompose_CX_directed(grid).apply(output.first);
      REQUIRE(respects_connectivity_constraints(output.first, grid, true));
    }
  }
}

SCENARIO("Test redirect_CX_gates pass", "[routing]") {
  Architecture arc({{1, 0}, {1, 2}});
  GIVEN("A circuit that requires no redirection.") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{1, 0}, {1, 2}});
    reassign_boundary(circ);
    Transform::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A circuit that requires redirection.") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {2, 1}});
    reassign_boundary(circ);
    Transform::decompose_CX_directed(arc).apply(circ);
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
    Transform::decompose_SWAP_to_CX(arc).apply(circ);
    Transform::decompose_CX_directed(arc).apply(circ);
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
    Transform::decompose_SWAP_to_CX(arc).apply(circ);
    Transform::decompose_CX_directed(arc).apply(circ);
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
    Routing route(circ, grid);
    std::pair<Circuit, bool> outs = route.solve();
    REQUIRE(outs.second == true);
    circ = outs.first;
    Transform::decompose_BRIDGE_to_CX().apply(circ);
    Transform::decompose_SWAP_to_CX(arc).apply(circ);
    Transform::decompose_CX_directed(grid).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, grid, true));
  }
}

SCENARIO("Test RoutingFrontiers and interaction vectors", "[routing]") {
  GIVEN("A simple circuit") {
    Circuit incirc(4);
    Vertex v1 = incirc.add_op<unsigned>(OpType::X, {0});
    Vertex v8 = incirc.add_op<unsigned>(OpType::S, {3});
    Vertex v9 = incirc.add_op<unsigned>(OpType::T, {3});
    Vertex v2 = incirc.add_op<unsigned>(OpType::CX, {0, 1});
    Vertex v3 = incirc.add_op<unsigned>(OpType::CY, {2, 3});
    Vertex v4 = incirc.add_op<unsigned>(OpType::H, {0});
    Vertex v10 = incirc.add_op<unsigned>(OpType::X, {0});
    Vertex v11 = incirc.add_op<unsigned>(OpType::S, {1});
    Vertex v12 = incirc.add_op<unsigned>(OpType::Z, {3});
    Vertex v13 = incirc.add_op<unsigned>(OpType::Y, {2});
    Vertex v14 = incirc.add_op<unsigned>(OpType::T, {1});
    Vertex v5 = incirc.add_op<unsigned>(OpType::CZ, {0, 2});
    Vertex v6 = incirc.add_op<unsigned>(OpType::Y, {0});
    Vertex v7 = incirc.add_op<unsigned>(OpType::CX, {3, 1});

    // Ring of size 4
    RingArch arc(4);
    node_vector_t ring_nodes = RingArch::get_nodes_canonical_order(4);
    // Create Routing Object
    Routing router(incirc, arc);
    RoutingTester tester(&router);
    Circuit *circ = tester.get_circ();
    RoutingFrontier sf1 = router.get_slicefrontier();
    Qubit qb0(0);
    Qubit qb1(1);
    Qubit qb2(2);
    Qubit qb3(3);
    qubit_bimap_t qm;
    for (unsigned i = 0; i < 4; ++i) {
      qm.insert({Qubit(i), ring_nodes[i]});
    }
    tester.set_qmap(qm);
    WHEN("First interaction vector is generated") {
      Interactions inte = tester.get_interaction(sf1);
      THEN("Interaction vector is correct") {
        CHECK(inte[ring_nodes.at(0)] == ring_nodes.at(1));
        CHECK(inte[ring_nodes.at(1)] == ring_nodes.at(0));
        CHECK(inte[ring_nodes.at(3)] == ring_nodes.at(2));
        CHECK(inte[ring_nodes.at(2)] == ring_nodes.at(3));
        REQUIRE(inte.size() == 4);
      }
    }
    WHEN("One operation is completed") {
      Edge new_0 = circ->skip_irrelevant_edges(circ->get_all_out_edges(v2)[0]);
      Edge new_1 = circ->skip_irrelevant_edges(circ->get_all_out_edges(v2)[1]);
      sf1.quantum_in_edges->replace(
          sf1.quantum_in_edges->find(qb0), {qb0, new_0});
      sf1.quantum_in_edges->replace(
          sf1.quantum_in_edges->find(qb1), {qb1, new_1});
      CutFrontier next_cut = circ->next_cut(
          sf1.quantum_in_edges, std::make_shared<b_frontier_t>());

      sf1.slice = next_cut.slice;
      sf1.quantum_out_edges = next_cut.u_frontier;
      Interactions inte = tester.get_interaction(sf1);
      THEN("Interaction vector is updated") {
        CHECK(inte[ring_nodes.at(0)] == ring_nodes.at(0));
        CHECK(inte[ring_nodes.at(1)] == ring_nodes.at(1));
        CHECK(inte[ring_nodes.at(3)] == ring_nodes.at(2));
        CHECK(inte[ring_nodes.at(2)] == ring_nodes.at(3));
        REQUIRE(inte.size() == 4);
      }
    }

    WHEN("Next RoutingFrontier is generated") {
      sf1.next_slicefrontier();
      THEN("The RoutingFrontier is correct") {
        REQUIRE(sf1.slice->size() == 2);
        CHECK(
            circ->get_Op_ptr_from_Vertex(sf1.slice->at(0)) ==
            incirc.get_Op_ptr_from_Vertex(v5));
        CHECK(
            circ->get_Op_ptr_from_Vertex(sf1.slice->at(1)) ==
            incirc.get_Op_ptr_from_Vertex(v7));

        CHECK(
            sf1.quantum_in_edges->find(qb1)->second !=
            circ->get_nth_out_edge(v2, 1));
        CHECK(
            sf1.quantum_in_edges->find(qb2)->second ==
            circ->get_nth_in_edge(sf1.slice->at(0), 1));

        CHECK(
            sf1.quantum_out_edges->find(qb0)->second !=
            circ->get_nth_in_edge(v6, 0));
        CHECK(
            sf1.quantum_out_edges->find(qb3)->second ==
            circ->get_nth_out_edge(sf1.slice->at(1), 0));
      }
      sf1.next_slicefrontier();
      REQUIRE(sf1.slice->empty());
    }
  }
}

SCENARIO(
    "Check that an already solved routing problem will not add unecessary "
    "swaps",
    "[routing]") {
  GIVEN("A solved problem") {
    // Test Circuit, sequential cxs on a ring, requires no routing
    Circuit test_circuit;
    test_circuit.add_blank_wires(4);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});

    // Ring of size 4
    RingArch arc(4);
    // Create Routing Object
    Routing router(test_circuit, arc);
    std::pair<Circuit, bool> post_c = router.solve();
    REQUIRE(post_c.second == true);
    REQUIRE(post_c.first.n_gates() == 4);
  }
  GIVEN("A solved problem supplied with map and custom architecture") {
    Circuit test_circuit;
    test_circuit.add_blank_wires(4);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});

    Architecture test_arc({{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    Placement test_p(test_arc);

    qubit_mapping_t map_;
    for (unsigned nn = 0; nn <= 3; ++nn) {
      map_[Qubit(nn)] = Node(nn);
    }
    test_p.place_with_map(test_circuit, map_);
    qubit_vector_t all_qs_post_place = test_circuit.all_qubits();
    Routing router(test_circuit, test_arc);
    std::pair<Circuit, bool> result = router.solve();
    qubit_vector_t all_qs_post_solve = test_circuit.all_qubits();

    REQUIRE(all_qs_post_place == all_qs_post_solve);
    REQUIRE(result.second == false);
    REQUIRE(result.first.n_gates() == 4);
  }
}

SCENARIO(
    "If a circuit has fewer qubits than the architecture has nodes, is a "
    "correct sub-architecture made",
    "[routing]") {
  GIVEN("A circuit and architecture obeying said scenario") {
    // 5 wires, all used
    Circuit test_circuit(5);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 4}, {2, 3}, {1, 4}});

    SquareGrid arc(3, 3);
    Routing route(test_circuit, arc);
    route.solve();
    node_vector_t nodes = route.get_active_nodes();

    REQUIRE(nodes.size() == 5);

    // 5 wires, 4 used
    Circuit test_circuit2(5);
    add_2qb_gates(test_circuit2, OpType::CX, {{0, 3}, {1, 2}});

    Routing route2(test_circuit2, arc);
    route2.solve();
    node_vector_t nodes2 = route.get_active_nodes();

    REQUIRE(nodes2.size() == 5);
  }
}

SCENARIO("Qubit activating edge case", "[routing]") {
  GIVEN("A node line with only 3 qubits line placed") {
    Circuit circ;
    circ.add_blank_wires(4);
    add_2qb_gates(
        circ, OpType::CX, {{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {1, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.25, {2, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {2, 1});
    circ.add_op<unsigned>(OpType::CU1, 0.125, {3, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.25, {3, 1});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {3, 2});
    Transform::rebase_tket().apply(circ);
    Architecture arc({{0, 1}, {1, 2}, {2, 3}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> c = router.solve();
    REQUIRE(respects_connectivity_constraints(c.first, arc, false, true));
    REQUIRE(c.second);
  }
}

SCENARIO("Empty Circuit test", "[routing]") {
  GIVEN("An Empty Circuit") {
    Circuit circ;
    circ.add_blank_wires(4);
    Architecture arc({{0, 1}, {1, 2}, {2, 3}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.first.n_gates() == 0);
    REQUIRE(result.second == true);
    REQUIRE(respects_connectivity_constraints(result.first, arc, true));
  }
}

SCENARIO("Routing on circuit with no multi-qubit gates", "[routing]") {
  GIVEN("A circuit with no multi-qubit gates") {
    Circuit circ;
    circ.add_blank_wires(4);
    add_1qb_gates(circ, OpType::X, {0, 2});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::Y, {1});
    // circ.add_op<unsigned>(OpType::Y,{3});
    Architecture arc({{0, 1}, {1, 2}, {2, 3}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(circ.n_vertices() - 8 == result.first.n_gates());
    REQUIRE(result.second == true);
    REQUIRE(respects_connectivity_constraints(result.first, arc, true));
  }
}

SCENARIO("Test routing for other multi-qubit ops", "[routing]") {
  GIVEN("Failed qft circuit") {
    Circuit circ(4, 4);
    add_1qb_gates(circ, OpType::X, {0, 2});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {1, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {0, 1});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CU1, 0.25, {2, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {2, 1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::CU1, 0.125, {3, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.25, {3, 1});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {3, 2});
    circ.add_op<unsigned>(OpType::H, {3});
    for (unsigned nn = 0; nn <= 3; ++nn) {
      circ.add_measure(nn, nn);
    }
    Transform::rebase_tket().apply(circ);
    Architecture arc({{0, 1}, {1, 2}, {2, 3}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();

    REQUIRE(respects_connectivity_constraints(result.first, arc, false, true));
    REQUIRE(result.second);
  }
}

SCENARIO(
    "Test routing on a directed architecture with bidirectional edges",
    "[routing]") {
  GIVEN("A simple two-qubit circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Architecture arc({{0, 1}, {1, 0}});
    Architecture arc2(std::vector<std::pair<unsigned, unsigned>>{{0, 1}});

    // routing ignored bi directional edge and solves correctly
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.first.n_gates() == 2);
    CHECK(respects_connectivity_constraints(result.first, arc, false));
    REQUIRE(result.second == true);
  }
}

SCENARIO(
    "Test routing on a directed architecture doesn't throw an error if "
    "non-cx optype is presented",
    "[routing]") {
  GIVEN(
      "A simple two-qubit circuit with non-cx multi-qubit gates and a "
      "directed architecture") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CU1, 0.5, {1, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {0, 1});
    circ.add_op<unsigned>(OpType::CY, {1, 0});
    circ.add_op<unsigned>(OpType::CY, {0, 1});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::CRz, 0.5, {1, 0});
    circ.add_op<unsigned>(OpType::CRz, 0.5, {0, 1});

    Architecture arc(std::vector<std::pair<unsigned, unsigned>>{{0, 1}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.second == true);
    REQUIRE(result.first.n_gates() == 8);
  }
}

SCENARIO("Dense CX circuits route succesfully", "[routing]") {
  GIVEN(
      "Complex CX circuits for large directed architecture based off "
      "IBMTokyo") {
    Circuit circ(20);
    for (unsigned x = 0; x < 17; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {  // swap the way directed chain runs each time
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    Architecture arc(
        {{0, 1},   {1, 2},   {2, 3},   {3, 4},   {0, 5},   {1, 6},   {1, 7},
         {2, 6},   {2, 7},   {3, 8},   {3, 9},   {4, 8},   {4, 9},   {5, 6},
         {5, 10},  {5, 11},  {6, 10},  {6, 11},  {6, 7},   {7, 12},  {7, 13},
         {7, 8},   {8, 12},  {8, 13},  {8, 9},   {10, 11}, {11, 16}, {11, 17},
         {11, 12}, {12, 16}, {12, 17}, {12, 13}, {13, 18}, {13, 19}, {13, 14},
         {14, 18}, {14, 19}, {15, 16}, {16, 17}, {17, 18}, {18, 19}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.second);
    (Transform::decompose_SWAP_to_CX() >> Transform::decompose_BRIDGE_to_CX())
        .apply(result.first);
    Transform::decompose_CX_directed(arc).apply(result.first);
    REQUIRE(respects_connectivity_constraints(result.first, arc, true));
  }
}

SCENARIO(
    "Dense CX circuits route succesfully on undirected Ring with "
    "placement.",
    "[routing]") {
  GIVEN("Complex CX circuits, big ring") {
    Circuit circ(29);
    for (unsigned x = 0; x < 29; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    RingArch arc(29);
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.second);
    Transform::decompose_SWAP_to_CX().apply(result.first);
    REQUIRE(respects_connectivity_constraints(result.first, arc, false, true));
  }
}

SCENARIO(
    "Dense CX circuits route succesfully on smart placement unfriendly "
    "architecture.",
    "[routing]") {
  GIVEN("Complex CX circuits, big ring") {
    Circuit circ(13);
    for (unsigned x = 0; x < 13; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    Architecture arc(
        {{0, 1},
         {2, 0},
         {2, 4},
         {6, 4},
         {8, 6},
         {8, 10},
         {12, 10},
         {3, 1},
         {3, 5},
         {7, 5},
         {7, 9},
         {11, 9},
         {11, 13},
         {12, 13},
         {6, 7}});
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.second);
    REQUIRE(respects_connectivity_constraints(result.first, arc, false, true));
  }
}

SCENARIO("Empty circuits, with and without blank wires", "[routing]") {
  GIVEN("An empty circuit with some qubits") {
    Circuit circ(6);
    RingArch arc(6);
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.first.depth() == 0);
    REQUIRE(result.first.n_gates() == 0);
    REQUIRE(result.first.n_qubits() == 6);
    REQUIRE(result.second == true);
    REQUIRE(respects_connectivity_constraints(result.first, arc, true));
  }
  GIVEN("An empty circuit with no qubits") {
    Circuit circ(0);
    RingArch arc(6);
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();
    REQUIRE(result.second == false);
    REQUIRE(result.first.depth() == 0);
    REQUIRE(result.first.n_gates() == 0);
    REQUIRE(result.first.n_qubits() == 0);
  }

  GIVEN("An empty circuit with no qubits, and empty architecture") {
    Circuit circ(0);
    std::vector<std::pair<Node, Node>> cons = {};
    Architecture arc(cons);
    REQUIRE_THROWS_AS(
        [&]() { Routing router(circ, arc); }(), ArchitectureMismatch);
  }
  GIVEN("An a mismatch") {
    Circuit circ(5);
    RingArch arc(4);
    REQUIRE_THROWS_AS(
        [&]() { Routing router(circ, arc); }(), ArchitectureMismatch);
  }
}

/* METHODS TO COVER IN TESTING: */
/* Routing class: */

// Routing::increment_distance
SCENARIO("Does increment distance work?", "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {2, 3}, {4, 5}});
  SquareGrid test_architecture(2, 3);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(2, 3);
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  GIVEN("Suitable Distance vector, Swap and increment.") {
    unsigned diameter = test_architecture.get_diameter();
    graphs::dist_vec test_distance(diameter, 2);
    Swap test_swap = {square_nodes[0], square_nodes[1]};
    int increment = 2;
    unsigned distance_index = diameter - test_architecture.get_distance(
                                             test_swap.first, test_swap.second);
    int pre_increment_val = test_distance[distance_index];
    routing_tester.increment_distance(test_distance, test_swap, increment);
    REQUIRE(pre_increment_val + increment == test_distance[distance_index]);
  }
  GIVEN("Realistic Distance Vector, non_adjacent Swap, absurd increment.") {
    unsigned diameter = test_architecture.get_diameter();
    graphs::dist_vec test_distance(diameter, 2);
    Swap test_swap = {square_nodes[0], square_nodes[5]};
    int increment = 30;
    unsigned distance_index = diameter - test_architecture.get_distance(
                                             test_swap.first, test_swap.second);
    int pre_increment_val = test_distance[distance_index];
    routing_tester.increment_distance(test_distance, test_swap, increment);
    REQUIRE(pre_increment_val + increment == test_distance[distance_index]);
  }
}

// Routing::generate_distance_vector
SCENARIO("Does generate_distance_vector work suitably?", "[routing]") {
  GIVEN("A realistic small interaction vector and architecture") {
    // Creating RoutingTester object
    Circuit test_circuit(6);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {2, 3}, {4, 5}});
    SquareGrid test_architecture(3, 2);
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
    // 0 -- 1
    // |    |
    // 2 -- 3
    // |    |
    // 4 -- 5
    Routing test_router(test_circuit, test_architecture);
    RoutingTester routing_tester(&test_router);

    std::array<unsigned, 6> inte_pattern = {1, 0, 5, 3, 4, 2};
    Interactions test_interaction;
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    // no placement invoked, should be 0 at diameter distance, 1 at distance 2,
    // 1 at distance 1. i.e. {0,2}
    graphs::dist_vec out_distances =
        routing_tester.generate_distance_vector(test_interaction);
    REQUIRE(
        out_distances[0] ==
        0);  // 0 entries at distance diameter away from each other
    REQUIRE(
        out_distances[1] ==
        2);  // 2 entries at distance diameter - 1  away from each other
  }
  GIVEN("A realistic large interaction vector and architecture") {
    // Creating larger RoutingTester object
    Circuit test_circuit(10);
    SquareGrid test_architecture(2, 5);
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(2, 5);
    // 0 -- 1 -- 2 -- 3 -- 4
    // |    |    |    |    |
    // 5 -- 6 -- 7 -- 8 -- 9
    Routing test_router(test_circuit, test_architecture);
    RoutingTester routing_tester(&test_router);

    unsigned ind = 0;
    Interactions test_interaction;
    std::array<unsigned, 10> inte_pattern{9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    // Expected distances:
    // 9-0 -> 5
    // 8-1 -> 3
    // 7-2 -> 1
    // 6-3 -> 3
    // 5-4 -> 5
    // i.e.
    graphs::dist_vec expected_distances = {
        4, 0, 4, 0};  // 4 qubits at diameter, 0 at diameter-1, 4 qubits at
                      // diameter-2, 0 at diameter-3
    graphs::dist_vec out_distances =
        routing_tester.generate_distance_vector(test_interaction);
    REQUIRE(out_distances == expected_distances);
  }
}

// Routing::update_distance_vector
SCENARIO("Does update_distance_vector update as intended?", "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  SquareGrid test_architecture(3, 2);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  // update_distance_vector is four indiviudal increment_distances
  GIVEN("Realistic Distance vector, Swap and Interaction vector.") {
    unsigned diameter = test_architecture.get_diameter();
    graphs::dist_vec test_distance = {0, 2};
    unsigned ind = 0;
    Interactions test_interaction;
    std::array<unsigned, 6> inte_pattern{1, 0, 5, 3, 4, 2};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    graphs::dist_vec quick_compare_distance =
        routing_tester.generate_distance_vector(test_interaction);
    REQUIRE(quick_compare_distance == test_distance);

    Swap test_swap = {square_nodes[2], square_nodes[4]};

    // Distances from full method
    graphs::dist_vec out_distance = routing_tester.update_distance_vector(
        test_swap, test_distance, test_interaction);
    // Forming Distances from individual steps:
    // (1) 2 in test_swap is interacting with qubit 5, a distance of 2 away
    // this swap brings the two qubits adjacent
    unsigned distance_index_1 =
        diameter - test_architecture.get_distance(
                       test_swap.first, test_interaction[test_swap.first]);
    int pre_increment_val_1 = test_distance[distance_index_1];
    routing_tester.increment_distance(
        test_distance, {test_swap.first, test_interaction[test_swap.first]},
        -2);
    REQUIRE(pre_increment_val_1 - 2 == test_distance[distance_index_1]);
    // (2), 4 in test_swap is not interacting, test_distances won't change
    REQUIRE(
        test_architecture.get_distance(
            test_swap.second, test_interaction[test_swap.second]) == 0);
    // (3), 4 in test_swap and the qubit 2 is interacting with 5 are adjacent,
    // test_distances won't change
    REQUIRE(
        test_architecture.get_distance(
            test_swap.second, test_interaction[test_swap.first]) == 1);
    // (4), 2 in test swap and the qubit 4 is interacting with 0 are adjacent,
    // test_distances won't change
    REQUIRE(
        test_architecture.get_distance(
            test_swap.first, test_interaction[test_swap.second]) == 1);

    REQUIRE(out_distance[0] == test_distance[0]);
    REQUIRE(out_distance[1] == test_distance[1]);
  }
}
// Routing::pair_dists
SCENARIO(
    "Does pair_dists return the correct distances, in the correct order?",
    "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  SquareGrid test_architecture(3, 2);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  GIVEN(
      "Realistic architecture nodes. Distance between pair_1 less than "
      "between pair_2.") {
    std::pair<Node, Node> pair_1 = {square_nodes[0], square_nodes[3]};
    std::pair<Node, Node> pair_2 = {square_nodes[1], square_nodes[4]};
    unsigned dist_1 =
        test_architecture.get_distance(pair_1.first, pair_1.second);
    REQUIRE(dist_1 == 2);
    unsigned dist_2 =
        test_architecture.get_distance(pair_2.first, pair_2.second);
    REQUIRE(dist_2 == 3);
    std::pair<unsigned, unsigned> pair_dists_results =
        routing_tester.pair_dists(
            pair_1.first, pair_1.second, pair_2.first, pair_2.second);
    REQUIRE(pair_dists_results.first == dist_2);
    REQUIRE(pair_dists_results.second == dist_1);
  }
  GIVEN(
      "Realistic architecture nodes. Distance between pair_1 greater than "
      "between pair_2.") {
    std::pair<Node, Node> pair_1 = {square_nodes[4], square_nodes[3]};
    std::pair<Node, Node> pair_2 = {square_nodes[0], square_nodes[2]};
    unsigned dist_1 =
        test_architecture.get_distance(pair_1.first, pair_1.second);
    REQUIRE(dist_1 == 2);
    unsigned dist_2 =
        test_architecture.get_distance(pair_2.first, pair_2.second);
    REQUIRE(dist_2 == 1);
    std::pair<unsigned, unsigned> pair_dists_results =
        routing_tester.pair_dists(
            pair_1.first, pair_1.second, pair_2.first, pair_2.second);
    REQUIRE(pair_dists_results.first == dist_1);
    REQUIRE(pair_dists_results.second == dist_2);
  }
}

// Routing::swap_decreases
SCENARIO(
    "Does swap_decreases correctly determine between two placements?",
    "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  SquareGrid test_architecture(3, 2);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  GIVEN("A swap that improves placement for given interaction vector.") {
    // only nodes 0 and 5 have an interacting pair of qubits between them
    Interactions test_interaction;
    unsigned ind = 0;
    std::array<unsigned, 6> inte_pattern{5, 1, 2, 3, 4, 0};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    Swap test_swap = {square_nodes[0], square_nodes[2]};
    Swap test_swap_interactions = {square_nodes[5], square_nodes[2]};
    // Confirm swap_decreases functions as expected
    REQUIRE(routing_tester.swap_decreases(test_swap, test_interaction) == true);
    // Confirm working elements of swap_decreases does also
    unsigned dist_1 = test_architecture.get_distance(
        test_swap.first, test_swap_interactions.first);
    REQUIRE(dist_1 == 3);
    unsigned dist_2 = test_architecture.get_distance(
        test_swap.second, test_swap_interactions.second);
    REQUIRE(dist_2 == 0);
    unsigned dist_3 = test_architecture.get_distance(
        test_swap.second, test_swap_interactions.first);
    REQUIRE(dist_3 == 2);
    unsigned dist_4 = test_architecture.get_distance(
        test_swap.first, test_swap_interactions.second);
    REQUIRE(dist_4 == 1);

    std::pair<unsigned, unsigned> old_dists = {dist_1, dist_2};
    std::pair<unsigned, unsigned> new_dists = {dist_3, dist_4};
    REQUIRE(new_dists < old_dists);
  }
  GIVEN("A swap containing non-interacting nodes.") {
    unsigned ind = 0;
    Interactions test_interaction;
    // only nodes 0 and 5 have an interacting pair of qubits between them
    std::array<unsigned, 6> inte_pattern{5, 1, 2, 3, 4, 0};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    Swap test_swap = {square_nodes[1], square_nodes[3]};
    Swap test_swap_interactions = {square_nodes[1], square_nodes[3]};
    // Confirm swap_decreases functions as expected
    REQUIRE(
        routing_tester.swap_decreases(test_swap, test_interaction) == false);
  }
}

// Routing::candidate_swaps
SCENARIO("Does candidate swaps return all suitable edges?", "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  SquareGrid test_architecture(3, 2);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5
  std::vector<Connection> test_arc = test_architecture.get_connections_vec();
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  GIVEN("One pair of interacting qubits, four suitable edges between them.") {
    unsigned ind = 0;
    Interactions test_interaction;
    std::array<unsigned, 6> inte_pattern{3, 1, 2, 0, 4, 5};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    std::vector<Swap> correct_swaps = {
        {square_nodes[0], square_nodes[1]},
        {square_nodes[0], square_nodes[2]},
        {square_nodes[1], square_nodes[3]},
        {square_nodes[2], square_nodes[3]}};
    std::vector<Swap> test_swaps =
        routing_tester.candidate_swaps(test_arc, test_interaction);
    REQUIRE(test_swaps.size() == 4);
    REQUIRE(test_swaps == correct_swaps);
  }
  GIVEN("A case wherein no edges are suitable.") {
    unsigned ind = 0;
    Interactions test_interaction;
    // easiest to replicate this case by making all interactions adjacent
    std::array<unsigned, 6> inte_pattern{1, 0, 3, 2, 5, 4};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    std::vector<Swap> test_swaps =
        routing_tester.candidate_swaps(test_arc, test_interaction);
    REQUIRE(test_swaps.size() == 0);
  }
  GIVEN(
      "A case with all qubits interacting, 5 suitable edges between "
      "them.") {
    unsigned ind = 0;
    Interactions test_interaction;
    std::array<unsigned, 6> inte_pattern{5, 2, 1, 4, 3, 0};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    std::vector<Swap> correct_swaps = {
        {square_nodes[0], square_nodes[1]},
        {square_nodes[0], square_nodes[2]},
        {square_nodes[2], square_nodes[3]},
        {square_nodes[3], square_nodes[5]},
        {square_nodes[4], square_nodes[5]}};
    std::vector<Swap> test_swaps =
        routing_tester.candidate_swaps(test_arc, test_interaction);
    REQUIRE(test_swaps.size() == 5);
    REQUIRE(test_swaps == correct_swaps);
  }
}

// Routing::cowtan_et_al_heuristic
SCENARIO(
    "Does implementation of heuristic outlined in paper work as expected?",
    "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  SquareGrid test_architecture(3, 2);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  GIVEN("One pair of interacting qubits, four suitable swap gates.") {
    std::vector<Swap> test_swaps = {
        {square_nodes[0], square_nodes[1]},
        {square_nodes[0], square_nodes[2]},
        {square_nodes[1], square_nodes[3]},
        {square_nodes[2], square_nodes[3]},
        {square_nodes[3], square_nodes[5]}};
    graphs::dist_vec test_distances = {0, 2};
    Interactions test_interaction;
    unsigned ind = 0;
    std::array<unsigned, 6> inte_pattern{3, 1, 2, 0, 4, 5};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    std::vector<Swap> output_swaps = routing_tester.cowtan_et_al_heuristic(
        test_swaps, test_distances, test_interaction);
    std::vector<Swap> expected_output = {
        {square_nodes[0], square_nodes[1]},
        {square_nodes[0], square_nodes[2]},
        {square_nodes[1], square_nodes[3]},
        {square_nodes[2], square_nodes[3]}};
    REQUIRE(output_swaps == expected_output);
  }
  GIVEN("Two pairs of interacting qubits, two suitable swap gates.") {
    std::vector<Swap> test_swaps = {
        {square_nodes[0], square_nodes[1]}, {square_nodes[0], square_nodes[2]},
        {square_nodes[1], square_nodes[3]}, {square_nodes[2], square_nodes[3]},
        {square_nodes[2], square_nodes[4]}, {square_nodes[3], square_nodes[5]},
        {square_nodes[4], square_nodes[5]}};
    unsigned ind = 0;
    Interactions test_interaction;
    std::array<unsigned, 6> inte_pattern{3, 4, 2, 0, 1, 5};
    for (unsigned i = 0; i < inte_pattern.size(); ++i) {
      test_interaction.insert({square_nodes[i], square_nodes[inte_pattern[i]]});
    }
    graphs::dist_vec test_distances = {2, 2};
    std::vector<Swap> output_swaps = routing_tester.cowtan_et_al_heuristic(
        test_swaps, test_distances, test_interaction);
    std::vector<Swap> expected_output = {
        {square_nodes[0], square_nodes[1]}, {square_nodes[1], square_nodes[3]}};
    REQUIRE(output_swaps == expected_output);
  }
}

// Routing::update_qmap
SCENARIO("Does update qmap correctly update mapping from swap?", "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(2);
  RingArch test_architecture(2);
  node_vector_t ring_nodes = RingArch::get_nodes_canonical_order(2);
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);
  Qubit qb0(0);
  Qubit qb1(1);
  qubit_bimap_t test_map;
  test_map.left.insert({qb0, ring_nodes[0]});
  test_map.left.insert({qb1, ring_nodes[1]});

  routing_tester.update_qmap(test_map, {ring_nodes[0], ring_nodes[1]});
  REQUIRE(test_map.right.at(ring_nodes[0]) == qb1);
  REQUIRE(test_map.right.at(ring_nodes[1]) == qb0);
}

// Routing::solve_furthest interior functions
SCENARIO(
    "Do solve_furthest interior methods find and swap along the expected "
    "path?",
    "[routing]") {
  // Creating RoutingTester object
  Circuit test_circuit(6);
  add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {2, 3}, {4, 5}});
  SquareGrid test_architecture(3, 2);
  node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5
  Routing test_router(test_circuit, test_architecture);
  RoutingTester routing_tester(&test_router);

  unsigned node0, node1;
  node_vector_t expected_path;
  std::vector<Swap> expected_swaps;
  GIVEN("An expected path with an even number of nodes.") {
    node0 = 0, node1 = 5;
    expected_path = {
        square_nodes[5], square_nodes[3], square_nodes[1], square_nodes[0]};
    expected_swaps = {
        {square_nodes[5], square_nodes[3]}, {square_nodes[3], square_nodes[1]}};
  }
  GIVEN("An expected path with an odd number of nodes.") {
    node0 = 0, node1 = 3;
    expected_path = {square_nodes[3], square_nodes[1], square_nodes[0]};
    expected_swaps = {{square_nodes[3], square_nodes[1]}};
  }
  GIVEN("An adjacent path doesn't fail.") {
    node0 = 0, node1 = 1;
    expected_path = {square_nodes[1], square_nodes[0]};
    expected_swaps = {};
  }
  // Collect path from architecture
  const node_vector_t test_path =
      test_architecture.get_path(square_nodes[node0], square_nodes[node1]);
  REQUIRE(test_path == expected_path);
  qubit_bimap_t test_map;
  for (unsigned i = 0; i < 6; i++) {
    test_map.left.insert({Qubit(i), square_nodes[i]});
  }
  routing_tester.set_qmap(test_map);
  const std::vector<Swap> path_swaps = routing_tester.path_to_swaps(test_path);
  REQUIRE(path_swaps == expected_swaps);
}

// generate_test_interaction_graph and qubit_lines
SCENARIO("Test interaction graph and line generation", "[routing]") {
  // 0 -- 1
  // |    |
  // 2 -- 3
  // |    |
  // 4 -- 5

  GIVEN("A small test circuit with 1 layer, all qubits in 2qb gates.") {
    Circuit test_circuit(6);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {2, 3}, {4, 5}});
    QubitGraph test_qubit_graph = generate_interaction_graph(test_circuit);

    REQUIRE(test_qubit_graph.n_connections() == 3);
    REQUIRE(test_qubit_graph.connection_exists(Qubit(0), Qubit(1)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(2), Qubit(3)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(4), Qubit(5)));

    QubitLineList qlines = qubit_lines(test_circuit);
    QubitLineList correct_lines = {
        {Qubit(0), Qubit(1)}, {Qubit(2), Qubit(3)}, {Qubit(4), Qubit(5)}};
    REQUIRE(qlines == correct_lines);
  }
  GIVEN("A small test circuit with 1 layer, not all qubits in 2qb gates.") {
    Circuit test_circuit(6);
    test_circuit.add_op<unsigned>(OpType::CX, {0, 1});
    test_circuit.add_op<unsigned>(OpType::H, {5});
    test_circuit.add_op<unsigned>(OpType::H, {3});
    test_circuit.add_op<unsigned>(OpType::CX, {2, 4});

    QubitGraph test_qubit_graph = generate_interaction_graph(test_circuit);

    REQUIRE(test_qubit_graph.n_connections() == 2);
    REQUIRE(test_qubit_graph.connection_exists(Qubit(0), Qubit(1)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(2), Qubit(4)));

    QubitLineList qlines = qubit_lines(test_circuit);
    QubitLineList correct_lines = {
        {Qubit(0), Qubit(1)},
        {Qubit(2), Qubit(4)},
        {Qubit(3)},
        {Qubit(5)}};  // It is not guaranteed to match qubit numbers as qubits
                      // are not unsigneds
    REQUIRE(qlines == correct_lines);
  }
  GIVEN("A small test circuit with 2 layers.") {
    Circuit test_circuit(6);
    add_2qb_gates(
        test_circuit, OpType::CX,
        {{0, 1}, {2, 3}, {4, 5}, {2, 1}, {4, 3}, {5, 1}});

    QubitGraph test_qubit_graph = generate_interaction_graph(test_circuit);

    REQUIRE(test_qubit_graph.n_connections() == 5);
    REQUIRE(test_qubit_graph.connection_exists(Qubit(0), Qubit(1)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(2), Qubit(3)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(4), Qubit(5)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(2), Qubit(1)));
    REQUIRE(test_qubit_graph.connection_exists(Qubit(4), Qubit(3)));

    QubitLineList qlines = qubit_lines(test_circuit);
    QubitLineList correct_lines = {
        {Qubit(0), Qubit(1), Qubit(2), Qubit(3), Qubit(4), Qubit(5)}};
    REQUIRE(qlines == correct_lines);
  }
}

// solve_with_map
SCENARIO("Test routing with partial map provided", "[routing]") {
  GIVEN("A partial map where no node should be removed.") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    Architecture arc({{0, 1}, {1, 2}});

    // force a partial map which requires unused node to solve
    Placement pl(arc);
    qubit_mapping_t map_ = {{Qubit(0), Node(0)}, {Qubit(2), Node(2)}};
    pl.place_with_map(circ, map_);
    Routing router(circ, arc);
    std::pair<Circuit, bool> result = router.solve();

    REQUIRE(result.second == true);
    // check solution is valid and respects map
    std::vector<Command> test_coms = result.first.get_commands();
    REQUIRE(test_coms.size() == 3);
    bool oph = (*test_coms[0].get_op_ptr() == *get_op_ptr(OpType::H));
    oph &= (test_coms[0].get_args()[0] == Node(0));
    REQUIRE(oph);
    REQUIRE(*test_coms[1].get_op_ptr() == *get_op_ptr(OpType::SWAP));
    REQUIRE(*test_coms[2].get_op_ptr() == *get_op_ptr(OpType::CX));
  }

  GIVEN("A mapped set of nodes") {
    Circuit circ(4);
    Qubit qb0(0);
    Qubit qb1(1);
    Qubit qb2(2);
    Qubit qb3(3);
    // test removal only happpens if subgraph remains connected
    SquareGrid test_architecture(3, 2);
    Architecture subarc = test_architecture;
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 2);
    // 0 -- 1
    // |    |
    // 2 -- 3
    // |    |
    // 4 -- 5
    // subarc = {0, 1, 3}
    subarc.remove_uid(square_nodes[5]);
    subarc.remove_uid(square_nodes[4]);
    subarc.remove_uid(square_nodes[3]);
    REQUIRE(subgraph_remove_if_connected(
        test_architecture, subarc, square_nodes[3]));
    REQUIRE(!subgraph_remove_if_connected(
        test_architecture, subarc, square_nodes[1]));
    REQUIRE(subgraph_remove_if_connected(
        test_architecture, subarc, square_nodes[4]));
    REQUIRE(subgraph_remove_if_connected(
        test_architecture, subarc, square_nodes[5]));
    REQUIRE(test_architecture.n_connections() == 2);

    SquareGrid test_architecture2(3, 2);
    qubit_bimap_t map;
    map.left.insert({qb0, square_nodes[0]});
    map.left.insert({qb1, square_nodes[1]});
    map.left.insert({qb2, square_nodes[2]});

    remove_unmapped_nodes(test_architecture2, map, circ);
    REQUIRE(test_architecture2.n_connections() == 2);
    REQUIRE(
        test_architecture2.connection_exists(square_nodes[0], square_nodes[1]));
    REQUIRE(
        test_architecture2.connection_exists(square_nodes[0], square_nodes[2]));

    // test unmapped nodes which cannot be removed are mapped to a qubit
    map.left.erase(qb0);

    remove_unmapped_nodes(test_architecture2, map, circ);
    REQUIRE(map.left.find(qb0)->second == square_nodes[0]);
    REQUIRE(test_architecture2.n_connections() == 2);
    REQUIRE(
        test_architecture2.connection_exists(square_nodes[0], square_nodes[1]));
    REQUIRE(
        test_architecture2.connection_exists(square_nodes[0], square_nodes[2]));

    // test when an unmapped node is mapped, the most connected is chosen
    // (i.e. least connected nodes are removed first)
    Architecture test_architecture3({{0, 1}, {0, 2}, {1, 3}, {2, 3}, {2, 4}});
    // 0 -- 1
    // |    |
    // 2 -- 3
    // |
    // 4

    qubit_bimap_t map2;
    map2.left.insert({qb0, Node(0)});
    map2.left.insert({qb3, Node(3)});
    remove_unmapped_nodes(test_architecture3, map2, circ);
    REQUIRE(map2.right.find(Node(2))->second == qb1);
    bool no_4 = map2.right.find(Node(4)) == map2.right.end();
    REQUIRE(no_4);
  }
}

// Every command in the circuit with a specified optype
// must have a specified single qubit argument.
static void require_arguments_for_specified_commands(
    const Circuit &circ, const std::map<OpType, Qubit> &the_map) {
  for (Command com : circ) {
    const auto type = com.get_op_ptr()->get_type();
    const auto citer = the_map.find(type);
    if (citer != the_map.cend()) {
      unit_vector_t comp = {citer->second};
      REQUIRE(com.get_args() == comp);
    }
  }
}

SCENARIO(
    "Does shifting single qubit gates through SWAP gates to get find nodes "
    "with better fidelity work?",
    "[routing]") {
  Architecture arc({{0, 1}, {1, 2}});
  gate_error_t ge_0(0.3);
  gate_error_t ge_1(0.2);
  gate_error_t ge_2(0.1);
  op_node_errors_t nec;
  op_errors_t gec_0({{OpType::H, ge_0}, {OpType::X, ge_1}});
  op_errors_t gec_1({{OpType::H, ge_1}, {OpType::X, ge_2}});
  op_errors_t gec_2({{OpType::H, ge_2}, {OpType::X, ge_0}});

  nec.insert({Node(2), gec_2});
  nec.insert({Node(0), gec_0});
  nec.insert({Node(1), gec_1});

  GIVEN(
      "A simple two qubit circuit with clear difference between node "
      "fidelities") {
    Circuit circ(2);
    add_1qb_gates(circ, OpType::H, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    reassign_boundary(circ);
    Transform::commute_SQ_gates_through_SWAPS(nec).apply(circ);
    require_arguments_for_specified_commands(
        circ, {{OpType::H, circ.all_qubits().at(1)}});
  }
  GIVEN(
      "A simple two qubit circuit with multiple single qubit operations "
      "requiring movement before a SWAP.") {
    Circuit circ(2);
    add_1qb_gates(circ, OpType::H, {0, 0, 0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    reassign_boundary(circ);
    Transform::commute_SQ_gates_through_SWAPS(nec).apply(circ);
    require_arguments_for_specified_commands(
        circ, {{OpType::H, circ.all_qubits().at(1)}});
  }
  GIVEN("Multiple SWAP gates, multiple single qubit gates.") {
    Circuit circ(3);
    add_1qb_gates(circ, OpType::H, {0, 0, 0, 1});
    add_2qb_gates(circ, OpType::SWAP, {{0, 1}, {1, 2}, {0, 1}, {1, 2}, {1, 2}});

    reassign_boundary(circ);
    Transform::commute_SQ_gates_through_SWAPS(nec).apply(circ);
    require_arguments_for_specified_commands(
        circ, {{OpType::H, circ.all_qubits().at(2)}});
  }
  GIVEN(
      "Multiple SWAP gates, multiple single qubit gates of various "
      "OpType.") {
    Circuit circ(3);
    add_1qb_gates(circ, OpType::X, {0, 0, 1, 1});
    add_1qb_gates(circ, OpType::H, {0, 0, 0, 1});
    add_2qb_gates(circ, OpType::SWAP, {{0, 1}, {1, 2}, {0, 1}, {1, 2}, {1, 2}});

    reassign_boundary(circ);
    Transform::commute_SQ_gates_through_SWAPS(nec).apply(circ);
    const qubit_vector_t qbs = circ.all_qubits();
    require_arguments_for_specified_commands(
        circ, {{OpType::H, qbs.at(2)}, {OpType::X, qbs.at(1)}});
  }
  GIVEN(
      "A large circuit of CX gates, H gates and X gates, routed and "
      "shifted.") {
    Circuit circ(9);
    for (unsigned x = 0; x < circ.n_qubits(); ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {
          circ.add_op<unsigned>(OpType::SWAP, {x, y});
          circ.add_op<unsigned>(OpType::X, {x});
          circ.add_op<unsigned>(OpType::H, {x});
          circ.add_op<unsigned>(OpType::SWAP, {y + 1, y});
        } else {
          circ.add_op<unsigned>(OpType::SWAP, {y, x});
          circ.add_op<unsigned>(OpType::H, {y});
          circ.add_op<unsigned>(OpType::X, {y});
          circ.add_op<unsigned>(OpType::SWAP, {y, y + 1});
        }
      }
    }
    SquareGrid arc(3, 3);
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 3);

    const std::vector<gate_error_t> gate_errors{
        0.3, 0.2, 0.1, 0.02, 0.22, 0.46, 0.18, 1.0 - 0.907, 1.0 - 0.7241};
    REQUIRE(arc.get_columns() * arc.get_rows() == gate_errors.size());
    REQUIRE(gate_errors.size() == square_nodes.size());
    REQUIRE(circ.n_qubits() == gate_errors.size());
    op_node_errors_t nec;
    unsigned ind = 0;
    for (unsigned nn = 0; nn < square_nodes.size(); ++nn) {
      nec[square_nodes[nn]] = op_errors_t(
          {{OpType::H, gate_errors[nn]},
           {OpType::X, gate_errors[(nn + 3) % gate_errors.size()]}});
    }
    DeviceCharacterisation characterisation(nec);

    Circuit test_0 = circ;
    reassign_boundary(test_0, square_nodes);
    Transform::decompose_SWAP_to_CX().apply(test_0);
    const auto sv0 = tket_sim::get_statevector(test_0);
    double pre_aggregate = 0;

    qubit_bimap_t qmap;
    qubit_vector_t free_qs = test_0.all_qubits();
    for (unsigned u = 0; u < free_qs.size(); u++) {
      qmap.insert({free_qs[u], square_nodes[u]});
    }

    for (Command com : test_0) {
      OpType ot = com.get_op_ptr()->get_type();
      if (ot == OpType::X || ot == OpType::H) {
        Node n = qmap.left.at(Qubit(com.get_args()[0]));
        pre_aggregate += 1.0 - characterisation.get_error(Node(n), ot);
      }
    }
    reassign_boundary(circ, square_nodes);
    Transform::commute_SQ_gates_through_SWAPS(nec).apply(circ);
    Circuit test_1 = circ;
    Transform::decompose_SWAP_to_CX().apply(test_1);
    const auto sv1 = tket_sim::get_statevector(test_1);
    double post_aggregate = 0;
    for (Command com : test_1) {
      OpType ot = com.get_op_ptr()->get_type();
      if (ot == OpType::X || ot == OpType::H) {
        Node n = qmap.left.at(Qubit(com.get_args()[0]));
        post_aggregate += 1.0 - characterisation.get_error(Node(n), ot);
      }
    }
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(sv0, sv1));
    REQUIRE(post_aggregate > pre_aggregate);
  }
}
SCENARIO("Test barrier is ignored by routing") {
  GIVEN("Circuit with 1qb barrier") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_barrier(uvec{0});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    SquareGrid test_architecture(1, 3);
    GraphPlacement gp(test_architecture);
    gp.place(circ);
    Routing router(circ, test_architecture);
    Circuit pc = router.solve().first;
    REQUIRE(pc.depth() == 2);
    check_command_types(
        pc, {OpType::CX, OpType::Rz, OpType::CX, OpType::Barrier});
  }
  GIVEN("Circuit with 2 qb barrier") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_barrier({0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    SquareGrid test_architecture(1, 2);
    Routing router(circ, test_architecture);
    check_command_types(
        router.solve().first, {OpType::CX, OpType::Barrier, OpType::CX});
  }
  GIVEN("Circuit with 4 qb barrier, using gen_full_mapping_pass.") {
    const std::vector<Node> nums = {
        Node("rig", 21), Node("rig", 22), Node("rig", 25), Node("rig", 35),
        Node("rig", 36)};
    const std::vector<std::pair<unsigned, unsigned>> coupling_list_indices = {
        {0, 1}, {0, 4}, {1, 0}, {1, 3}, {4, 0}, {4, 3}, {3, 1}, {3, 4}};

    std::vector<std::pair<Node, Node>> coupling_list;
    for (const auto &pair : coupling_list_indices) {
      coupling_list.push_back(
          std::make_pair(nums[pair.first], nums[pair.second]));
    }
    Architecture arc(coupling_list);
    Circuit circ(4);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    add_2qb_gates(circ, OpType::CZ, {{1, 2}, {3, 2}, {3, 1}});
    circ.add_barrier({0, 1, 2, 3});

    std::vector<std::reference_wrapper<RoutingMethod>> config = {
        LexiRouteRoutingMethod(100)};
    PlacementPtr pp = std::make_shared<GraphPlacement>(arc);
    PassPtr p = gen_full_mapping_pass(arc, pp, config);
    CompilationUnit cu(circ);
    p->apply(cu);
    REQUIRE(
        respects_connectivity_constraints(cu.get_circ_ref(), arc, false, true));
  }
  GIVEN("Check Circuit with 2qb barrier does not add swaps for the barrier") {
    Circuit circ(3);
    Architecture line({{0, 1}, {1, 2}});
    GraphPlacement gp(line);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    circ.add_barrier({0, 2});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    gp.place(circ);

    Routing router(circ, line);
    qubit_vector_t all_qs_pre = circ.all_qubits();
    std::pair<Circuit, bool> pc = router.solve();
    qubit_vector_t all_qs_post = circ.all_qubits();
    REQUIRE(all_qs_pre == all_qs_post);
    REQUIRE(pc.second == false);
    REQUIRE(pc.first.depth() == 4);
  }
  GIVEN(
      "Check Circuit with 2qb barrier does not add swaps for the barrier, "
      "with no placement.") {
    Circuit circ(3);
    Architecture line({{0, 1}, {1, 2}});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    circ.add_barrier({0, 2});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    GraphPlacement gp(line);
    gp.place(circ);
    Routing router(circ, line);
    unsigned pre_depth = circ.depth();
    std::pair<Circuit, bool> pc = router.solve({});
    REQUIRE(pc.second == false);
    unsigned post_depth = pc.first.depth();
    REQUIRE(post_depth == pre_depth);
    REQUIRE(post_depth == 4);
  }
  GIVEN("Circuit with 3 qb barrier") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    circ.add_barrier({0, 1, 2});
    Architecture line({{0, 1}, {1, 2}});
    GraphPlacement gp(line);
    gp.place(circ);
    Routing router(circ, line);
    qubit_vector_t all_qs_pre = circ.all_qubits();
    std::pair<Circuit, bool> pc = router.solve();
    qubit_vector_t all_qs_post = circ.all_qubits();
    REQUIRE(all_qs_pre == all_qs_post);
    REQUIRE(pc.first.depth() == 2);
    REQUIRE(pc.second == false);
  }
}

SCENARIO(
    "Does identification and insertion of bridge circuit do as expected?") {
  GIVEN(
      "A proposed SWAP which will act detrimentally for the next timestep, "
      "i.e. a bridge should be inserted.") {
    Circuit circ(9);
    for (unsigned i = 0; i < 9; i++) {
      circ.add_op<unsigned>(OpType::H, {i});
    }
    add_2qb_gates(circ, OpType::CX, {{0, 4}, {3, 8}, {4, 7}, {3, 6}});

    SquareGrid arc(3, 3);
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 3);
    Routing router(circ, arc);
    RoutingTester test_router(&router);
    RoutingConfig new_config(50, 0, 0, 0);
    test_router.set_config(new_config);
    qubit_bimap_t qmap = test_router.set_default_initial_map(square_nodes);
    test_router.initialise_slicefrontier();
    test_router.set_interaction();
    std::pair<std::pair<bool, Node>, std::pair<bool, Node>> output =
        test_router.check_distributed_cx({square_nodes[1], square_nodes[4]});
    std::pair<std::pair<bool, Node>, std::pair<bool, Node>> expected = {
        {false, Node(0)}, {false, Node(0)}};
    REQUIRE(output == expected);
  }
  GIVEN(
      "A proposed SWAP which will act better for the next timestep, i.e. a "
      "bridge should not be inserted.") {
    Circuit circ(9);
    for (unsigned i = 0; i < 9; i++) {
      circ.add_op<unsigned>(OpType::H, {i});
    }
    add_2qb_gates(circ, OpType::CX, {{0, 4}, {3, 8}, {4, 7}, {3, 6}});

    SquareGrid arc(3, 3);
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(3, 3);
    Routing router(circ, arc);
    RoutingTester test_router(&router);
    qubit_bimap_t qmap = test_router.set_default_initial_map(square_nodes);
    RoutingConfig new_config(50, 0, 0, 0);
    test_router.set_config(new_config);
    test_router.initialise_slicefrontier();
    test_router.set_interaction();
    std::pair<std::pair<bool, Node>, std::pair<bool, Node>> output =
        test_router.check_distributed_cx({square_nodes[3], square_nodes[4]});
    std::pair<std::pair<bool, Node>, std::pair<bool, Node>> expected{
        {false, Node(0)}, {false, Node(0)}};
    REQUIRE(output == expected);
  }
  GIVEN("Multiple bridges to be inserted.") {
    Circuit circ(6);
    SquareGrid arc(6, 1);
    add_2qb_gates(circ, OpType::CX, {{0, 2}, {3, 5}, {1, 3}});
    node_vector_t square_nodes = SquareGrid::get_nodes_canonical_order(6, 1);
    Routing router(circ, arc);
    RoutingTester test_router(&router);
    qubit_bimap_t qmap = test_router.set_default_initial_map(square_nodes);
    test_router.initialise_slicefrontier();
    test_router.set_interaction();
    test_router.add_distributed_cx(
        square_nodes[3], square_nodes[5], square_nodes[4]);
    test_router.add_distributed_cx(
        square_nodes[0], square_nodes[2], square_nodes[1]);
    REQUIRE(test_router.get_circ()->n_gates() == 3);
    test_router.advance_frontier();
  }
  GIVEN("Consecutive CX edge case") {
    Circuit circ(5);
    Architecture arc({{0, 1}, {1, 2}, {0, 3}, {1, 4}, {3, 4}});
    add_2qb_gates(circ, OpType::CX, {{1, 2}, {0, 2}, {0, 1}});
    Routing router(circ, arc);
    RoutingTester test_router(&router);
    qubit_bimap_t qmap = test_router.set_default_initial_map();
    test_router.initialise_slicefrontier();
    test_router.advance_frontier();
    test_router.set_interaction();
    test_router.add_distributed_cx(Node(0), Node(2), Node(1));
    test_router.advance_frontier();
  }
}

SCENARIO(
    "Do Placement and Routing work if the given graph perfectly solves the "
    "problem?") {
  GIVEN("A perfect example without clifford_simp") {
    Circuit circ(5);
    add_2qb_gates(
        circ, OpType::CX,
        {{1, 2},
         {0, 3},
         {1, 4},
         {1, 2},
         {0, 1},
         {2, 0},
         {2, 1},
         {0, 1},
         {2, 0},
         {1, 4},
         {1, 3},
         {1, 0}});

    Architecture arc({{1, 0}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {4, 3}});
    RoutingConfig default_config(50, 0, 0, 0);
    QubitGraph q_graph =
        monomorph_interaction_graph(circ, arc.n_connections(), 5);
    std::vector<qubit_bimap_t> potential_maps =
        monomorphism_edge_break(arc, q_graph, 10000, 60000);

    qubit_mapping_t init_map = bimap_to_map(potential_maps[0].left);
    Placement pl(arc);
    pl.place_with_map(circ, init_map);
    Routing router(circ, arc);
    std::pair<Circuit, bool> out_circ = router.solve(default_config);
    REQUIRE(
        respects_connectivity_constraints(out_circ.first, arc, false) == true);
    REQUIRE(out_circ.second);
  }
  GIVEN(
      "The circuit left after clifford_simp, without clifford simp "
      "applied") {
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
    RoutingConfig default_config(50, 0, 0, 0);
    QubitGraph q_graph =
        monomorph_interaction_graph(circ, arc.n_connections(), 5);
    std::vector<qubit_bimap_t> potential_maps =
        monomorphism_edge_break(arc, q_graph, 10000, 60000);
    qubit_mapping_t init_map = bimap_to_map(potential_maps[0].left);
    Placement pl(arc);
    pl.place_with_map(circ, init_map);
    Routing router(circ, arc);
    std::pair<Circuit, bool> out_circ = router.solve(default_config);
    REQUIRE(
        respects_connectivity_constraints(out_circ.first, arc, false) == true);
    REQUIRE(out_circ.second);
  }
  GIVEN(
      "A smaller circuit that once had a segmentation fault when iterating "
      "through commands after clifford_simp() is applied and routing "
      "completed.") {
    Circuit circ(5);
    add_2qb_gates(
        circ, OpType::CX,
        {{1, 2}, {0, 3}, {1, 4}, {0, 1}, {2, 0}, {0, 1}, {1, 0}});
    Transform::clifford_simp().apply(circ);
    Architecture arc({{1, 0}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {4, 3}});
    RoutingConfig default_config(50, 0, 0, 0);
    QubitGraph q_graph =
        monomorph_interaction_graph(circ, arc.n_connections(), 5);
    std::vector<qubit_bimap_t> potential_maps =
        monomorphism_edge_break(arc, q_graph, 10000, 60000);

    qubit_mapping_t init_map = bimap_to_map(potential_maps[0].left);

    Placement pl(arc);
    pl.place_with_map(circ, init_map);
    Routing router(circ, arc);
    Circuit out_circ = router.solve(default_config).first;
    REQUIRE(respects_connectivity_constraints(out_circ, arc, false) == true);
  }
  GIVEN("The circuit that dies with clifford_simp") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}, {1, 0}, {2, 1}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.replace_SWAPs();
    Architecture arc({{1, 0}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {4, 3}});
    RoutingConfig default_config(50, 0, 0, 0);
    GraphPlacement pl(arc);
    qubit_mapping_t pl_map = pl.get_placement_map(circ);
    pl.place(circ);
    Routing router(circ, arc);
    Circuit out_circ = router.solve(default_config).first;
    qubit_mapping_t map = router.return_final_map();
    Vertex x = out_circ.add_op<Qubit>(OpType::X, {map.at(pl_map.at(Qubit(4)))});
    Vertex pred = out_circ.get_predecessors(x).front();
    REQUIRE(out_circ.get_OpType_from_Vertex(pred) == OpType::Z);
    REQUIRE(NoWireSwapsPredicate().verify(out_circ));
    REQUIRE(respects_connectivity_constraints(out_circ, arc, false) == true);
  }
}

SCENARIO(
    "Does the decompose_BRIDGE_gates function correctly decompose the "
    "BRIDGE Op, and pick the correct decomposition given the structure of "
    "surrounding CX gates?") {
  GIVEN("A single BRIDGE gate to be decomposed.") {
    Architecture test_arc({{0, 1}, {1, 2}});
    Circuit test_pc(3);
    test_pc.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
    Transform::decompose_BRIDGE_to_CX().apply(test_pc);
    auto it = test_pc.begin();
    unit_vector_t opt1 = {Qubit(0), Qubit(1)};
    unit_vector_t opt2 = {Qubit(1), Qubit(2)};
    CHECK(it->get_op_ptr()->get_type() == OpType::CX);
    CHECK(it->get_args() == opt2);
    ++it;
    CHECK(it->get_op_ptr()->get_type() == OpType::CX);
    CHECK(it->get_args() == opt1);
    ++it;
    CHECK(it->get_op_ptr()->get_type() == OpType::CX);
    CHECK(it->get_args() == opt2);
    ++it;
    CHECK(it->get_op_ptr()->get_type() == OpType::CX);
    CHECK(it->get_args() == opt1);
  }
  GIVEN("MultpleBRIDGE gate to be decomposed.") {
    Architecture test_arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    Circuit test_circuit(6);
    test_circuit.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
    test_circuit.add_op<unsigned>(OpType::BRIDGE, {1, 2, 3});
    test_circuit.add_op<unsigned>(OpType::BRIDGE, {2, 1, 0});
    test_circuit.add_op<unsigned>(OpType::BRIDGE, {2, 3, 4});
    test_circuit.add_op<unsigned>(OpType::BRIDGE, {3, 4, 5});
    Circuit test_pc(test_circuit);
    Transform::decompose_BRIDGE_to_CX().apply(test_pc);
    REQUIRE(test_pc.n_gates() == 20);
  }
}

SCENARIO("Does the rerouting of a solved circuit return 'false'?") {
  GIVEN("A simple circuit using default solve.") {
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
    Routing router(circ, arc);
    std::pair<Circuit, bool> out_circ = router.solve();
    REQUIRE(out_circ.second == true);
    Routing router2(out_circ.first, arc);
    std::pair<Circuit, bool> test_out2 = router2.solve();
    REQUIRE(test_out2.second == false);
    Routing router3(test_out2.first, arc);
    std::pair<Circuit, bool> test_out3 = router3.solve();
    REQUIRE(test_out3.second == false);
  }
  GIVEN("A simple circuit, but using a custom map for finding a solution.") {
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
    RoutingConfig default_config(50, 0, 0, 0);
    QubitGraph q_graph =
        monomorph_interaction_graph(circ, arc.n_connections(), 5);
    std::vector<qubit_bimap_t> potential_maps =
        monomorphism_edge_break(arc, q_graph, 10000, 60000);
    qubit_mapping_t init_map = bimap_to_map(potential_maps[0].left);
    Placement pl(arc);
    pl.place_with_map(circ, init_map);
    Routing router(circ, arc);
    std::pair<Circuit, bool> out_circ = router.solve(default_config);
    REQUIRE(
        respects_connectivity_constraints(out_circ.first, arc, false) == true);
    REQUIRE(out_circ.second);

    // Now try repeating it, making sure returned bool changes
    // make a LinePlacement plaer for the architecture
    LinePlacement lp_d(arc);

    Circuit c0 = out_circ.first;
    qubit_mapping_t m_0 = lp_d.get_placement_map(c0);
    lp_d.place_with_map(c0, m_0);
    Routing router2(c0, arc);
    std::pair<Circuit, bool> test_out2 = router2.solve();

    Circuit c1 = test_out2.first;
    REQUIRE(test_out2.second == true);
    Routing router3(c1, arc);
    qubit_vector_t pre_c1 = c1.all_qubits();
    std::pair<Circuit, bool> test_out3 = router3.solve();
    qubit_vector_t post_c1 = test_out3.first.all_qubits();
    REQUIRE(test_out3.second == false);
    Circuit c2 = test_out3.first;
    Routing router4(c2, arc);
    std::pair<Circuit, bool> test_out4 = router4.solve();
    REQUIRE(test_out4.second == false);
  }
}
SCENARIO("Routing on architecture with non-contiguous qubit labels") {
  GIVEN("A 2-qubit architecture with a gap") {
    Architecture arc(std::vector<std::pair<unsigned, unsigned>>{{0, 2}});
    PassPtr pass = gen_default_mapping_pass(arc);
    Circuit circ(2);
    CompilationUnit cu(circ);
    pass->apply(cu);
  }
  GIVEN("A 2-qubit architecture with a gap and some two-qubit gates") {
    Architecture arc(std::vector<std::pair<unsigned, unsigned>>{{0, 2}});
    PassPtr pass = gen_default_mapping_pass(arc);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    CompilationUnit cu(circ);
    pass->apply(cu);
  }
}

SCENARIO("Routing of aas example") {
  GIVEN("aas routing - simple example") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example II") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example III") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example IV") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example V") {
    Architecture arc(std::vector<Connection>{{Node(0), Node(1)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();
    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example VI") {
    Architecture arc(std::vector<Connection>{{Node(0), Node(2)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));

    const auto s = tket_sim::get_unitary(circ);
    const auto s1 = tket_sim::get_unitary(result);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        s, s1, tket_sim::MatrixEquivalence::EQUAL));
  }
  GIVEN("aas routing - simple example VII") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(2)}, {Node(2), Node(4)}, {Node(4), Node(6)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));

    const auto s = tket_sim::get_unitary(circ);
    const auto s1 = tket_sim::get_unitary(result);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        s, s1, tket_sim::MatrixEquivalence::EQUAL));
  }
  GIVEN("aas routing - simple example VIII") {
    Architecture arc(std::vector<Connection>{
        {Node(1000), Node(10)}, {Node(10), Node(100)}, {Node(100), Node(1)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {3});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - simple example IX, other gate set") {
    Architecture arc(std::vector<Connection>{
        {Node(1000), Node(10)}, {Node(10), Node(100)}, {Node(100), Node(1)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::X, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {3});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::X, {3});

    CompilationUnit cu(circ);

    REQUIRE(pass->apply(cu));

    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing with measure") {
    Architecture arc(std::vector<Connection>{{Node(0), Node(2)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(2, 2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    for (unsigned mes = 0; mes < 2; ++mes) {
      circ.add_measure(mes, mes);
    }

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
  }
  GIVEN("aas routing - circuit with fewer qubits then nodes in the arch") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("aas routing - circuit with fewer qubits then nodes in the arch II") {
    Architecture arc(std::vector<Connection>{
        {Node(0), Node(1)},
        {Node(1), Node(2)},
        {Node(2), Node(3)},
        {Node(3), Node(4)}});
    PassPtr pass = gen_full_mapping_pass_phase_poly(arc);
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});

    CompilationUnit cu(circ);
    REQUIRE(pass->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
}

SCENARIO("Routing preserves the number of qubits") {
  std::vector<std::pair<Node, Node>> cons;
  cons.push_back({Node("x", 1), Node("x", 0)});
  cons.push_back({Node("x", 2), Node("x", 1)});
  Architecture arc(
      std::vector<std::pair<Node, Node>>(cons.begin(), cons.end()));
  PassPtr pass = gen_default_mapping_pass(arc);
  Circuit c(3);
  c.add_op<unsigned>(OpType::CnX, {2, 1});
  CompilationUnit cu(c);
  bool applied = pass->apply(cu);
  const Circuit &c1 = cu.get_circ_ref();
  REQUIRE(c.n_qubits() == c1.n_qubits());
}

SCENARIO(
    "Methods related to correct routing and decomposition of circuits with "
    "classical wires.") {
  GIVEN(
      "A circuit with classical wires on CX gates. No Bridge gate "
      "allowed.") {
    Architecture test_arc({{0, 1}, {1, 2}});
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {2, 1}, {0, 1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 2);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {2, 1}, {1, 0}, 3);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 2}, {0, 1}, 0);
    Routing test_router(circ, test_arc);
    std::pair<Circuit, bool> output = test_router.solve({50, 0, 0, 0});
    Transform::decompose_SWAP_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(
        output.first, test_arc, false, false));
    Transform::decompose_BRIDGE_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(
        output.first, test_arc, false, false));
  }
  GIVEN(
      "A circuit that requires modification to satisfy architecture "
      "constraints.") {
    Architecture sg({{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    Circuit circ(5, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {1, 3}, {1, 4}, {0, 1}});
    Routing test_router(circ, sg);
    std::pair<Circuit, bool> output = test_router.solve({50, 0, 0, 0});
    Transform::decompose_SWAP_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, sg, false, false));
    Transform::decompose_BRIDGE_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, sg, false, false));
    Command classical_com = output.first.get_commands()[0];
    REQUIRE(classical_com.get_args()[0] == output.first.all_bits()[0]);
  }
  GIVEN("A single Bridge gate with multiple classical wires, decomposed.") {
    Architecture arc({{0, 1}, {1, 2}});
    Circuit circ(3, 3);
    circ.add_conditional_gate<unsigned>(
        OpType::BRIDGE, {}, {0, 1, 2}, {0, 1, 2}, 1);
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
    Transform::decompose_BRIDGE_to_CX().apply(circ);
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
    Transform::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true, false));
    std::vector<Command> all_coms = circ.get_commands();
    REQUIRE(all_coms[0].get_args()[0] == circ.all_bits()[1]);
    REQUIRE(all_coms[0].get_args()[1] == circ.all_bits()[0]);
    REQUIRE(all_coms[1].get_args()[0] == circ.all_bits()[0]);
    REQUIRE(all_coms[1].get_args()[1] == circ.all_bits()[1]);
  }
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
    Routing router(circ, arc);
    std::pair<Circuit, bool> output = router.solve();
    Transform::decompose_SWAP_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, arc, false, true));
    Transform::decompose_BRIDGE_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, arc, false, true));
  }
  GIVEN(
      "A large circuit, with a mixture of conditional CX and CX gates with "
      "multiple classical wires, non conditional CX and, single qubit "
      "gates, and a directed architecture.") {
    SquareGrid arc(10, 4, 2);
    Circuit circ(60, 10);
    for (unsigned i = 0; i < 58; i++) {
      circ.add_op<unsigned>(OpType::CX, {i, i + 1});
      circ.add_conditional_gate<unsigned>(
          OpType::CX, {}, {i + 2, i}, {0, 2, 3, 5}, 1);
      circ.add_conditional_gate<unsigned>(OpType::H, {}, {i}, {0, 7}, 1);
      circ.add_conditional_gate<unsigned>(
          OpType::CX, {}, {i + 2, i + 1}, {1, 2, 3, 5, 9}, 0);
      circ.add_conditional_gate<unsigned>(OpType::S, {}, {i + 1}, {1, 2, 7}, 1);
      circ.add_conditional_gate<unsigned>(
          OpType::CX, {}, {i, i + 1}, {4, 6, 8, 7, 9}, 0);
      circ.add_conditional_gate<unsigned>(OpType::X, {}, {i + 2}, {0, 3}, 0);
    }
    Routing router(circ, arc);
    std::pair<Circuit, bool> output = router.solve();
    Transform::decompose_SWAP_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, arc, false, true));
    Transform::decompose_BRIDGE_to_CX().apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, arc, false, true));
    Transform::decompose_CX_directed(arc).apply(output.first);
    REQUIRE(respects_connectivity_constraints(output.first, arc, true, true));
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
    Routing router(circ, arc);
    Circuit c = router.solve().first;
    Transform T_1 = Transform::decompose_SWAP_to_CX();
    T_1.apply(c);
    REQUIRE(c.count_gates(OpType::SWAP) == 0);
  }
}

SCENARIO("Does add_distributed_cx account for incorrect BRIDGE nodes?") {
  GIVEN("An incorrect and a correct BRIDGE orientation.") {
    Architecture a({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    Circuit c(6);
    c.add_op<unsigned>(OpType::CX, {3, 5});
    c.add_op<unsigned>(OpType::CX, {2, 0});

    Placement placer(a);
    qubit_vector_t c_qubits = c.all_qubits();
    node_vector_t a_nodes = a.get_all_uids_vec();

    qubit_mapping_t initial_map = {
        {c_qubits[0], a_nodes[0]}, {c_qubits[1], a_nodes[1]},
        {c_qubits[2], a_nodes[2]}, {c_qubits[3], a_nodes[3]},
        {c_qubits[4], a_nodes[4]}, {c_qubits[5], a_nodes[5]}};

    placer.place_with_map(c, initial_map);

    Routing r(c, a);
    RoutingTester rt(&r);

    rt.initialise_slicefrontier();
    qubit_bimap_t qbm;
    for (unsigned nn = 0; nn <= 5; ++nn) {
      qbm.insert({a_nodes[nn], Node(nn)});
    }

    rt.set_qmap(qbm);

    rt.add_distributed_cx(Node(5), Node(3), Node(4));
    rt.add_distributed_cx(Node(2), Node(0), Node(1));

    std::vector<Command> bridge_commands = rt.get_circ()->get_commands();
    qubit_vector_t com_0_qubits = {a_nodes[2], a_nodes[1], a_nodes[0]};
    qubit_vector_t com_1_qubits = {a_nodes[3], a_nodes[4], a_nodes[5]};
    REQUIRE(bridge_commands[0].get_qubits() == com_0_qubits);
    REQUIRE(bridge_commands[1].get_qubits() == com_1_qubits);
  }
  GIVEN("An invalid BRIDGE.") {
    Architecture a({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    Circuit c(6);
    c.add_op<unsigned>(OpType::CX, {2, 5});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    Placement placer(a);
    qubit_vector_t c_qubits = c.all_qubits();
    node_vector_t a_nodes = a.get_all_uids_vec();

    qubit_mapping_t initial_map = {
        {c_qubits[0], a_nodes[0]}, {c_qubits[1], a_nodes[1]},
        {c_qubits[2], a_nodes[2]}, {c_qubits[3], a_nodes[3]},
        {c_qubits[4], a_nodes[4]}, {c_qubits[5], a_nodes[5]}};

    placer.place_with_map(c, initial_map);

    Routing r(c, a);
    RoutingTester rt(&r);

    rt.initialise_slicefrontier();
    qubit_bimap_t qbm;
    for (unsigned nn = 0; nn <= 5; ++nn) {
      qbm.insert({a_nodes[nn], Node(nn)});
    }

    rt.set_qmap(qbm);

    REQUIRE_THROWS_AS(
        rt.add_distributed_cx(Node(2), Node(4), Node(5)), BridgeInvalid);
    REQUIRE_THROWS_AS(
        rt.add_distributed_cx(Node(0), Node(1), Node(3)), BridgeInvalid);
    REQUIRE_THROWS_AS(
        rt.add_distributed_cx(Node(0), Node(1), Node(2)), BridgeInvalid);
  }
}

}  // namespace test_Routing
}  // namespace tket
