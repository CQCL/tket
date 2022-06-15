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
#include <random>

#include "../testutil.hpp"
#include "Placement/Placement.hpp"

namespace tket {
namespace test_Placement {

using Connection = Architecture::Connection;

// confirm all weight i+1 edges are broken before any weight i edge
bool check_edge_break_order(
    const Architecture& arc, const QubitGraph& qg, const qubit_bimap_t& map) {
  std::vector<unsigned> mapped_edge_weights, unmapped_edge_weights;
  for (auto [q1, q2] : qg.get_all_edges_vec()) {
    const unsigned weight = qg.get_connection_weight(q1, q2);
    if (weight >= unmapped_edge_weights.size()) {
      unmapped_edge_weights.resize(weight + 1, 0);
    }
    if (map.left.find(q1) == map.left.end() ||
        map.left.find(q2) == map.left.end()) {
      unmapped_edge_weights[weight]++;
    } else {
      const Node n1 = map.left.at(q1), n2 = map.left.at(q2);
      if (!(arc.edge_exists(n1, n2) || arc.edge_exists(n2, n1))) {
        unmapped_edge_weights[weight]++;
      }
    }
  }

  bool zeros_finished = false;
  for (unsigned x : unmapped_edge_weights) {
    if (zeros_finished) {
      if (x == 0) return false;
    } else {
      if (x > 0) zeros_finished = true;
    }
  }
  return true;
}
// monomorphism_edge_break
SCENARIO("Small monomorphisms. place_qubit") {
  GIVEN("2 qubit graphs") {
    std::vector<std::pair<unsigned, unsigned>> edges{{0, 1}};
    Architecture arc(edges);
    Circuit circ(2);
    qubit_vector_t qbs = circ.all_qubits();
    std::map<Node, Qubit> node_to_qubit;
    for (unsigned i = 0; i < 2; ++i) {
      node_to_qubit.insert({Node(i), qbs[i]});
    }
    QubitGraph qg(qbs);
    qg.add_connection(qbs[0], qbs[1], 1);
    std::vector<qubit_bimap_t> result =
        monomorphism_edge_break(arc, qg, 10, 60000);
    REQUIRE(result[0].size() == 2);
    // REQUIRE(result.first.size()==boost::num_vertices(inter));
    for (const auto& [qb, n] : result[0].left) {
      REQUIRE(qb == node_to_qubit[n]);
    }
    for (const qubit_bimap_t& map : result) {
      REQUIRE(check_edge_break_order(arc, qg, map));
    }

    WHEN("The architecture has no edges, but two nodes") {
      arc.remove_connection({Node(0), Node(1)});
      std::vector<qubit_bimap_t> result2 =
          monomorphism_edge_break(arc, qg, 10, 60000);
      THEN(
          "Interaction edge broken, both nodes removed, and an empty "
          "map returned.") {
        for (const qubit_bimap_t& map : result2) {
          REQUIRE(map.empty());
        }
      }
    }
  }
  GIVEN("4 qubit graphs") {
    Architecture arc({{0, 1}, {1, 2}, {2, 0}, {1, 3}});
    Circuit circ(4);
    qubit_vector_t qbs = circ.all_qubits();
    std::map<Node, Qubit> node_to_qubit;
    for (unsigned i = 0; i < 4; ++i) {
      node_to_qubit.insert({Node(i), qbs[i]});
    }
    QubitGraph qg(qbs);
    qg.add_connection(qbs[0], qbs[1], 1);
    qg.add_connection(qbs[1], qbs[2], 1);
    qg.add_connection(qbs[2], qbs[0], 1);
    qg.add_connection(qbs[1], qbs[3], 1);

    std::vector<qubit_bimap_t> result =
        monomorphism_edge_break(arc, qg, 10, 60000);
    // REQUIRE(result.second==0);
    REQUIRE(result[0].size() == qg.n_nodes());
    for (const auto& [qb, n] : result[0].left) {
      REQUIRE(qb == node_to_qubit[n]);
    }
    for (const qubit_bimap_t& map : result) {
      REQUIRE(check_edge_break_order(arc, qg, map));
    }

    WHEN("Remove an edge") {
      arc.remove_connection({Node(1), Node(2)});
      std::vector<qubit_bimap_t> result2 =
          monomorphism_edge_break(arc, qg, 10, 60000);
      for (const qubit_bimap_t& map : result2) {
        REQUIRE(check_edge_break_order(arc, qg, map));
      }

      THEN(
          "can still find mapping, but requires edge removal from "
          "interaction graph") {
        REQUIRE(
            result2[0].size() == 4);  // can still find mapping, but requires
                                      // edge removal from interaction graph;
      }
    }
    WHEN("Remove a different edge") {
      qg.remove_connection({qbs[2], qbs[0]});
      std::vector<qubit_bimap_t> result2 =
          monomorphism_edge_break(arc, qg, 10, 60000);
      // REQUIRE(result2.second==0);
      REQUIRE(result2[0].size() == qg.n_nodes());
      for (const qubit_bimap_t& map : result2) {
        REQUIRE(check_edge_break_order(arc, qg, map));
      }
      for (const auto& [qb, n] : result2[0].left) {
        REQUIRE(qb == node_to_qubit[n]);
      }
    }
  }
  GIVEN("Interaction graphs that don't fit on architecture") {
    SquareGrid arc(3, 4);
    Circuit circ(10);
    qubit_vector_t qbs = circ.all_qubits();
    QubitGraph qg(qbs);

    for (unsigned slice = 1; slice <= 4; slice++) {
      std::shuffle(qbs.begin(), qbs.end(), std::default_random_engine(slice));
      for (unsigned i = 0; i < 6; i++) {
        qg.add_connection(qbs[i], qbs[i + 1], slice);
      }
    }

    std::vector<qubit_bimap_t> result =
        monomorphism_edge_break(arc, qg, 10, 60000);
    for (const qubit_bimap_t& map : result) {
      REQUIRE(check_edge_break_order(arc, qg, map));
    }
  }
}

// Monomorpher
SCENARIO("Check Monomorpher satisfies correct placement conditions") {
  GIVEN("A simple architecture.") {
    Architecture arc({{0, 1}, {1, 2}});
    WHEN("A depth 1 circuit which fits is placed") {
      Circuit test_circ(3);
      test_circ.add_op<unsigned>(OpType::T, {1});
      add_2qb_gates(test_circ, OpType::CX, {{2, 0}, {0, 1}});

      Monomorpher morph(test_circ, arc, {}, {3, arc.n_connections()});
      qubit_mapping_t map = morph.place(1)[0].map;
      // qubit_mapping_t expected_map = {{2, 0}, {0, 1}, {1, 2}};
      // print_map(map);
      THEN(
          "All qubits are placed, with the connecting qubit on "
          "connecting node") {
        REQUIRE(map[Qubit(q_default_reg(), 0)] == Node(1));
        REQUIRE(map.size() == 3);
      }
    }
    WHEN("A depth 2 circuit which fits is placed") {
      Circuit test_circ(3);
      test_circ.add_op<unsigned>(OpType::T, {1});
      test_circ.add_op<unsigned>(OpType::CX, {2, 0});
      test_circ.add_op<unsigned>(OpType::CX, {0, 2});
      test_circ.add_op<unsigned>(OpType::S, {0});
      test_circ.add_op<unsigned>(OpType::CX, {2, 1});
      // test_circ.add_op<unsigned>(OpType::CX, {0, 1});

      Monomorpher morph(test_circ, arc, {}, {3, arc.n_connections()});
      qubit_mapping_t map = morph.place(1)[0].map;
      // qubit_mapping_t expected_map = {{2, 0}, {0, 1}, {1, 2}};
      // print_map(map);
      THEN(
          "All qubits are placed, with the connecting qubit on "
          "connecting node") {
        REQUIRE(map[Qubit(q_default_reg(), 2)] == Node(1));
        REQUIRE(map.size() == 3);
      }
    }
  }
  GIVEN("A linear architecture") {
    std::vector<Connection> edges = {
        {Node(1), Node(2)}, {Node(0), Node(1)}, {Node(2), Node(3)}};
    Architecture arc(edges);
    WHEN("A node needs to be removed for placement.") {
      Circuit test_circ(4);
      test_circ.add_op<unsigned>(OpType::T, {1});
      test_circ.add_op<unsigned>(OpType::CX, {0, 1});
      test_circ.add_op<unsigned>(OpType::CX, {1, 3});
      test_circ.add_op<unsigned>(OpType::S, {0});
      test_circ.add_op<unsigned>(OpType::CX, {3, 0});
      test_circ.add_op<unsigned>(OpType::CX, {2, 1});
      // test_circ.add_op<unsigned>(OpType::CX, {0, 1});

      Monomorpher morph(test_circ, arc, {}, {4, 5});
      qubit_mapping_t map = morph.place(1)[0].map;
      // qubit_mapping_t expected_map = {{2, 0}, {0, 1}, {1, 2}};
      THEN("Only 3 qubits are placed, and the correct one is removed") {
        REQUIRE(map.size() == 3);
        REQUIRE(map.find(Qubit(q_default_reg(), 2)) == map.end());
      }
    }
    WHEN("Directness is specified via edge error rate.") {
      Circuit test_circ(4);
      add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {1, 3}, {3, 0}});
      avg_link_errors_t link_errors;
      gate_error_t CX_error_good(0.1);
      gate_error_t CX_error_bad(1.0 - (0.9 * 0.99 * 0.99 * 0.99 * 0.99));
      for (unsigned i = 0; i < edges.size(); i++) {
        link_errors.insert({{edges[i].first, edges[i].second}, CX_error_good});
        link_errors.insert({{edges[i].second, edges[i].first}, CX_error_bad});
      }

      Monomorpher morph(test_circ, arc, {{}, link_errors}, {4, 5});
      qubit_mapping_t map = morph.place(1)[0].map;
      REQUIRE(map.size() == 3);

      THEN("The chosen map satisfies directionality") {
        qubit_vector_t qbs = test_circ.all_qubits();
        REQUIRE(arc.edge_exists(map[qbs[0]], map[qbs[1]]));
        REQUIRE(!arc.edge_exists(map[qbs[1]], map[qbs[0]]));
        REQUIRE(arc.edge_exists(map[qbs[1]], map[qbs[3]]));
        REQUIRE(!arc.edge_exists(map[qbs[3]], map[qbs[1]]));
      }
    }

    WHEN(
        "The circuit is two qubits and there is a preferred edge "
        "fidelity.") {
      Circuit test_circ(2);
      test_circ.add_op<unsigned>(OpType::CX, {0, 1});
      avg_link_errors_t link_errors;
      for (unsigned i = 0; i < edges.size(); i++) {
        gate_error_t CX_error_good(0.1 - i * 0.01);
        gate_error_t CX_error_bad(0.1 + i * 0.01);
        link_errors.insert({{edges[i].first, edges[i].second}, CX_error_good});
        link_errors.insert({{edges[i].second, edges[i].first}, CX_error_bad});
      }

      Monomorpher morph(test_circ, arc, {{}, link_errors}, {4, 5});
      qubit_mapping_t map = morph.place(1)[0].map;
      THEN("The circuit is placed on the best edge.") {
        qubit_vector_t qbs = test_circ.all_qubits();
        REQUIRE(map.size() == 2);
        REQUIRE(map[qbs[0]] == edges[edges.size() - 1].first);
        REQUIRE(map[qbs[1]] == edges[edges.size() - 1].second);
      }
    }
    WHEN(
        "The circuit is two qubits and there is a preferred edge fidelity "
        "direction.") {
      Circuit test_circ(2);
      test_circ.add_op<unsigned>(OpType::CX, {0, 1});
      avg_link_errors_t link_errors;
      for (unsigned i = 0; i < edges.size() - 1; i++) {
        gate_error_t CX_error_good(0.1 - 1 * 0.01);
        gate_error_t CX_error_bad(0.1 + i * 0.01);
        link_errors.insert({{edges[i].first, edges[i].second}, CX_error_good});
        link_errors.insert({{edges[i].second, edges[i].first}, CX_error_bad});
      }

      Monomorpher morph(test_circ, arc, {{}, link_errors}, {4, 5});
      qubit_mapping_t map = morph.place(1)[0].map;
      THEN("The circuit is placed on the best edge.") {
        REQUIRE(map.size() == 2);
        qubit_vector_t qbs = test_circ.all_qubits();
        const Architecture::Connection mapped_edge = {map[qbs[0]], map[qbs[1]]};
        const Architecture::Connection reversed_edge = {
            mapped_edge.second, mapped_edge.first};
        const Architecture::Connection preferred_edge = {
            edges[edges.size() - 1].first, edges[edges.size() - 1].second};
        bool edge_equal = mapped_edge == preferred_edge;
        edge_equal |= (reversed_edge == preferred_edge);
        REQUIRE(edge_equal);
      }
    }

    WHEN(
        "The circuit is two qubits and there is a preferred edge by "
        "single qubit error.") {
      Circuit test_circ(2);
      test_circ.add_op<unsigned>(OpType::CX, {0, 1});
      avg_link_errors_t link_errors;
      gate_error_t CX_error_good(0.05);
      gate_error_t CX_error_bad(0.09);

      for (unsigned i = 0; i < edges.size(); i++) {
        link_errors.insert({{edges[i].first, edges[i].second}, CX_error_good});
        link_errors.insert({{edges[i].second, edges[i].first}, CX_error_bad});
      }

      gate_error_t single_error(0.01);
      gate_error_t gd_single_error(0.001);

      avg_node_errors_t node_errors;

      node_errors.insert({Node(0), single_error});
      node_errors.insert({Node(2), gd_single_error});
      node_errors.insert({Node(3), single_error});
      node_errors.insert({Node(1), gd_single_error});

      Monomorpher morph(test_circ, arc, {node_errors, link_errors}, {4, 5});
      qubit_mapping_t map = morph.place(1)[0].map;
      THEN("The circuit is placed on the best edge.") {
        qubit_vector_t qbs = test_circ.all_qubits();
        REQUIRE(map.size() == 2);
        REQUIRE(map[qbs[0]] == Node(1));
        REQUIRE(map[qbs[1]] == Node(2));
      }
    }

    WHEN(
        "The circuit is two qubits, shallow and there is a preferred edge "
        "by readout error.") {
      Circuit test_circ(3);
      const std::vector<std::pair<unsigned, unsigned>> inters{{0, 1}, {2, 0}};
      add_2qb_gates(test_circ, OpType::CX, inters);
      avg_link_errors_t link_errors;
      std::vector<double> cx_errs = {0.013, 0.008, 0.01};
      for (unsigned i = 0; i < edges.size(); i++) {
        gate_error_t error_cont = cx_errs[i];
        link_errors.insert({{edges[i].first, edges[i].second}, error_cont});
        link_errors.insert({{edges[i].second, edges[i].first}, error_cont});
      }

      gate_error_t bd_single_gate_error = 0.12;
      gate_error_t gd_single_gate_error = 0.013;
      avg_readout_errors_t readout_errors;

      readout_errors.insert({Node(0), bd_single_gate_error});
      readout_errors.insert({Node(1), bd_single_gate_error});
      readout_errors.insert({Node(2), gd_single_gate_error});
      readout_errors.insert({Node(3), gd_single_gate_error});

      THEN("The circuit is placed on the best edge.") {
        Monomorpher morph(
            test_circ, arc, {{}, link_errors, readout_errors}, {4, 5});
        qubit_mapping_t map = morph.place(1)[0].map;
        qubit_vector_t qbs = test_circ.all_qubits();

        REQUIRE(map.size() == 3);
        CHECK(map[qbs[0]] == Node(2));
        CHECK(map[qbs[1]] == Node(1));
        CHECK(map[qbs[2]] == Node(3));
      }

      AND_WHEN("The circuit is made deeper") {
        for (unsigned i = 0; i < 20; i++) {
          add_2qb_gates(test_circ, OpType::CX, inters);
        }

        THEN("CX errors are preferentially considered.") {
          Monomorpher morph(
              test_circ, arc, {{}, link_errors, readout_errors}, {4, 5});
          qubit_mapping_t map = morph.place(1)[0].map;
          qubit_vector_t qbs = test_circ.all_qubits();

          REQUIRE(map.size() == 3);
          CHECK(map[qbs[0]] == Node(1));
          CHECK(map[qbs[1]] == Node(2));
          CHECK(map[qbs[2]] == Node(0));
        }
      }

      AND_WHEN("Readout error differences are small") {
        gate_error_t bd_single_gate_error2 = 0.014;
        avg_readout_errors_t readout_errors2;

        readout_errors2.insert({Node(0), bd_single_gate_error2});
        readout_errors2.insert({Node(1), bd_single_gate_error2});
        readout_errors2.insert({Node(2), gd_single_gate_error});
        readout_errors2.insert({Node(3), gd_single_gate_error});

        THEN("CX errors are preferentially considered.") {
          Monomorpher morph(
              test_circ, arc, {{}, link_errors, readout_errors2}, {4, 5});
          qubit_mapping_t map = morph.place(1)[0].map;
          qubit_vector_t qbs = test_circ.all_qubits();

          REQUIRE(map.size() == 3);
          CHECK(map[qbs[0]] == Node(1));
          CHECK(map[qbs[1]] == Node(2));
          CHECK(map[qbs[2]] == Node(0));
        }
      }
    }

    WHEN(
        "A deep circuit is placed on an architecture with a highly "
        "connected region.") {
      arc = SquareGrid(4, 4);
      Circuit test_circ(4);
      std::vector<unsigned> qblist;
      qblist.resize(4);
      std::iota(qblist.begin(), qblist.end(), 0);
      // make circuit big enouogh to cause architecture constriction
      for (unsigned slice = 1; slice <= 21; slice++) {
        // WARNING: the detailed shuffle can differ across
        // different platforms, compilers, AND compiler versions!
        std::shuffle(
            qblist.begin(), qblist.end(), std::default_random_engine(slice));
        for (unsigned i = 0; i < 2; i++) {
          test_circ.add_op<unsigned>(OpType::CX, {qblist[i], qblist[i + 1]});
        }
      }

      Monomorpher morph(test_circ, arc, {}, {10, arc.n_connections()});
    }
  }
}

SCENARIO(
    "Does 'noise aware placement' deal with an undirected architecture "
    "(i.e. a coupling list with {0,1} and {1,0})?") {
  GIVEN(
      "A small undirected architecture, the graph placement method, a "
      "basic CX circuit.") {
    Circuit test_circ(2);
    add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {1, 0}});

    Architecture test_arc({{0, 1}, {1, 0}});
    QubitGraph q_graph =
        monomorph_interaction_graph(test_circ, test_arc.n_connections(), 5);
    std::vector<qubit_bimap_t> potential_maps =
        monomorphism_edge_break(test_arc, q_graph, 10000, 60000);
    REQUIRE(potential_maps.size() > 0);
  }
  GIVEN("A much larger example.") {
    Circuit test_circ(10);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1},
         {1, 0},
         {2, 0},
         {3, 5},
         {7, 8},
         {7, 9},
         {6, 1},
         {4, 6},
         {0, 7},
         {1, 5},
         {2, 4},
         {9, 5},
         {7, 6},
         {1, 9},
         {0, 4},
         {3, 4}});

    Architecture test_arc({{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2},
                           {3, 4}, {4, 3}, {4, 5}, {5, 4}, {5, 6}, {6, 5},
                           {6, 7}, {7, 6}, {7, 8}, {8, 7}, {8, 9}, {9, 8},
                           {2, 4}, {4, 2}, {2, 6}, {6, 2}, {7, 1}, {1, 7},
                           {9, 2}, {2, 9}, {7, 9}, {9, 7}});
    QubitGraph q_graph =
        monomorph_interaction_graph(test_circ, test_arc.n_connections(), 5);
    std::vector<qubit_bimap_t> potential_maps =
        monomorphism_edge_break(test_arc, q_graph, 10000, 60000);
    REQUIRE(potential_maps.size() > 0);
  }
}
SCENARIO("Test NaivePlacement class") {
  Architecture test_arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}});
  GIVEN(
      "No Qubits placed in Circuit, same number of qubits and architecture "
      "nodes.") {
    Circuit test_circ(7);
    NaivePlacement np(test_arc);
    qubit_mapping_t p = np.get_placement_map(test_circ);
    REQUIRE(p[Qubit(0)] == Node(0));
    REQUIRE(p[Qubit(1)] == Node(1));
    REQUIRE(p[Qubit(2)] == Node(2));
    REQUIRE(p[Qubit(3)] == Node(3));
    REQUIRE(p[Qubit(4)] == Node(4));
    REQUIRE(p[Qubit(5)] == Node(5));
    REQUIRE(p[Qubit(6)] == Node(6));
  }
  GIVEN("No Qubits placed in Circuit, less qubits than architecture nodes.") {
    Circuit test_circ(6);
    NaivePlacement np(test_arc);
    qubit_mapping_t p = np.get_placement_map(test_circ);
    REQUIRE(p[Qubit(0)] == Node(0));
    REQUIRE(p[Qubit(1)] == Node(1));
    REQUIRE(p[Qubit(2)] == Node(2));
    REQUIRE(p[Qubit(3)] == Node(3));
    REQUIRE(p[Qubit(4)] == Node(4));
    REQUIRE(p[Qubit(5)] == Node(5));
  }
  GIVEN(
      "Some Qubits placed in Circuit, same number of qubits and architecture "
      "nodes.") {
    Circuit test_circ(4);
    test_circ.add_qubit(Node(0));
    test_circ.add_qubit(Node(1));
    test_circ.add_qubit(Node(2));
    NaivePlacement np(test_arc);
    qubit_mapping_t p = np.get_placement_map(test_circ);

    REQUIRE(p[Qubit(0)] == Node(3));
    REQUIRE(p[Qubit(1)] == Node(4));
    REQUIRE(p[Qubit(2)] == Node(5));
    REQUIRE(p[Qubit(3)] == Node(6));
    REQUIRE(p[Node(0)] == Node(0));
    REQUIRE(p[Node(1)] == Node(1));
    REQUIRE(p[Node(2)] == Node(2));
  }
  GIVEN("Some Qubits placed in Circuit, less qubits than architecture nodes.") {
    Circuit test_circ(2);
    test_circ.add_qubit(Node(0));
    test_circ.add_qubit(Node(1));
    test_circ.add_qubit(Node(2));
    NaivePlacement np(test_arc);
    qubit_mapping_t p = np.get_placement_map(test_circ);

    REQUIRE(p[Qubit(0)] == Node(3));
    REQUIRE(p[Qubit(1)] == Node(4));
    REQUIRE(p[Node(0)] == Node(0));
    REQUIRE(p[Node(1)] == Node(1));
    REQUIRE(p[Node(2)] == Node(2));
  }
}

// Tests for new placement method wrappers
SCENARIO(
    "Does the base Placement class correctly modify Circuits and return "
    "maps?") {
  Architecture test_arc({{0, 1}, {1, 2}, {2, 3}});
  Placement test_p(test_arc);
  Qubit uid0 = Qubit("unplaced", 0);
  Qubit uid1 = Qubit("unplaced", 1);
  Qubit uid2 = Qubit("unplaced", 2);
  Qubit uid3 = Qubit("unplaced", 3);
  const std::array<Qubit, 4> expected_qubits{uid0, uid1, uid2, uid3};

  GIVEN("A basic circuit and architecture. place method.") {
    Circuit test_circ(4);
    add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {2, 1}, {3, 1}});

    test_p.place(test_circ);
    const qubit_vector_t all_qs = test_circ.all_qubits();
    for (unsigned nn = 0; nn < expected_qubits.size(); ++nn) {
      REQUIRE(all_qs[nn] == expected_qubits[nn]);
    }
  }
  GIVEN("A basic circuit and architecture. get_placement_map method.") {
    Circuit test_circ(4);
    add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {2, 1}, {3, 1}});

    const qubit_mapping_t test_m = test_p.get_placement_map(test_circ);
    const qubit_vector_t all_qs = test_circ.all_qubits();
    for (unsigned nn = 0; nn < expected_qubits.size(); ++nn) {
      REQUIRE(test_m.at(all_qs[nn]) == expected_qubits[nn]);
    }
  }
}

SCENARIO(
    "Does the LinePlacement class correctly modify Circuits and return "
    "maps?") {
  Architecture test_arc({{0, 1}, {1, 2}, {2, 3}});
  LinePlacement test_p(test_arc);
  Qubit uid0 = Qubit("unplaced", 0);
  Qubit uid1 = Qubit("unplaced", 1);
  Qubit uid2 = Qubit("unplaced", 2);
  Qubit uid3 = Qubit("unplaced", 3);
  GIVEN("A basic circuit and architecture. place method.") {
    Circuit test_circ(4);
    add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {2, 1}, {3, 1}});

    test_p.place(test_circ);

    const qubit_vector_t all_qs = test_circ.all_qubits();
    for (unsigned nn = 0; nn <= 2; ++nn) {
      REQUIRE(all_qs[nn] == Node(nn + 1));
    }
    REQUIRE(all_qs[3] == uid0);
  }
  GIVEN("A basic circuit and architecture. get_placement_map method.") {
    Circuit test_circ(4);
    add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {2, 1}, {3, 1}});

    const qubit_mapping_t test_m = test_p.get_placement_map(test_circ);
    const qubit_vector_t all_qs = test_circ.all_qubits();
    for (unsigned nn = 0; nn <= 2; ++nn) {
      REQUIRE(test_m.at(all_qs[nn]) == Node(nn + 1));
    }
    REQUIRE(test_m.at(all_qs[3]) == uid0);
  }
}

SCENARIO(
    "Does the GraphPlacement class correctly modify Circuits and return "
    "maps?") {
  Architecture test_arc({{0, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 5}});
  GraphPlacement test_p(test_arc);
  Qubit uid0 = Qubit("unplaced", 0);
  Qubit uid1 = Qubit("unplaced", 1);
  Qubit uid2 = Qubit("unplaced", 2);
  Qubit uid3 = Qubit("unplaced", 3);
  Qubit uid4 = Qubit("unplaced", 4);
  Qubit uid5 = Qubit("unplaced", 5);

  GIVEN("A basic circuit and architecture. place method.") {
    Circuit test_circ(6);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1}, {2, 1}, {3, 1}, {2, 5}, {3, 4}, {0, 5}});

    test_p.place(test_circ);

    const qubit_vector_t all_qs = test_circ.all_qubits();

    const std::array<unsigned, 5> indices{0, 1, 2, 3, 5};
    for (unsigned nn = 0; nn < indices.size(); ++nn) {
      REQUIRE(all_qs[nn] == Node(indices[nn]));
    }
    REQUIRE(all_qs[5] == uid0);
  }
  GIVEN("A basic circuit and architecture. get_placement_map method.") {
    Circuit test_circ(6);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1}, {2, 1}, {3, 1}, {2, 5}, {3, 4}, {0, 5}});

    const qubit_mapping_t test_m = test_p.get_placement_map(test_circ);
    const qubit_vector_t all_qs = test_circ.all_qubits();
    REQUIRE(test_m.at(all_qs[4]) == uid0);
    const std::vector<std::pair<unsigned, unsigned>> mapping{
        {0, 0}, {1, 1}, {2, 2}, {3, 3}, {5, 5}};
    for (const auto& pair : mapping) {
      REQUIRE(test_m.at(all_qs[pair.first]) == Node(pair.second));
    }
  }
}

SCENARIO(
    "Does the timeout config option work as expected with monomorpher "
    "place?",
    "[.long]") {
  GIVEN("A large architecture, qubit graph and small timeout") {
    const SquareGrid arc(10, 10, 5);
    Circuit circ(40);
    for (unsigned i = 1; i < 39; i++) {
      add_2qb_gates(circ, OpType::CX, {{i, i + 1}, {i, i - 1}});
    }
    for (unsigned i = 3; i < 35; i++) {
      add_2qb_gates(circ, OpType::CX, {{i - 1, i + 2}, {i, i + 2}});
    }
    PlacementConfig pc;
    pc.depth_limit = 5;
    pc.max_interaction_edges = arc.n_connections();
    pc.vf2_max_matches = 10000000;
    pc.arc_contraction_ratio = 10;
    pc.timeout = 1000;

    GraphPlacement placer(arc);
    placer.set_config(pc);
    std::vector<qubit_mapping_t> all_maps = placer.get_all_placement_maps(circ);
    REQUIRE(all_maps.size() < pc.vf2_max_matches);
  }
}

SCENARIO(
    "Does the NoiseAwarePlacement class correctly modify Circuits and "
    "return "
    "maps?") {
  Architecture test_arc({{0, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 5}});
  NoiseAwarePlacement test_p(test_arc);
  Qubit uid0 = Qubit("unplaced", 0);
  Qubit uid1 = Qubit("unplaced", 1);
  Qubit uid2 = Qubit("unplaced", 2);
  Qubit uid3 = Qubit("unplaced", 3);
  Qubit uid4 = Qubit("unplaced", 4);
  Qubit uid5 = Qubit("unplaced", 5);

  GIVEN("A basic circuit and architecture. place method.") {
    Circuit test_circ(6);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1}, {2, 1}, {3, 1}, {2, 5}, {3, 4}, {0, 5}});

    qubit_vector_t pre_place = test_circ.all_qubits();
    test_p.place(test_circ);
    qubit_vector_t all_qs = test_circ.all_qubits();
    REQUIRE(pre_place != all_qs);
  }
  GIVEN("A basic circuit and architecture. get_placement_map method.") {
    Circuit test_circ(6);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1}, {2, 1}, {3, 1}, {2, 5}, {3, 4}, {0, 5}});

    qubit_vector_t pre_place = test_circ.all_qubits();
    qubit_mapping_t test_m = test_p.get_placement_map(test_circ);
    qubit_vector_t all_qs = test_circ.all_qubits();
    REQUIRE(pre_place == all_qs);
  }
}

}  // namespace test_Placement
}  // namespace tket
