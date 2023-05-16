// Copyright 2019-2023 Cambridge Quantum Computing
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
#include "tket/Placement/Placement.hpp"

namespace tket {

SCENARIO("Base GraphPlacement class") {
  GIVEN("Empty Architecture, GraphPlacement::GraphPlacement.") {
    Architecture architecture;
    REQUIRE_THROWS_AS(GraphPlacement(architecture), std::logic_error);
  }
  GIVEN("Empty circuit, two qubit Architecture, GraphPlacement::Place.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit;
    GraphPlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.n_qubits() == 0);
  }
  GIVEN("Single qubit circuit, two qubit Architecture, GraphPlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(1);
    GraphPlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Qubit(Placement::unplaced_reg(), 0));
  }
  GIVEN(
      "Two qubit unconnected circuit, two qubit Architecture, "
      "GraphPlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(2);
    GraphPlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Qubit(Placement::unplaced_reg(), 0));
    REQUIRE(circuit.all_qubits()[1] == Qubit(Placement::unplaced_reg(), 1));
  }
  GIVEN(
      "Three qubit unconnected circuit, two qubit Architecture, "
      "GraphPlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(3);
    GraphPlacement placement(architecture);
    REQUIRE_THROWS_AS(placement.place(circuit), std::invalid_argument);
  }
  GIVEN(
      "Two qubit connected circuit, three qubit Architecture, "
      "GraphPlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    GraphPlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Node(1));
    REQUIRE(circuit.all_qubits()[1] == Node(2));
  }
  GIVEN(
      "Five qubit connected CX circuit, five qubit Architecture, many relevant "
      "isomorphisms, GraphPlacement::get_all_placement_maps") {
    /**
     * Architecture graph:
     *      4
     *      |
     * 2 -- 0 -- 1
     *      |
     *      3
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1}, {0, 2}, {0, 3}, {0, 4}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph:
     *      4
     *      |
     * 2 -- 0 -- 1
     *      |
     *      3
     */
    Circuit circuit(5);
    add_2qb_gates(circuit, OpType::CX, {{0, 1}, {0, 2}, {0, 3}, {0, 4}});
    GraphPlacement placement(architecture);
    std::vector<std::map<Qubit, Node>> placement_maps =
        placement.get_all_placement_maps(circuit, 25);
    // any permutation of Qubits 1,2,3,4 on Nodes 1,2,3,4 give identical results
    // n.b. this is fewer than 25, as there should only be 24 matches in this
    // case
    REQUIRE(placement_maps.size() == 24);
    for (const std::map<Qubit, Node>& map : placement_maps) {
      REQUIRE(map.at(Qubit(0)) == Node(0));
    }
  }
  GIVEN(
      "Six qubit connected CX circuit, six qubit Architecture, exact "
      "isomorphism, GraphPlacement::get_placement_map") {
    /**
     * Architecture graph:
     * 5    4
     * |    |
     * 2 -- 1 -- 0
     *   \  |
     *      3
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}, {1, 3},
                                                        {1, 4}, {2, 3}, {2, 5}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph:
     * 5    4
     * |    |
     * 2 -- 1 -- 0
     *   \  |
     *      3
     */
    Circuit circuit(6);
    add_2qb_gates(
        circuit, OpType::CX, {{0, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 5}});
    GraphPlacement placement(architecture);
    std::vector<std::map<Qubit, Node>> placement_maps =
        placement.get_all_placement_maps(circuit, 3);
    // 0 and 4 can be swapped without impacting results, giving two maps
    // n.b. this is fewer than 3 as there should only be 2 matches in this case
    REQUIRE(placement_maps.size() == 2);
    std::map<Qubit, Node> placement_map = placement_maps[0];
    REQUIRE(placement_map[Qubit(0)] == Node(4));
    REQUIRE(placement_map[Qubit(1)] == Node(1));
    REQUIRE(placement_map[Qubit(2)] == Node(2));
    REQUIRE(placement_map[Qubit(3)] == Node(3));
    REQUIRE(placement_map[Qubit(4)] == Node(0));
    REQUIRE(placement_map[Qubit(5)] == Node(5));
    REQUIRE(placement_maps[1][Qubit(0)] == Node(0));
    REQUIRE(placement_maps[1][Qubit(4)] == Node(4));
  }

  GIVEN(
      "Nine qubit disconnected CX circuit, 15 qubit Architecture with multiple "
      "mappings, no exact isomorphism, GraphPlacement::get_placement_map, one "
      "unplaced qubit.") {
    /**
     * Architecture graph:
     * 0 -- 1 -- 2 -- 3 -- 4 -- 5
     * |                   |
     * 10                  11
     * |                   |
     * 13-- 14-- 15-- 16-- 17-- 18
     *           |
     *           19
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1},   {1, 2},   {2, 3},   {3, 4},   {4, 5},
        {0, 10},  {10, 13}, {4, 11},  {11, 17}, {13, 14},
        {14, 15}, {15, 16}, {16, 17}, {17, 18}, {15, 19}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph 1:
     * 5 -- 1 -- 3
     */
    /**
     * Qubit Interaction graph 2:
     *           2
     *           |
     * 4 -- 7 -- 0 -- 8
     *           |
     *           6
     */
    Circuit circuit(9);
    add_2qb_gates(
        circuit, OpType::CX,
        {{8, 0}, {5, 1}, {4, 7}, {0, 6}, {1, 3}, {0, 2}, {7, 0}});
    GraphPlacement placement(architecture, 2000, 200000);
    std::map<Qubit, Node> placement_map = placement.get_placement_map(circuit);
    REQUIRE(placement_map[Qubit(0)] == Node(17));
    REQUIRE(placement_map[Qubit(1)] == Node(13));
    REQUIRE(placement_map[Qubit(2)] == Node(18));
    REQUIRE(placement_map[Qubit(3)] == Node(10));
    REQUIRE(placement_map[Qubit(4)] == Node(15));
    REQUIRE(placement_map[Qubit(5)] == Node(14));
    REQUIRE(placement_map[Qubit(6)] == Node("unplaced", 0));
    REQUIRE(placement_map[Qubit(7)] == Node(16));
    REQUIRE(placement_map[Qubit(8)] == Node(11));
  }
  GIVEN(
      "Large Architecture, large complex circuit, small timeout, "
      "Placement::get_all_placement_maps runtime_error.") {
    SquareGrid architecture(20, 20);
    Circuit circuit(10);
    std::vector<std::pair<unsigned, unsigned>> cx_gates;
    for (unsigned i = 0; i < 10; i++) {
      for (unsigned j = 0; j < 10; j++) {
        if (i != j) {
          cx_gates.push_back({i, j});
        }
      }
    }
    add_2qb_gates(circuit, OpType::CX, cx_gates);
    GraphPlacement placement(architecture, 100000, 1);
    REQUIRE_THROWS_AS(
        placement.get_all_placement_maps(circuit, 100000), std::runtime_error);
  }
  GIVEN(
      "9 Qubit Circuit, 9 Qubit Architecture, "
      "GraphPlacement::get_placement_map, increasing number of edges in "
      "pattern graph until subgraph monomorphism found.") {
    /**
     * Architecture graph:
     *      2         6
     *      |         |
     * 0 -- 1 -- 4 -- 5 -- 8
     *      |         |
     *      3         7
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1}, {1, 2}, {1, 3}, {1, 4}, {4, 5}, {5, 6}, {5, 7}, {5, 8}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph:
     *      2         6
     *      |         |
     * 0 -- 1 -- 4 -- 5 -- 8
     *      |         |
     *      3         7
     */
    Circuit circuit(9);
    add_2qb_gates(circuit, OpType::CX, edges);

    // only allow 1 edge in pattern graph, don't find solutions
    GraphPlacement placement(architecture, 1000, 100, 1, 1);
    std::vector<std::map<Qubit, Node>> placement_maps =
        placement.get_all_placement_maps(circuit, 1000);
    REQUIRE(placement_maps.size() == 16);
    std::map<Qubit, Node> map = placement_maps[0];
    REQUIRE(map[Qubit(0)] == Node(2));
    REQUIRE(map[Qubit(1)] == Node(1));
    REQUIRE(map.find(Qubit(2)) == map.end());
    REQUIRE(map.find(Qubit(3)) == map.end());
    REQUIRE(map.find(Qubit(4)) == map.end());
    REQUIRE(map.find(Qubit(5)) == map.end());
    REQUIRE(map.find(Qubit(6)) == map.end());
    REQUIRE(map.find(Qubit(7)) == map.end());
    REQUIRE(map.find(Qubit(8)) == map.end());

    // allow more edges in pattern graph, find better solutions
    placement = GraphPlacement(architecture, 1000, 100, 3, 3);
    placement_maps = placement.get_all_placement_maps(circuit, 1000);
    REQUIRE(placement_maps.size() == 48);
    map = placement_maps[0];
    REQUIRE(map[Qubit(0)] == Node(7));
    REQUIRE(map[Qubit(1)] == Node(5));
    REQUIRE(map[Qubit(2)] == Node(6));
    REQUIRE(map[Qubit(3)] == Node(4));
    REQUIRE(map.find(Qubit(4)) == map.end());
    REQUIRE(map.find(Qubit(5)) == map.end());
    REQUIRE(map.find(Qubit(6)) == map.end());
    REQUIRE(map.find(Qubit(7)) == map.end());
    REQUIRE(map.find(Qubit(8)) == map.end());

    // allow 9 edges in pattern graph, find full solutions
    placement = GraphPlacement(architecture, 1000, 100, 9, 9);
    placement_maps = placement.get_all_placement_maps(circuit, 1000);
    REQUIRE(placement_maps.size() == 72);
    map = placement_maps[0];
    REQUIRE(map[Qubit(0)] == Node(6));
    REQUIRE(map[Qubit(1)] == Node(5));
    REQUIRE(map[Qubit(2)] == Node(8));
    REQUIRE(map[Qubit(3)] == Node(7));
    REQUIRE(map[Qubit(4)] == Node(4));
    REQUIRE(map[Qubit(5)] == Node(1));
    REQUIRE(map[Qubit(6)] == Node(0));
    REQUIRE(map[Qubit(7)] == Node(3));
    REQUIRE(map[Qubit(8)] == Node(2));
    for (std::map<Qubit, Node>& pmap : placement_maps) {
      REQUIRE(pmap[Qubit(4)] == Node(4));
    }
  }
  GIVEN("A Circuit with a Barrier.") {
    Circuit circuit(3, 3);
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    circuit.add_op<unsigned>(OpType::H, {1});
    circuit.add_op<unsigned>(OpType::CX, {1, 2});
    circuit.add_measure(Qubit(0), Bit(0));
    circuit.add_measure(Qubit(1), Bit(1));
    circuit.add_barrier({Qubit(0), Qubit(1), Qubit(2)});
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    circuit.add_op<unsigned>(OpType::H, {0});
    circuit.add_measure(Qubit(2), Bit(2));

    GraphPlacement placement(architecture);
    std::map<Qubit, Node> placement_map = placement.get_placement_map(circuit);
    std::map<Qubit, Node> comparison_map = {
        {Qubit(0), Node(2)}, {Qubit(1), Node(1)}, {Qubit(2), Node(0)}};
    REQUIRE(placement_map == comparison_map);
  }
}

}  // namespace tket
