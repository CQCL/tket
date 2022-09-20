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

SCENARIO("Base GraphPlacement class") {
  GIVEN("Empty Circuit, Empty Architecture, GraphPlacement::Place.") {
    Architecture architecture;
    Circuit circuit;
    GraphPlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.n_qubits() == 0);
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
    REQUIRE(circuit.all_qubits()[0] == Qubit(0));
  }
  GIVEN(
      "Two qubit unconnected circuit, two qubit Architecture, "
      "GraphPlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(2);
    GraphPlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Qubit(0));
    REQUIRE(circuit.all_qubits()[1] == Qubit(1));
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
        placement.get_all_placement_maps(circuit);
    // any permutation of Qubits 1,2,3,4 on Nodes 1,2,3,4 give identical results
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
        placement.get_all_placement_maps(circuit);
    // 0 and 4 can be swapped without impacting results, giving two maps
    // REQUIRE(placement_maps.size() == 2);
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
      "mappings, no exact isomorphism, GraphPlacement::get_placement_map.") {
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
    GraphPlacement placement(architecture);
    std::map<Qubit, Node> placement_map = placement.get_placement_map(circuit);
    std::cout << "Return Placement Map! " << std::endl;
    for (auto x : placement_map) {
      std::cout << x.first.repr() << " " << x.second.repr() << std::endl;
    }
  }
}

}  // namespace tket
