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

#include "tket/Placement/Placement.hpp"

namespace tket {
SCENARIO("Base Placement class") {
  GIVEN("Empty Circuit, Empty Architecture, Placement::Place.") {
    Architecture architecture;
    Circuit circuit;
    Placement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.n_qubits() == 0);
  }
  GIVEN("Empty circuit, two qubit Architecture, Placement::Place.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit;
    Placement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.n_qubits() == 0);
  }
  GIVEN("Single qubit circuit, two qubit Architecture, Placement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(1);
    Placement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
  }
  GIVEN(
      "Two qubit unconnected circuit, two qubit Architecture, "
      "Placement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(2);
    Placement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
    REQUIRE(circuit.all_qubits()[1] == Node(1));
  }
  GIVEN(
      "Three qubit unconnected circuit, two qubit Architecture, "
      "Placement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(3);
    Placement placement(architecture);
    REQUIRE_THROWS_AS(placement.place(circuit), std::invalid_argument);
  }
  GIVEN(
      "Three qubit unconnected circuit, two qubit Architecture, "
      "Placement::get_all_placement_maps") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(3);
    Placement placement(architecture);
    REQUIRE_THROWS_AS(
        placement.get_all_placement_maps(circuit, 100), std::invalid_argument);
  }
  GIVEN(
      "Two qubit connected circuit, three qubit Architecture, "
      "Placement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    Placement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
    REQUIRE(circuit.all_qubits()[1] == Node(1));
  }
  GIVEN(
      "Two qubit circuit, three qubit Architecture, "
      "Placement::place_with_map, valid map.") {
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    std::map<Qubit, Node> placement_map = {
        {Qubit(0), Node(2)}, {Qubit(1), Node(0)}};
    Placement::place_with_map(circuit, placement_map);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
    REQUIRE(circuit.all_qubits()[1] == Node(2));
  }
  GIVEN(
      "Two qubit circuit, three qubit Architecture, "
      "Placement::place_with_map, valid map with extra Qubit.") {
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    std::map<Qubit, Node> placement_map = {
        {Qubit(0), Node(2)}, {Qubit(1), Node(0)}, {Qubit(3), Node(1)}};
    Placement::place_with_map(circuit, placement_map);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
    REQUIRE(circuit.all_qubits()[1] == Node(2));
  }
  GIVEN(
      "Seven qubit unconnected circuit, seven Qubit architecture, no Qubits "
      "pre-assigned, Placement::get_placement_map") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}, {2, 3},
                                                        {3, 4}, {4, 5}, {5, 6}};
    Architecture architecture(edges);
    Circuit circuit(7);
    Placement placement(architecture);
    std::map<Qubit, Node> placement_map = placement.get_placement_map(circuit);
    REQUIRE(placement_map[Qubit(0)] == Node(0));
    REQUIRE(placement_map[Qubit(1)] == Node(1));
    REQUIRE(placement_map[Qubit(2)] == Node(2));
    REQUIRE(placement_map[Qubit(3)] == Node(3));
    REQUIRE(placement_map[Qubit(4)] == Node(4));
    REQUIRE(placement_map[Qubit(5)] == Node(5));
    REQUIRE(placement_map[Qubit(6)] == Node(6));
  }
  GIVEN(
      "Seven qubit unconnected circuit, seven qubit Architecture, some Qubits "
      "pre-assigned, Placement::get_placement_map and "
      "Placment::place_with_map.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}, {2, 3},
                                                        {3, 4}, {4, 5}, {5, 6}};
    Architecture architecture(edges);
    Circuit circuit(4);
    circuit.add_qubit(Node(0));
    circuit.add_qubit(Node(1));
    circuit.add_qubit(Node(2));
    Placement placement(architecture);
    std::map<Qubit, Node> placement_map = placement.get_placement_map(circuit);
    REQUIRE(placement_map[Qubit(0)] == Node(3));
    REQUIRE(placement_map[Qubit(1)] == Node(4));
    REQUIRE(placement_map[Qubit(2)] == Node(5));
    REQUIRE(placement_map[Qubit(3)] == Node(6));
    REQUIRE(placement_map[Node(0)] == Node(0));
    REQUIRE(placement_map[Node(1)] == Node(1));
    REQUIRE(placement_map[Node(2)] == Node(2));
    Placement::place_with_map(circuit, placement_map);
    std::vector<Qubit> comparison = {Node(0), Node(1), Node(2), Node(3),
                                     Node(4), Node(5), Node(6)};
    REQUIRE(circuit.all_qubits() == comparison);
  }
}
}  // namespace tket
