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

SCENARIO("Base NoiseAwarePlacement class") {
    GIVEN("Empty Architecture, NoiseAwarePlacement::NoiseAwarePlacement."){
    Architecture architecture;
    REQUIRE_THROWS_AS(NoiseAwarePlacement(architecture), std::logic_error);
    }
  GIVEN("Empty circuit, two qubit Architecture, NoiseAwarePlacement::Place.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit;
    NoiseAwarePlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.n_qubits() == 0);
  }
  GIVEN("Single qubit circuit, two qubit Architecture, NoiseAwarePlacement::Place.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(1);
    NoiseAwarePlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Qubit(0));
  }
  GIVEN(
      "Two qubit unconnected circuit, two qubit Architecture, "
      "NoiseAwarePlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(2);
    NoiseAwarePlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Qubit(0));
    REQUIRE(circuit.all_qubits()[1] == Qubit(1));
  }
  GIVEN(
      "Three qubit unconnected circuit, two qubit Architecture, "
      "NoiseAwarePlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(3);
    NoiseAwarePlacement placement(architecture);
    REQUIRE_THROWS_AS(placement.place(circuit), std::invalid_argument);
  }
  GIVEN(
      "Two qubit connected circuit, three qubit Architecture, "
      "NoiseAwarePlacement::get_placement_map, noise.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit copy = circuit;
    NoiseAwarePlacement placement(architecture);
    std::map<Qubit, Node> map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(2));
    REQUIRE(map[Qubit(1)] == Node(1));

    gate_error_t double_gate_error_0(0.2), double_gate_error_1(0.8);
    avg_link_errors_t op_link_errors;
    op_link_errors[{Node(0), Node(1)}] = double_gate_error_0;
    op_link_errors[{Node(1), Node(2)}] = double_gate_error_1;
    DeviceCharacterisation characterisation({}, op_link_errors);

    placement.set_characterisation(characterisation);
    map = placement.get_placement_map(copy);

    REQUIRE(map[Qubit(0)] == Node(0));
    REQUIRE(map[Qubit(1)] == Node(1));
  }
  GIVEN("Four qubit connected circuit, 8 qubit Architecture, NoiseAwarePlacement::get_placement_map, noise."){
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}, {2,3}, {3,0}, {1,4}, {4,5}, {5,6}, {6,7}, {2,7}, {4,7}};
    Architecture architecture(edges);
    Circuit circuit(4);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    circuit.add_op<unsigned>(OpType::CX, {1, 2});
    circuit.add_op<unsigned>(OpType::CX, {2, 3});
    circuit.add_op<unsigned>(OpType::CX, {0, 3});
    Circuit copy = circuit;
    
    NoiseAwarePlacement placement(architecture);
    std::map<Qubit,Node> map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(3));
    REQUIRE(map[Qubit(1)] == Node(0));
    REQUIRE(map[Qubit(2)] == Node(1));
    REQUIRE(map[Qubit(3)] == Node(2));


    gate_error_t double_gate_error_1(0.1), double_gate_error_2(0.2), double_gate_error_3(0.3), double_gate_error_4(0.4), double_gate_error_5(0.5), double_gate_error_6(0.6), double_gate_error_7(0.7); 
    avg_link_errors_t op_link_errors;
    op_link_errors[{Node(0), Node(3)}] = double_gate_error_7;
    op_link_errors[{Node(0), Node(1)}] = double_gate_error_6;
    op_link_errors[{Node(2), Node(3)}] = double_gate_error_6;
    op_link_errors[{Node(1), Node(2)}] = double_gate_error_5;
    op_link_errors[{Node(1), Node(4)}] = double_gate_error_4;
    op_link_errors[{Node(2), Node(7)}] = double_gate_error_4;
    op_link_errors[{Node(4), Node(7)}] = double_gate_error_3;
    op_link_errors[{Node(4), Node(5)}] = double_gate_error_2;
    op_link_errors[{Node(7), Node(6)}] = double_gate_error_2;
    op_link_errors[{Node(5), Node(6)}] = double_gate_error_1;
    DeviceCharacterisation characterisation({}, op_link_errors);

    placement.set_characterisation(characterisation);
    map = placement.get_placement_map(copy);
    REQUIRE(map[Qubit(0)] == Node(7));
    REQUIRE(map[Qubit(1)] == Node(6));
    REQUIRE(map[Qubit(2)] == Node(5));
    REQUIRE(map[Qubit(3)] == Node(4));
  }
}
} // namespace tket