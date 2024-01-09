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

SCENARIO("Base NoiseAwarePlacement class") {
  GIVEN("Empty Architecture, NoiseAwarePlacement::NoiseAwarePlacement.") {
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
  GIVEN(
      "Single qubit circuit, two qubit Architecture, "
      "NoiseAwarePlacement::Place.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(1);
    NoiseAwarePlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
  }
  GIVEN(
      "Two qubit unconnected circuit, two qubit Architecture, "
      "NoiseAwarePlacement::Place") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}};
    Architecture architecture(edges);
    Circuit circuit(2);
    NoiseAwarePlacement placement(architecture);
    placement.place(circuit);
    REQUIRE(circuit.all_qubits()[0] == Node(0));
    REQUIRE(circuit.all_qubits()[1] == Node(1));
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
      "NoiseAwarePlacement::get_placement_map, single qubit noise.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit copy = circuit;
    NoiseAwarePlacement placement(architecture);
    std::map<Qubit, Node> map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(2));
    REQUIRE(map[Qubit(1)] == Node(1));

    gate_error_t single_gate_error_0(0.2), single_gate_error_1(0.3),
        single_gate_error_2(0.5);
    avg_node_errors_t op_node_errors;
    op_node_errors[Node(0)] = single_gate_error_0;
    op_node_errors[Node(1)] = single_gate_error_1;
    op_node_errors[Node(2)] = single_gate_error_2;
    DeviceCharacterisation characterisation(op_node_errors, {});

    placement.set_characterisation(characterisation);
    map = placement.get_placement_map(copy);
    REQUIRE(map[Qubit(0)] == Node(0));
    REQUIRE(map[Qubit(1)] == Node(1));
  }
  GIVEN(
      "Two qubit connected circuit, three qubit Architecture, "
      "NoiseAwarePlacement::get_placement_map, single qubit and readout "
      "noise.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    Circuit circuit(2);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    Circuit copy = circuit;
    NoiseAwarePlacement placement(architecture);
    std::map<Qubit, Node> map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(2));
    REQUIRE(map[Qubit(1)] == Node(1));

    gate_error_t single_gate_error_0(0.2);
    avg_node_errors_t op_node_errors;
    op_node_errors[Node(0)] = single_gate_error_0;
    op_node_errors[Node(1)] = single_gate_error_0;
    op_node_errors[Node(2)] = single_gate_error_0;

    readout_error_t single_readout_error_0(0.2), single_readout_error_1(0.3),
        single_readout_error_2(0.5);
    avg_readout_errors_t readout_node_errors;
    readout_node_errors[Node(0)] = single_readout_error_0;
    readout_node_errors[Node(1)] = single_readout_error_1;
    readout_node_errors[Node(2)] = single_readout_error_2;
    DeviceCharacterisation characterisation(
        op_node_errors, {}, readout_node_errors);

    placement.set_characterisation(characterisation);
    map = placement.get_placement_map(copy);
    REQUIRE(map[Qubit(0)] == Node(0));
    REQUIRE(map[Qubit(1)] == Node(1));
  }
  GIVEN(
      "Two qubit connected circuit, three qubit Architecture, "
      "NoiseAwarePlacement::get_placement_map, two qubit noise.") {
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
  GIVEN(
      "Four qubit connected circuit, 8 qubit Architecture, "
      "NoiseAwarePlacement::get_placement_map, unhomogeneous two qubit "
      "noise.") {
    /**
     * Architecture graph:
     * 0 -- 1 -- 4 -- 5
     * |    |    |    |
     * 3 -- 2 -- 7 -- 6
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 4},
        {4, 5}, {5, 6}, {6, 7}, {2, 7}, {4, 7}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph:
     * 0 -- 1
     * |    |
     * 2 -- 3
     */
    Circuit circuit(4);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    circuit.add_op<unsigned>(OpType::CX, {1, 2});
    circuit.add_op<unsigned>(OpType::CX, {2, 3});
    circuit.add_op<unsigned>(OpType::CX, {0, 3});
    Circuit copy = circuit;

    NoiseAwarePlacement placement(architecture);
    std::map<Qubit, Node> map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(3));
    REQUIRE(map[Qubit(1)] == Node(0));
    REQUIRE(map[Qubit(2)] == Node(1));
    REQUIRE(map[Qubit(3)] == Node(2));

    gate_error_t double_gate_error_1(0.1), double_gate_error_2(0.2),
        double_gate_error_3(0.3), double_gate_error_4(0.4),
        double_gate_error_5(0.5), double_gate_error_6(0.6),
        double_gate_error_7(0.7);
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
  GIVEN(
      "Four qubit connected circuit, 8 qubit Architecture, "
      "NoiseAwarePlacement::get_placement_map, homogeneous two qubit noise, "
      "unhomogeneous single qubit noise, unhomogeneous readout noise.") {
    /**
     * Architecture graph:
     * 0 -- 1 -- 4 -- 5
     * |    |    |    |
     * 3 -- 2 -- 7 -- 6
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 4},
        {4, 5}, {5, 6}, {6, 7}, {2, 7}, {4, 7}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph:
     * 0 -- 1
     * |    |
     * 2 -- 3
     */
    Circuit circuit(4);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    circuit.add_op<unsigned>(OpType::CX, {1, 2});
    circuit.add_op<unsigned>(OpType::CX, {2, 3});
    circuit.add_op<unsigned>(OpType::CX, {0, 3});
    // In this case there are many valid placements, it happens to return
    // this one
    NoiseAwarePlacement placement(architecture);
    std::map<Qubit, Node> map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(3));
    REQUIRE(map[Qubit(1)] == Node(0));
    REQUIRE(map[Qubit(2)] == Node(1));
    REQUIRE(map[Qubit(3)] == Node(2));

    gate_error_t double_gate_error(0.1);
    avg_link_errors_t op_link_errors;
    op_link_errors[{Node(0), Node(3)}] = double_gate_error;
    op_link_errors[{Node(0), Node(1)}] = double_gate_error;
    op_link_errors[{Node(2), Node(3)}] = double_gate_error;
    op_link_errors[{Node(1), Node(2)}] = double_gate_error;
    op_link_errors[{Node(1), Node(4)}] = double_gate_error;
    op_link_errors[{Node(2), Node(7)}] = double_gate_error;
    op_link_errors[{Node(4), Node(7)}] = double_gate_error;
    op_link_errors[{Node(4), Node(5)}] = double_gate_error;
    op_link_errors[{Node(7), Node(6)}] = double_gate_error;
    op_link_errors[{Node(5), Node(6)}] = double_gate_error;
    DeviceCharacterisation characterisation_link({}, op_link_errors);

    // similarly as all gate errors are identical, all maps are valid
    placement.set_characterisation(characterisation_link);
    map = placement.get_placement_map(circuit);
    REQUIRE(map[Qubit(0)] == Node(2));
    REQUIRE(map[Qubit(1)] == Node(7));
    REQUIRE(map[Qubit(2)] == Node(4));
    REQUIRE(map[Qubit(3)] == Node(1));

    gate_error_t single_gate_error_0(0.3), single_gate_error_1(0.4);
    avg_node_errors_t op_node_errors;
    op_node_errors[Node(0)] = single_gate_error_1;
    op_node_errors[Node(1)] = single_gate_error_1;
    op_node_errors[Node(2)] = single_gate_error_1;
    op_node_errors[Node(3)] = single_gate_error_1;
    op_node_errors[Node(4)] = single_gate_error_0;
    op_node_errors[Node(5)] = single_gate_error_0;
    op_node_errors[Node(6)] = single_gate_error_0;
    op_node_errors[Node(7)] = single_gate_error_0;

    DeviceCharacterisation characterisation_link_node(
        op_node_errors, op_link_errors);

    // Here the difference in single qubit error rates makes this placement
    // (or a rotation of) best
    placement.set_characterisation(characterisation_link_node);
    map = placement.get_placement_map(circuit);

    REQUIRE(map[Qubit(0)] == Node(7));
    REQUIRE(map[Qubit(1)] == Node(6));
    REQUIRE(map[Qubit(2)] == Node(5));
    REQUIRE(map[Qubit(3)] == Node(4));

    op_node_errors[Node(0)] = single_gate_error_0;
    op_node_errors[Node(1)] = single_gate_error_0;
    op_node_errors[Node(2)] = single_gate_error_0;
    op_node_errors[Node(3)] = single_gate_error_0;

    readout_error_t single_readout_error_0(0.05), single_readout_error_1(0.9);
    avg_readout_errors_t readout_node_errors;
    readout_node_errors[Node(0)] = single_readout_error_1;
    readout_node_errors[Node(1)] = single_readout_error_0;
    readout_node_errors[Node(2)] = single_readout_error_0;
    readout_node_errors[Node(3)] = single_readout_error_1;
    readout_node_errors[Node(4)] = single_readout_error_0;
    readout_node_errors[Node(5)] = single_readout_error_1;
    readout_node_errors[Node(6)] = single_readout_error_1;
    readout_node_errors[Node(7)] = single_readout_error_0;

    DeviceCharacterisation characterisation_link_node_readout(
        op_node_errors, op_link_errors, readout_node_errors);

    // Here the readout errors are more potent than the single qubit errors, so
    // it now assigns to a different qubit subset
    placement.set_characterisation(characterisation_link_node_readout);
    map = placement.get_placement_map(circuit);

    REQUIRE(map[Qubit(0)] == Node(2));
    REQUIRE(map[Qubit(1)] == Node(7));
    REQUIRE(map[Qubit(2)] == Node(4));
    REQUIRE(map[Qubit(3)] == Node(1));
  }
  GIVEN(
      "Six qubit connected circuit, 10 qubit Architecture, "
      "NoiseAwarePlacement::get_all_placement_maps, homogeneous two qubit "
      "noise, "
      "unhomogeneous single qubit noise") {
    /**
     * Architecture Graph:
     *   1 -- 2
     *  /      \
     * 0        3
     *  \      /
     *   5 -- 4
     *  /      \
     * 6        9
     *  \      /
     *   7 -- 8
     */
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {0, 5},
        {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 4}};
    Architecture architecture(edges);
    /**
     * Qubit interaction graph:
     *   1 -- 2
     *  /      \
     * 0        3
     *  \      /
     *   5 -- 4
     */
    Circuit circuit(6);
    circuit.add_op<unsigned>(OpType::CX, {1, 0});
    circuit.add_op<unsigned>(OpType::CX, {2, 3});
    circuit.add_op<unsigned>(OpType::CX, {4, 5});
    circuit.add_op<unsigned>(OpType::CX, {1, 2});
    circuit.add_op<unsigned>(OpType::CX, {4, 3});
    circuit.add_op<unsigned>(OpType::CX, {0, 5});

    NoiseAwarePlacement placement(architecture);
    // note we allow for more matches than should be returned
    // as noise aware placement returns equal best weighted results
    // as it does additional cost with device characteristics
    std::vector<std::map<Qubit, Node>> maps =
        placement.get_all_placement_maps(circuit, 100);
    REQUIRE(maps.size() == 24);

    // now add noise, making the upper hexagon better, such that it returns
    // fewer maps
    gate_error_t double_gate_error_0(0.2), double_gate_error_1(0.3);
    avg_link_errors_t op_link_errors;
    op_link_errors[{Node(0), Node(1)}] = double_gate_error_0;
    op_link_errors[{Node(2), Node(1)}] = double_gate_error_0;
    op_link_errors[{Node(2), Node(3)}] = double_gate_error_0;
    op_link_errors[{Node(4), Node(3)}] = double_gate_error_0;
    op_link_errors[{Node(0), Node(5)}] = double_gate_error_0;
    op_link_errors[{Node(4), Node(5)}] = double_gate_error_0;
    op_link_errors[{Node(4), Node(9)}] = double_gate_error_1;
    op_link_errors[{Node(9), Node(8)}] = double_gate_error_1;
    op_link_errors[{Node(7), Node(8)}] = double_gate_error_1;
    op_link_errors[{Node(6), Node(7)}] = double_gate_error_1;
    op_link_errors[{Node(5), Node(6)}] = double_gate_error_1;

    DeviceCharacterisation characterisation_link({}, op_link_errors);
    placement.set_characterisation(characterisation_link);

    maps = placement.get_all_placement_maps(circuit, 100);
    REQUIRE(maps.size() == 6);
    // there are 6 returned maps as the direction is considered
    // we check one to confirm it's the correct orientation and
    // side of the hexagon, but assume others are suitable ortations.
    std::map<Qubit, Node> map = maps[0];
    REQUIRE(map[Qubit(0)] == Node(3));
    REQUIRE(map[Qubit(1)] == Node(2));
    REQUIRE(map[Qubit(2)] == Node(1));
    REQUIRE(map[Qubit(3)] == Node(0));
    REQUIRE(map[Qubit(4)] == Node(5));
    REQUIRE(map[Qubit(5)] == Node(4));

    // now make the middle segment better, such that there is a single best map
    gate_error_t double_gate_error_2(0.05);
    op_link_errors[{Node(4), Node(5)}] = double_gate_error_2;
    DeviceCharacterisation characterisation_link_middle({}, op_link_errors);
    placement.set_characterisation(characterisation_link_middle);

    maps = placement.get_all_placement_maps(circuit, 100);
    REQUIRE(maps.size() == 1);
    map = maps[0];
    REQUIRE(map[Qubit(0)] == Node(3));
    REQUIRE(map[Qubit(1)] == Node(4));
    REQUIRE(map[Qubit(2)] == Node(5));
    REQUIRE(map[Qubit(3)] == Node(0));
    REQUIRE(map[Qubit(4)] == Node(1));
    REQUIRE(map[Qubit(5)] == Node(2));
  }
  GIVEN(
      "A circuit with only single-qubit gates, assigns Qubits to Nodes with "
      "lowest single qubit error rates.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1}, {1, 2}, {0, 2}, {2, 3}};
    Architecture architecture(edges);
    Circuit circuit(3);
    circuit.add_op<unsigned>(OpType::H, {0});
    circuit.add_op<unsigned>(OpType::H, {1});
    circuit.add_op<unsigned>(OpType::H, {2});

    avg_node_errors_t op_node_errors;
    op_node_errors[Node(0)] = 0.25;
    op_node_errors[Node(1)] = 0.01;
    op_node_errors[Node(2)] = 0.01;
    op_node_errors[Node(3)] = 0.05;

    DeviceCharacterisation characterisation(op_node_errors, {}, {});
    NoiseAwarePlacement placement(architecture);
    placement.set_characterisation(characterisation);

    std::map<Qubit, Node> placement_map = placement.get_placement_map(circuit);
    REQUIRE(placement_map[Qubit(0)] == Node(1));
    REQUIRE(placement_map[Qubit(1)] == Node(2));
    REQUIRE(placement_map[Qubit(2)] == Node(3));
  }
}
}  // namespace tket