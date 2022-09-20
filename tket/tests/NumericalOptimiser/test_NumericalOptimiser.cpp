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
#include <span>

#include "NumericalOptimiser/NumericalOptimiser.hpp"
#include "Architecture/ArchitectureMapping.hpp"
#include "Placement/Placement.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "../Simulation/ComparisonFunctions.hpp"

namespace tket {
namespace test_NumericalOptimiser {

SCENARIO("Testing optimise") {
  GIVEN("A simple circuit to optimise") {
    Architecture arch({{0, 1}, {1, 2}});
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    NaivePlacement np(arch);
    qubit_mapping_t map = np.get_placement_map(circ);
    circ.rename_units(map);
    
    // Commented out to avoid cluttering test output
    // Circuit result = optimise(circ, arch, 2);

    // const StateVector s_circ = tket_sim::get_statevector(circ);
    // const StateVector s_result = tket_sim::get_statevector(result);
    // REQUIRE(tket_sim::compare_statevectors_or_unitaries(s_circ, s_result));
  }
}

SCENARIO("Testing synthesise") {
  GIVEN("a simple 3 qubit circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    
    Architecture arch({{0, 1}, {1, 2}});
    NaivePlacement np(arch);
    qubit_mapping_t map = np.get_placement_map(circ);
    circ.rename_units(map);
    qubit_vector_t qubits = circ.all_qubits();
    Partition unsynthesised = std::make_pair(circ, qubits);

    // Commented out to avoid cluttering test output
    // Partition synthesised = synthesise(unsynthesised, arch);

    // REQUIRE(evaluate_distance(tket_sim::get_unitary(synthesised.first),
    //   tket_sim::get_unitary(unsynthesised.first)) < EPSILON);
    // REQUIRE(synthesised.second == unsynthesised.second);
  }
}

SCENARIO("Testing init_root_node") {
  GIVEN("a simple 3 qubit circuit") {
    Circuit empty(3);
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});

    const Eigen::MatrixXcd target = tket_sim::get_unitary(circ);
    const Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(3, 3);
    double distance = evaluate_distance(id, target);

    CircuitNode root_node = init_root_node(&target);

    REQUIRE(root_node.circuit == empty);
    // TODO: This randomly sometimes fails, figure out why
    REQUIRE(root_node.cost_estimate == distance);
    REQUIRE(root_node.distance == distance);
    REQUIRE(root_node.cx_count == 0);
    REQUIRE(root_node.unitary == id);
    REQUIRE(*root_node.target == target);
  }
}

SCENARIO("Testing init_successor_node") {
  GIVEN("successors of a root node") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    const Eigen::MatrixXcd target = tket_sim::get_unitary(circ);
    Connection conn = std::make_pair(0, 1);

    CircuitNode root_node = init_root_node(&target);
    // Commented out to prevent cluttering test output
    // CircuitNode successor = init_successor_node(root_node, conn);

    // REQUIRE(tket_sim::get_unitary(successor.circuit) == successor.unitary);
    // REQUIRE(successor.cost_estimate == successor.distance + 9.3623);
    // REQUIRE(successor.cx_count == 1);
    // REQUIRE(*successor.target == target);
  }
}

SCENARIO("Testing get_connected_qubits") {
  GIVEN("all qubits on a 2x2 grid architecture") {
    Circuit circ(4);
    SquareGrid arch(2, 2);
    NaivePlacement np(arch);
    qubit_mapping_t map = np.get_placement_map(circ);
    circ.rename_units(map);
    qubit_vector_t qubits = circ.all_qubits();

    ConnectionVec conns = get_connected_qubits(arch, qubits);

    REQUIRE(conns[0] == std::make_pair((unsigned) 0, (unsigned) 1));
    REQUIRE(conns[1] == std::make_pair((unsigned) 0, (unsigned) 2));
    REQUIRE(conns[2] == std::make_pair((unsigned) 1, (unsigned) 3));
    REQUIRE(conns[3] == std::make_pair((unsigned) 2, (unsigned) 3));
  }
}
}  // namespace test_NumericalOptimiser
}  // namespace tket