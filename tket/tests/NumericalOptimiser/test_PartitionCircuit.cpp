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

#include "NumericalOptimiser/PartitionCircuit.hpp"
#include "Placement/Placement.hpp"

namespace tket {
namespace test_PartitionCircuit {

SCENARIO("Testing partition") {
  GIVEN("A simple circuit to partition") {
    Circuit test_circ(3);
    test_circ.add_op<unsigned>(OpType::H, {0});
    test_circ.add_op<unsigned>(OpType::H, {1});
    test_circ.add_op<unsigned>(OpType::H, {2});
    test_circ.add_op<unsigned>(OpType::CX, {0, 1});
    test_circ.add_op<unsigned>(OpType::CX, {1, 2});
    test_circ.add_op<unsigned>(OpType::H, {0});
    test_circ.add_op<unsigned>(OpType::H, {1});
    test_circ.add_op<unsigned>(OpType::H, {2});

    Circuit partition_0(2);
    partition_0.add_op<unsigned>(OpType::H, {0});
    partition_0.add_op<unsigned>(OpType::H, {1});
    partition_0.add_op<unsigned>(OpType::CX, {0, 1});
    partition_0.add_op<unsigned>(OpType::H, {0});

    Circuit partition_1(2);
    partition_1.add_op<unsigned>(OpType::H, {1});
    partition_1.add_op<unsigned>(OpType::CX, {0, 1});
    partition_1.add_op<unsigned>(OpType::H, {0});
    partition_1.add_op<unsigned>(OpType::H, {1});
    
    Architecture arch({{0, 1}, {1, 2}});
    NaivePlacement np(arch);
    qubit_mapping_t map = np.get_placement_map(test_circ);
    test_circ.rename_units(map);
    
    PartitionVec partitions = partition(test_circ, arch, 2);

    // TODO: this comparison method occasionally fails
    // (seemingly randomly?) find out why
    // REQUIRE(partitions[0].first == partition_0);
    // REQUIRE(partitions[1].first == partition_1);
  }
}

SCENARIO("Testing get_connected_subarch") {
  GIVEN("A 2 by 2 grid architecture with partition size 3") {
    SquareGrid arch(2, 2);
    std::vector<Node> nodes = arch.get_all_nodes_vec();
    node_set_t node_set_1 = {nodes[0], nodes[1], nodes[2]};
    node_set_t node_set_2 = {nodes[0], nodes[1], nodes[3]};
    node_set_t node_set_3 = {nodes[0], nodes[2], nodes[3]};
    node_set_t node_set_4 = {nodes[1], nodes[2], nodes[3]};

    std::vector<node_set_t> result = get_connected_subarch(arch, 3);

    REQUIRE(std::find(result.begin(), result.end(), node_set_1) != result.end());
    REQUIRE(std::find(result.begin(), result.end(), node_set_2) != result.end());
    REQUIRE(std::find(result.begin(), result.end(), node_set_3) != result.end());
    REQUIRE(std::find(result.begin(), result.end(), node_set_4) != result.end());
  }
}

SCENARIO("Testing max_partition") {
  GIVEN("A single gate circuit") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::H, {0});
    qubit_vector_t qubits = circ.all_qubits();
    Subcircuit max_partition = get_max_partition(circ, qubits);
    REQUIRE(max_partition.verts.size() == 1);
    REQUIRE(circ.get_OpType_from_Vertex(*max_partition.verts.begin()) == OpType::H);
  }
  GIVEN("A simple 3 qubit gate with 2 qubit partition") {
    Circuit circ(3);
    Vertex v0 = circ.add_op<unsigned>(OpType::H, {0});
    Vertex v1 = circ.add_op<unsigned>(OpType::H, {2});
    Vertex v2 = circ.add_op<unsigned>(OpType::CX, {0, 1});
    Vertex v3 = circ.add_op<unsigned>(OpType::CX, {1, 2});
    Vertex v4 = circ.add_op<unsigned>(OpType::H, {0});
    Vertex v5 = circ.add_op<unsigned>(OpType::H, {2});
    qubit_vector_t qubits = circ.all_qubits();
    qubit_vector_t valid_qubits = {qubits.begin(), qubits.end() - 1};

    Subcircuit max_partition = get_max_partition(circ, valid_qubits);

    VertexSet vertices = max_partition.verts;
    REQUIRE(max_partition.verts.size() == 3);
    REQUIRE(vertices.find(v0) != vertices.end());
    REQUIRE(vertices.find(v2) != vertices.end());
    REQUIRE(vertices.find(v4) != vertices.end());
  }
}

SCENARIO("Testing get_all_predecessors") {
  GIVEN("A single qubit circuit") {
    VertexSet result;
    Circuit circ(1);
    Vertex child = circ.add_op<unsigned>(OpType::H, {0});
    Vertex root = circ.add_op<unsigned>(OpType::H, {0});

    get_all_predecessors(circ, root, result);

    REQUIRE(result.size() == 1);
    REQUIRE(result.find(child) != result.end());
  }
  GIVEN("A short 3 qubit circuit") {
    VertexSet result;
    Circuit circ(3);
    Vertex v0 = circ.add_op<unsigned>(OpType::H, {0});
    Vertex v1 = circ.add_op<unsigned>(OpType::H, {2});
    Vertex v2 = circ.add_op<unsigned>(OpType::CX, {0, 1});
    Vertex v3 = circ.add_op<unsigned>(OpType::CX, {1, 2});
    Vertex root = circ.add_op<unsigned>(OpType::H, {0});
    Vertex v4 = circ.add_op<unsigned>(OpType::H, {2});

    get_all_predecessors(circ, root, result);

    REQUIRE(result.size() == 2);
    REQUIRE(result.find(v0) != result.end());
    REQUIRE(result.find(v2) != result.end());
  }
}
}  // namespace test_PartitionCircuit
}  // namespace tket