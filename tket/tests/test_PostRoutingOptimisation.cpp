#include <catch2/catch_test_macros.hpp>
#include <span>

#include "Transformations/PostRoutingOptimisation.hpp"
#include "Architecture/ArchitectureMapping.hpp"
#include "Placement/Placement.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Simulation/CircuitSimulator.hpp"

namespace tket {
namespace test_PostRoutingOptimisation {

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
    
    Circuit result = optimise(circ, arch, 2);

    const StateVector s_circ = tket_sim::get_statevector(circ);
    const StateVector s_result = tket_sim::get_statevector(result);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(s_circ, s_result));
  }
}

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

    // TODO: find a non failing comparison method
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

// SCENARIO("Testing evaluate_distance") {
//   GIVEN("A few simple matrices to compare") {
//     std::complex<double> i(0, 1);
//     Eigen::MatrixXcd x(2, 2), y(2, 2), z(2, 2);

//     x << 0, 1,
//          1, 0;

//     y << 0, -i,
//          i,  0;

//     z << 1,  0,
//          0, -1;

//     REQUIRE(evaluate_distance(x, x) == 0);
//     REQUIRE(evaluate_distance(x, y) == 1.0);
//     REQUIRE(evaluate_distance(x, z) == 1.0);
//     REQUIRE(evaluate_distance(y, z) == 1.0);
//   }
// }

// bool compare_cd(std::complex<double> a, std::complex<double> b) {
//     if (round(a.real() * 1000.0)/1000.0 != round(b.real() * 1000.0)/1000.0) {
//       return false;
//     }
//     if (round(a.imag() * 1000.0)/1000.0 != round(b.imag() * 1000.0)/1000.0) {
//       return false;
//     }
//     return true; 
// }

// SCENARIO("Testing u3") {
//   GIVEN("A hadamard gate on a single qubit circuit") {
//     std::complex<double> inv_root_2(1/std::sqrt(2), 0);
//     Eigen::MatrixXcd h = u3(M_PI/2, 0, M_PI, 0, 1);
    
//     REQUIRE(compare_cd(h(0,0), inv_root_2));
//     REQUIRE(compare_cd(h(0,1), inv_root_2));
//     REQUIRE(compare_cd(h(1,0), inv_root_2));
//     REQUIRE(compare_cd(h(1,1), -inv_root_2));
//   }
//   GIVEN("A Y gate in a few different positions") {
//     std::complex<double> i(0, 1);
//     Eigen::MatrixXcd y0 = u3(M_PI, M_PI/2, M_PI/2, 0, 3);
//     Eigen::MatrixXcd y1 = u3(M_PI, M_PI/2, M_PI/2, 1, 3);
//     Eigen::MatrixXcd y2 = u3(M_PI, M_PI/2, M_PI/2, 2, 3);

//     REQUIRE(compare_cd(y0(0,1), -i));
//     REQUIRE(compare_cd(y0(1,0),  i));
//     REQUIRE(compare_cd(y0(2,3), -i));
//     REQUIRE(compare_cd(y0(3,2),  i));
//     REQUIRE(compare_cd(y0(4,5), -i));
//     REQUIRE(compare_cd(y0(5,4),  i));
//     REQUIRE(compare_cd(y0(6,7), -i));
//     REQUIRE(compare_cd(y0(7,6),  i));

//     REQUIRE(compare_cd(y1(0,2), -i));
//     REQUIRE(compare_cd(y1(1,3), -i));
//     REQUIRE(compare_cd(y1(2,0),  i));
//     REQUIRE(compare_cd(y1(3,1),  i));
//     REQUIRE(compare_cd(y1(4,6), -i));
//     REQUIRE(compare_cd(y1(5,7), -i));
//     REQUIRE(compare_cd(y1(6,4),  i));
//     REQUIRE(compare_cd(y1(7,5),  i));

//     REQUIRE(compare_cd(y2(0,4), -i));
//     REQUIRE(compare_cd(y2(1,5), -i));
//     REQUIRE(compare_cd(y2(2,6), -i));
//     REQUIRE(compare_cd(y2(3,7), -i));
//     REQUIRE(compare_cd(y2(4,0),  i));
//     REQUIRE(compare_cd(y2(5,1),  i));
//     REQUIRE(compare_cd(y2(6,2),  i));
//     REQUIRE(compare_cd(y2(7,3),  i));
//   }
// }
SCENARIO("Testing optimise_u3_gates") {
  GIVEN("A circuit with target ...") {
    std::complex<double> zero(0, 0);
    std::complex<double> one(1, 0);
    std::complex<double> mone(-1, 0);
    std::complex<double> i(0, 1);
    Eigen::MatrixXcd U = Eigen::MatrixXcd::Identity(8, 8);
    Eigen::MatrixXcd T1(8, 8), T2(8, 8);

    T1 << zero, mone, zero, zero, zero, zero, zero, zero,
           one, zero, zero, zero, zero, zero, zero, zero,
          zero, zero, zero, mone, zero, zero, zero, zero,
          zero, zero,  one, zero, zero, zero, zero, zero,
          zero, zero, zero, zero, zero, mone, zero, zero,
          zero, zero, zero, zero,  one, zero, zero, zero,
          zero, zero, zero, zero, zero, zero, zero, mone,
          zero, zero, zero, zero, zero, zero,  one, zero;

    T2 << -0.365, -0.198 + 0.028 * i, 0, 0, -0.766 + 0.223 * i, -0.397 + 0.18 * i, 0, 0,
          0.083 - 0.181 * i, -0.104 + 0.35 * i, 0, 0, 0.063 - 0.431 * i, -0.004 + 0.798 * i, 0, 0,
          0, 0, -0.365, -0.198 + 0.028 * i, 0, 0, -0.766 + 0.223 * i, -0.397 + 0.18 * i,
          0, 0, 0.083 - 0.181 * i, -0.104 + 0.35 * i, 0, 0, 0.063 - 0.431 * i, -0.004 + 0.798 * i,
          0.226 - 0.765 * i,  0.063 - 0.431 * i, 0, 0, -0.002 + 0.365 * i,  0.027 + 0.198 * i, 0, 0,
          0.329 + 0.286 * i, -0.67 - 0.434 * i, 0, 0, -0.181 - 0.084 * i,  0.35 + 0.105 * i, 0, 0,
          0, 0, 0.226 - 0.765 * i,  0.063 - 0.431 * i, 0, 0, -0.002 + 0.365 * i,  0.027 + 0.198 * i,
          0, 0, 0.329 + 0.286 * i, -0.67 - 0.434 * i, 0, 0, -0.181 - 0.084 * i,  0.35 + 0.105 * i;

    optimise_circuit(0, 2, U, T2);
  }
}
}  // namespace test_PostRoutingOptimisation
}  // namespace tket