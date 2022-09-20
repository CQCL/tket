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

#include "NumericalOptimiser.hpp"

namespace tket {

Circuit optimise(Circuit &circ, Architecture &arch, unsigned int k) {
  PartitionVec post_synthesis;
  PartitionVec pre_synthesis = partition(circ, arch, k);

  for (Partition partition : pre_synthesis) {
    post_synthesis.insert(post_synthesis.begin(), synthesise(partition, arch));
  }

  for (Partition partition : post_synthesis) {
    EdgeVec edges = {};
    for (Qubit qubit : partition.second) {
      edges.push_back(circ.get_nth_out_edge(circ.get_in(qubit), 0));
    }
    Subcircuit to_replace = {edges, edges};
    circ.substitute(partition.first, to_replace);
  }

  return circ;
}

Partition synthesise(Partition &partition, Architecture &arch) {
  auto compare_node = [](CircuitNode lhs, CircuitNode rhs) { 
    return lhs.cost_estimate > rhs.cost_estimate;
  };
  std::priority_queue<CircuitNode, std::vector<CircuitNode>, 
    decltype(compare_node)> queue(compare_node);

  Eigen::MatrixXcd target = tket_sim::get_unitary(partition.first);
  queue.push(init_root_node(&target));
  
  while(queue.top().distance > EPSILON) {
    CircuitNode current = queue.top();
    for (Connection conn : get_connected_qubits(arch, partition.second)) {
      queue.push(init_successor_node(current, conn));
    }
  }

  return make_pair(queue.top().circuit, partition.second);
}

CircuitNode init_root_node(const Eigen::MatrixXcd *target) {
  int n = log2(target->cols());
  Circuit circuit(n);
  Eigen::MatrixXcd unitary = Eigen::MatrixXcd::Identity(n, n);
  double distance = evaluate_distance(unitary, *target);
  const CircuitNode node {circuit, distance, distance, 0, unitary, target};
  return node;
}

CircuitNode init_successor_node(CircuitNode &node, Connection &conn) {
  Circuit circuit(node.circuit);
  unsigned int index_1 = conn.first;
  unsigned int index_2 = conn.second;
  
  std::vector<double> p1 = optimise_u3(index_1, node.unitary, *node.target);
  std::vector<double> p2 = optimise_u3(index_2, node.unitary, *node.target);
  
  // TODO: make this work instead of tket_sim::get_unitary
  // Figure out how to make a Gate out of CX (wrong number of parameters)
  // Eigen::MatrixXcd cx = Gate(OpType::CX, {}, {index_1, index_2}).get_unitary();
  // Eigen::MatrixXcd u3_1 = Gate(OpType::U3, {p1[0], p1[1], p1[2]}, index_1).get_unitary();
  // Eigen::MatrixXcd u3_2 = Gate(OpType::U3, {p2[0], p2[1], p2[2]}, index_2).get_unitary();
  
  circuit.add_op<unsigned>(OpType::CX, {index_1, index_2});
  circuit.add_op<unsigned>(OpType::U3, {p1[0], p1[1], p1[2]}, {index_1});
  circuit.add_op<unsigned>(OpType::U3, {p2[0], p2[1], p2[2]}, {index_2});

  Eigen::MatrixXcd unitary = tket_sim::get_unitary(circuit);
  double distance = evaluate_distance(unitary, *node.target);

  const CircuitNode successor {circuit, distance + (node.cx_count + 1) * 9.3623, 
    distance, node.cx_count+1, unitary, node.target};
  
  return successor;
}

ConnectionVec get_connected_qubits(
  Architecture &arch, qubit_vector_t &qubits) {
  ConnectionVec connections;
  for (unsigned i = 0; i < qubits.size(); i++) {
    for (unsigned j = i+1; j < qubits.size(); j++) {
      if (arch.edge_exists(Node(qubits[i]), Node(qubits[j])) || 
        arch.edge_exists(Node(qubits[j]), Node(qubits[i]))) {
        connections.push_back(std::make_pair(i, j));
      }
    }
  }
  return connections;
}
} // namespace tket