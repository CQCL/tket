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

#pragma once

#include "PartitionCircuit.hpp"
#include "CeresSolver.hpp"
#include "Gate/Gate.hpp"
#include "Simulation/CircuitSimulator.hpp"

namespace tket {

// maximum distance between original and synthesised circuit
static const double EPSILON = 0.001;

// denotes an edge between two qubits in the architecture
typedef std::pair<unsigned, unsigned> Connection;
typedef std::vector<Connection> ConnectionVec;

/**
 * @brief Optimises a circuit for a given architecture by partitioning it into
 * partitions of size k, resynthesising each partition and stitching the 
 * resynthesised circuits back together.
 * 
 * @param circ circuit to be synthesised
 * @param arch target archetecture to synthesise for
 * @param k maximum partition size
 * @return optimised circuit
 */
Circuit optimise(Circuit &circ, Architecture &arch, unsigned int k);

/**
 * @brief Re-synthesised a circuit using the method from:
 * "Towards Optimal Topology Aware Quantum Circuit Synthesis"
 * The aim of the method is to minimise the number of CX gates.
 * 
 * @param partition the circuit partition to be re-synthesised
 * @param arch the target archetechture to synthesise for
 * @return the optimised partition
 */
Partition synthesise(Partition &partition, Architecture &arch);

/** A stuct representing a node in a tree of all possible circuits
 * as described in page 226 of "Towards Optimal Topology Aware Quantum
 *  Circuit Synthesis" */
struct CircuitNode {
  Circuit circuit;
  double cost_estimate;
  double distance;
  int cx_count;
  Eigen::MatrixXcd unitary;
  const Eigen::MatrixXcd *target;

  CircuitNode(Circuit circuit, double cost_estimate, 
  double distance, int cx_count, Eigen::MatrixXcd unitary, 
  const Eigen::MatrixXcd *target) : circuit(circuit), cost_estimate(cost_estimate),
  distance(distance), cx_count(cx_count), unitary(unitary), target(target) {}
};

/**
 * @brief Initialises the root node (empty circuit) of the circuit tree.
 * 
 * @param target the target unitary to synthesise
 * @return the root node of the tree
 */
CircuitNode init_root_node(const Eigen::MatrixXcd *target);

/**
 * @brief Initialises the successor node of a given node. The successor
 * contains an additional CX followed by a U3 gate on each qubit in the
 * provided connection.
 * 
 * @param node parent node of the successor
 * @param conn connected qubits between which to add a CX gate
 * @return the successor node
 */
CircuitNode init_successor_node(CircuitNode &node, Connection &conn);

/**
 * @brief Creates a vector of the connected qubits in a partition. Also
 * Translates the qubit indices from the original circuit to their new
 * indices in the partition.
 * 
 * @param arch the target archetecture to synthesise for
 * @param qubits the qubits to check for connections
 * @return vector of connected (translated) qubits
 */
ConnectionVec get_connected_qubits(
  Architecture &arch, qubit_vector_t &qubits);

} // namespace tket