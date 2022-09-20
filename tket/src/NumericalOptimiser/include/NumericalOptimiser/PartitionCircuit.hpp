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

#include <vector>

#include "Circuit/Circuit.hpp"
#include "Architecture/Architecture.hpp"

namespace tket {

/** A partition is defined as the circuit it implements paired with an
  * ordered vector of the qubits from the original circuit it acts on */
typedef std::pair<Circuit, qubit_vector_t> Partition;
typedef std::vector<Partition> PartitionVec;

/**
 * @brief Partitions a circuit into subcircuits of maximum size k.
 * The qubits in each partition form a connected subgraph of the given
 * architecture to maximise the number of positions where a CX gate may
 * be placed during resynthesis.
 * 
 * @param circ the circuit to be partitioned
 * @param arch the target architecture the circuit will be run on
 * @param k the maximum number of qubits per partition
 * @return an ordered vector of the partitions
 */
PartitionVec partition(Circuit &circ, Architecture &arch, unsigned int k);

/**
 * @brief Enumerates all connected sub-architectures of order k using
 * the VSimple algorithm from arXiv:2112.07197
 * 
 * @param arch the architecture defining the graph connectivity
 * @param k the desired order of the sub-graphs
 * @return a vector containing all sets of nodes which form a valid sub-arch
 */
std::vector<node_set_t> get_connected_subarch(Architecture &arch, unsigned int k);

/* Auxilliary function for the VSimple algorithm */
bool expand(
  node_set_t &current, node_set_t &to_expand, node_set_t &to_ignore, 
  Architecture &arch, unsigned int k, std::vector<node_set_t> &result);

/**
 * @brief Recursively identifies all vertives which a given vertex depends on.
 * 
 * @param circ the circuit containing the vertex
 * @param vertex the chosen vertex
 * @param result accumulates dependent vertices
 */
void get_all_predecessors(Circuit &circ, Vertex &vertex, VertexSet &result);

/**
 * @brief Finds a subcircuit containing operations which are only 
 * dependent on the given set of qubit. Currently identifies the largest
 * possible subcircuit which starts at the beginning of the circuit.
 * (TODO: larger subcircuits may exist in the middle of the circuit,
 * potentially look into methods to identify these)
 * 
 * @param circ the circuit to partition
 * @param qubits the qubits which the subcircuit is allowed to depend on
 * @return the largest subcircuit
 */
Subcircuit get_max_partition(Circuit &circ, qubit_vector_t &qubits);

} // namespace tket