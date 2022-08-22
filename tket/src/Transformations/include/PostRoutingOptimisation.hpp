#pragma once

#include <vector>
#include "Circuit/Circuit.hpp"
#include "Architecture/Architecture.hpp"

namespace tket {

typedef std::pair<Circuit, qubit_vector_t> Partition;
typedef std::vector<Partition> PartitionVec;

Circuit optimise(Circuit &circ, Architecture &arch, unsigned int k);

PartitionVec partition(Circuit &circ, Architecture &arch, unsigned int k);

std::vector<node_set_t> get_connected_subarch(Architecture &arch, unsigned int k);

bool expand(
  node_set_t &current, node_set_t &to_expand, node_set_t &to_ignore, 
  Architecture &arch, unsigned int k, std::vector<node_set_t> &result);

void get_all_predecessors(Circuit &circ, Vertex &vertex, VertexSet &result);

Subcircuit get_max_partition(Circuit &circ, qubit_vector_t &qubits);

Partition synthesise(Partition &partition);

} // namespace tket