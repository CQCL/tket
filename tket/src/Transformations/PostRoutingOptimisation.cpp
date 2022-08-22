#include "PostRoutingOptimisation.hpp"
#include <algorithm>

namespace tket {

// Make into a class and store circuit, architecture and partition size as members?

Circuit optimise(Circuit &circ, Architecture &arch, unsigned int k) {
  PartitionVec post_synthesis;
  // Partition circuit
  PartitionVec pre_synthesis = partition(circ, arch, k);
  // Synthesise partitions
  for (Partition partition : pre_synthesis) {
    post_synthesis.insert(post_synthesis.begin(), synthesise(partition));
  }
  // Define empty circuit at the beginning of the circuit to replace with the
  // newly synthesised circuit
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

PartitionVec partition(Circuit &circ, Architecture &arch, unsigned int k) {
  PartitionVec partitions;
  while (circ.n_gates() != 0) {
    Partition max_partition = std::make_pair(Circuit(0), qubit_vector_t());
    Subcircuit max_subcircuit;
    for (node_set_t nodes : get_connected_subarch(arch, k)) {
      qubit_vector_t qubits(nodes.begin(), nodes.end());
      Subcircuit partition = get_max_partition(circ, qubits);
      if (partition.verts.size() > max_partition.first.n_gates()) {
        max_partition = make_pair(circ.subcircuit(partition), qubits);
        max_subcircuit = partition;
      }
    }
    partitions.push_back(max_partition);
    for (Vertex vertex : max_subcircuit.verts) {
      circ.remove_vertex(
        vertex, Circuit::GraphRewiring::Yes,
        Circuit::VertexDeletion::Yes);
    }
  }
  return partitions;
}

// arXiv:2112.07197 : VSimple algorithm for enumerating connected subgraphs of order k
std::vector<node_set_t> get_connected_subarch(Architecture &arch, unsigned int k) {
  std::vector<node_set_t> result;
  node_set_t to_ignore;
  for (Node node : arch.get_all_nodes_vec()) {
    node_set_t current, to_expand;
    current.insert(node);
    node_set_t neighbours = arch.get_neighbour_nodes(node);
    std::set_difference(
      neighbours.begin(), neighbours.end(), 
      to_ignore.begin(), to_ignore.end(), std::inserter(to_expand, to_expand.begin()));
    node_set_t new_to_ignore(to_ignore);
    expand(current, to_expand, to_ignore, arch, k, result);
    to_ignore = new_to_ignore;
    to_ignore.insert(node);
  }
  return result;
}

bool expand(
  node_set_t &current, node_set_t &to_expand, node_set_t &to_ignore, 
  Architecture &arch, unsigned int k, std::vector<node_set_t> &result) {
  // add current node group to result if it's the correct size
  if (current.size() == k) {
    result.push_back(current);
    return true;
  }
  // recursively expand tree of neighbouring nodes to find connected groups
  bool is_done = false;
  for (Node node : to_expand) {
    node_set_t new_current(current);
    new_current.insert(node);
    node_set_t new_to_expand(to_expand);
    node_set_t neighbours = arch.get_neighbour_nodes(node);
    new_to_expand.insert(neighbours.begin(), neighbours.end());
    node_set_t set1, set2;
    std::set_difference(
      new_to_expand.begin(), new_to_expand.end(), 
      new_current.begin(), new_current.end(), std::inserter(set1, set1.begin()));
    std::set_difference(
      set1.begin(), set1.end(),
      to_ignore.begin(), to_ignore.end(), std::inserter(set2, set2.begin()));
    node_set_t new_to_ignore(to_ignore);
    if (expand(new_current, set2, new_to_ignore, arch, k, result)) {
      is_done = true;
    } else { break; }
    to_ignore.insert(node);
    if (arch.n_nodes() - to_ignore.size() < k) { break; }
  }
  return is_done;
}

void get_all_predecessors(Circuit &circ, Vertex &vertex, VertexSet &result) {
  for (Vertex predecessor : circ.get_predecessors(vertex)) {
    if (!is_initial_q_type(circ.get_OpType_from_Vertex(predecessor))) {
      result.insert(predecessor);
      get_all_predecessors(circ, predecessor, result);
    }
  }
}

Subcircuit get_max_partition(Circuit &circ, qubit_vector_t &qubits) {
  VertexVec inputs;
  VertexSet invalid_vertices;
  VertexSet max_partition;
  EdgeVec in_edges;
  EdgeVec out_edges;

  for (Qubit qubit : qubits) {
    inputs.push_back(circ.get_in(qubit));
  }
  // add valid input edges to the subcircuit's input edges and add invalid
  // inputs to the invalid_vertices set
  for (Vertex input : circ.all_inputs()) {
    if (std::count(inputs.begin(), inputs.end(), input) != 0) {
      in_edges.push_back(circ.get_nth_out_edge(input, 0));
    } else {
      invalid_vertices.insert(input);
    }
  }

  VertexVec vertices = circ.vertices_in_order();
  for (Vertex &v : vertices) {
    EdgeVec current_out_edges;
    bool isValid = true;
    VertexVec preds = circ.get_predecessors(v);
    // ignore partitions defined by boundary vertices
    if (is_boundary_q_type(circ.get_OpType_from_Vertex(v))) { continue; }
    // invalidate partition if its predecessors lead to an invalid input
    for (const Vertex pred : preds) {
      if (invalid_vertices.find(pred) != invalid_vertices.end()) {
          isValid = false;
          invalid_vertices.insert(v);
      }
    }
    // assign new max partition if it's valid
    if (isValid) {
      get_all_predecessors(circ, v, max_partition);
      max_partition.insert(v);
      // find output edges to current partition
      for (Vertex vert : max_partition) {
        const EdgeVec edges = circ.get_all_out_edges(vert);
        for (const Edge edge : edges) {
          Vertex next_vertex = circ.target(edge);
          if (max_partition.find(next_vertex) == max_partition.end()) {
            current_out_edges.push_back(edge);
          }
        }
      }
      out_edges = current_out_edges;
    }
}
  Subcircuit sub = {in_edges, out_edges, max_partition};
  return sub;
}

Partition synthesise(Partition &partition) { return partition; }

} // namespace tket