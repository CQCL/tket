#include "PostRoutingOptimisation.hpp"
#include <algorithm>


#include <string.h>

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

std::vector<double> optimise_u3_gates(
  int indexA, int indexB, Eigen::MatrixXcd &U, Eigen::MatrixXcd &T) {
  // Set up problem
  f.open("output.txt");
  ceres::Problem problem;
  double parameters[6] = {0.1, 0.0, 0.0, 0.0, 0.0, 0.0};
  ceres::CostFunction* cost_function = new CircuitCostFunction(indexA, indexB, U, T);
  problem.AddResidualBlock(cost_function, nullptr, parameters);
  
  // Setting range of 0 to 2pi for possible angles
  for (int i=0; i<6; i++) {
    problem.SetParameterLowerBound(parameters, i, 0);
    problem.SetParameterUpperBound(parameters, i, 2*M_PI);
  }
  
  // Solve for optimal U3 gates
  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;
  Solve(options, &problem, &summary);

  f << summary.BriefReport() << "\n";
  f.close();

  std::vector<double> result(std::begin(parameters), std::end(parameters));
  return result;
}

CircuitCostFunction::CircuitCostFunction(int indexA, int indexB, 
  Eigen::MatrixXcd U, Eigen::MatrixXcd T) {
  this->size = log2(T.cols());
  this->indexA = indexA;
  this->indexB = indexB;
  this->U = U;
  this->T = T;
}

bool CircuitCostFunction::Evaluate(double const* const* parameters,
  double* residuals, double** jacobians) const {
  f << "x1: " << parameters[0][0] << "\n";
  f << "y1: " << parameters[0][1] << "\n";
  f << "z1: " << parameters[0][2] << "\n";
  f << "x2: " << parameters[0][3] << "\n";
  f << "y2: " << parameters[0][4] << "\n";
  f << "z2: " << parameters[0][5] << "\n";
  
  std::vector<double> jacs = evaluate_distance(parameters[0]);

  residuals[0] = jacs[6];

  f << "distance: " << residuals[0] << "\n";

  if (jacobians != nullptr && jacobians[0] != nullptr) {
    for (int i=0; i<6; i++) {
      jacobians[0][i] = jacs[i];
      f << "Jacobian " << i << ": " << jacobians[0][i] << "\n";
    }
  }
  
  return true;
}

std::vector<double> CircuitCostFunction::evaluate_distance(const double* p) const {
  std::vector<Eigen::MatrixXcd> A = u3_matrices(p[0], p[1], p[2]);
  std::vector<Eigen::MatrixXcd> B = u3_matrices(p[3], p[4], p[5]);
  Eigen::MatrixXcd C = A[0]*B[0]*U;

  std::complex<double> S = T.cwiseProduct(C.conjugate()).sum();
  double dsq = 1 - std::abs(S)/T.cols();
  
  f << "S: " << S << "\n";
  f << "dsq: " << dsq << "\n";

  std::vector<Eigen::MatrixXcd> ju;
  std::vector<std::complex<double>> jus;
  std::vector<double> jacs;

  for(int i=0; i<3; i++) {
    place(A[i+1], indexA);
    ju.push_back(T.cwiseProduct(A[i+1].conjugate()));
  }

  for(int i=0; i<3; i++) {
    place(B[i+1], indexB);
    ju.push_back(T.cwiseProduct(B[i+1].conjugate()));
  }

  for(int i=0; i<6; i++) {
    jus.push_back(ju[i].sum());
  }

  f << "JUS 0:\n" << jus[0] << "\n";
  f << "JUS 1:\n" << jus[1] << "\n";
  f << "JUS 2:\n" << jus[1] << "\n";

  for(int i=0; i<6; i++) {
    jacs.push_back(-(S.real()*jus[i].real() + S.imag()*jus[i].imag())*T.cols() / std::abs(S));
    f << "jacs[i]: " << jacs[i] << "\n";
  }

  jacs.push_back(dsq);
  
  return jacs;
}

std::vector<Eigen::MatrixXcd> CircuitCostFunction::u3_matrices(double x, double y, double z) const {
  std::complex<double> i(0, 1);
  double ct = std::cos(x/2);
  double st = std::sin(x/2);
  double cp = std::cos(y);
  double sp = std::sin(y);
  double cl = std::cos(z);
  double sl = std::sin(z);

  Eigen::MatrixXcd U(2,2), J1(2,2), J2(2,2), J3(2,2);
  
  U(0, 0) = ct;
  U(0, 1) = -st * (cl + i * sl);
  U(1, 0) = st * (cp + i * sp);
  U(1, 1) = ct * (cl * cp - sl * sp + i * cl * sp + i * sl * cp);

  J1(0, 0) = -0.5*st;
  J1(0, 1) = -0.5*ct * (cl + i * sl);
  J1(1, 0) = 0.5*ct * (cp + i * sp);
  J1(1, 1) = -0.5*st * (cl * cp - sl * sp + i * cl * sp + i * sl * cp);

  J2(0, 0) = 0;
  J2(0, 1) = 0;
  J2(1, 0) = st *(-sp + i * cp);
  J2(1, 1) = ct *(cl * -sp - sl * cp + i * cl * cp + i * sl * -sp);

  J3(0, 0) = 0;
  J3(0, 1) = -st *(-sl + i * cl);
  J3(1, 0) = 0;
  J3(1, 1) = ct *(-sl * cp - cl * sp + i * -sl * sp + i * cl * cp);

  std::vector<Eigen::MatrixXcd> v;
  v.push_back(U);
  v.push_back(J1);
  v.push_back(J2);
  v.push_back(J3);

  return v;
}

void CircuitCostFunction::place(Eigen::MatrixXcd u, int pos) const {
  Eigen::MatrixXcd below = Eigen::MatrixXcd::Identity(std::pow(2, (size-(pos+1))), std::pow(2, (size-(pos+1))));
  Eigen::MatrixXcd above = Eigen::MatrixXcd::Identity(std::pow(2, pos), std::pow(2, pos));

  kroneckerProduct(below, kroneckerProduct(u, above).eval()).eval();
}

} // namespace tket