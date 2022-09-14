#include "PostRoutingOptimisation.hpp"
#include <algorithm>
#include <queue>
#include <random>

#include <chrono>

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

Partition synthesise(Partition &partition) {
  // define queue for A* search
  auto compare_node = [](std::vector<node> lhs, std::vector<node> rhs) { 
    return lhs.cost_estimate > rhs.cost_estimate;
  };
  std::priority_queue<std::vector<node>, std::vector<std::vector<node>>, 
    decltype(compare_node)> param_queue(compare_node);
  // add initial u3 gates
  Cicuit initial(partition.second.size());
  // while head distance is less than threshold
  // for each connected pair
  // generate successor
  // evaluate successor
  // calculate f(n) = cnot count +  9.3623 * distance
  // insert node in queue

  return partition;
}

std::vector<double> optimise_circuit(
  int indexA, int indexB, Eigen::MatrixXcd &U, Eigen::MatrixXcd &T) {
  double lower_bound = 0;
  double upper_bound = 2*M_PI;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound), adj(0, 0.1);
  std::default_random_engine re;
  double best_cost = 10.0;
  std::vector<double> best_params;

  f.open("output.txt");
  // temporary algorithm: evaluate 1000 starting points, choose 10 best,
  // choose 10 points near each best and optimise each
  // TODO: add fun techniques to improve optimisation :)
  auto compare = [](std::vector<double> lhs, std::vector<double> rhs) { 
    return lhs[6] > rhs[6];
  };
  std::priority_queue<std::vector<double>, std::vector<std::vector<double>>, 
    decltype(compare)> param_queue(compare);
    
  auto start = std::chrono::high_resolution_clock::now();
  for(int i=0; i<1000; i++) {
    std::vector<double> p;
    for(int j=0; j<6; j++) {
      p.push_back(unif(re));
    }
    Eigen::MatrixXcd A = evaluate_u3(p[0], p[1], p[2], indexA, T.cols());
    Eigen::MatrixXcd B = evaluate_u3(p[3], p[4], p[5], indexB, T.cols());
    p.push_back(evaluate_distance(A*B*U, T));
    param_queue.push(p);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  f << "duration: " << duration.count()/1000000.0 << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for(int i=0; i<10; i++) {
    std::vector<double> initial = param_queue.top();
    param_queue.pop();
    for(int j=0; j<10; j++) {
      double parameters[6];
      for(int k=0; k<6; k++) {
        parameters[k] = initial[k] + (adj(re) - 0.05);
      }
      std::vector<double> result = optimise_u3_gates(indexA, indexB, U, T, parameters);
      if (result[6] < best_cost) { 
        best_cost = result[6];
        best_params = result;
      }
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  f << "duration: " << duration.count()/1000000.0 << std::endl;
  f.close();
  return best_params;
}

std::vector<double> optimise_u3_gates(
  int indexA, int indexB, Eigen::MatrixXcd &U, Eigen::MatrixXcd &T, double parameters[6]) {
  // Set up problem
  ceres::Problem problem;
  ceres::CostFunction* cost_function = new CircuitCostFunction(indexA, indexB, U, T);
  std::vector<ceres::ResidualBlockId> residual_block_ids;
  ceres::ResidualBlockId block_id = problem.AddResidualBlock(cost_function, nullptr, parameters);
  residual_block_ids.push_back(block_id);
  
  // Setting range of 0 to 2pi for possible angles
  for (int i=0; i<6; i++) {
    problem.SetParameterLowerBound(parameters, i, 0);
    problem.SetParameterUpperBound(parameters, i, 2*M_PI);
  }
  
  // Solve for optimal U3 gates
  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = false;
  options.linear_solver_type = ceres::DENSE_QR;
  options.function_tolerance = 5e-16;
  options.gradient_tolerance = 1e-15;
  // options.num_threads = num_threads;
  // options.max_num_iterations = max_iters;

  ceres::Solver::Summary summary;
  Solve(options, &problem, &summary);

  //f << summary.BriefReport() << "\n";

  std::vector<double> result;
  for(int i=0; i<6; i++)
    result.push_back(parameters[i]);
  
  double total_cost;
  std::vector<double> residuals;
  ceres::Problem::EvaluateOptions evaluate;
  evaluate.residual_blocks = residual_block_ids;
  problem.Evaluate(evaluate, &total_cost, &residuals, nullptr, nullptr);
  result.push_back(total_cost);
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
  std::vector<double> jacs = evaluate_jacs(parameters[0]);
  residuals[0] = jacs[6];
  
  if (jacobians != nullptr && jacobians[0] != nullptr) {
    for (int i=0; i<6; i++) {
      jacobians[0][i] = jacs[i];
    }
  }
  
  return true;
}

std::vector<double> CircuitCostFunction::evaluate_jacs(const double* p) const {
  std::vector<Eigen::MatrixXcd> A = jac_matrices(p[0], p[1], p[2]);
  std::vector<Eigen::MatrixXcd> B = jac_matrices(p[3], p[4], p[5]);
  
  for(int i=0; i<4; i++) {
    A[i] = place(A[i], indexA, size);
    B[i] = place(B[i], indexB, size);
  }

  Eigen::MatrixXcd C = A[0]*B[0]*U;

  std::complex<double> S = T.cwiseProduct(C.conjugate()).sum();
  double dsq = 1 - std::abs(S)/T.cols();

  std::vector<Eigen::MatrixXcd> ju;
  std::vector<std::complex<double>> jus;
  std::vector<double> jacs;

  for(int i=0; i<3; i++) {
    ju.push_back(T.cwiseProduct(A[i+1].conjugate()));
  }

  for(int i=0; i<3; i++) {
    ju.push_back(T.cwiseProduct(B[i+1].conjugate()));
  }

  for(int i=0; i<6; i++) {
    jus.push_back(ju[i].sum());
  }

  for(int i=0; i<6; i++) {
    jacs.push_back(-(S.real()*jus[i].real() + S.imag()*jus[i].imag())*T.cols() / std::abs(S));
  }

  jacs.push_back(dsq);
  
  return jacs;
}

std::vector<Eigen::MatrixXcd> CircuitCostFunction::jac_matrices(double x, double y, double z) const {
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
} // namespace tket