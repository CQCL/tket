#pragma once

#include <vector>
#include "Circuit/Circuit.hpp"
#include "Architecture/Architecture.hpp"

#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/sized_cost_function.h>
#include <ceres/solver.h>

#include <iostream>
#include <fstream>

namespace tket {

std::ofstream f;

typedef std::pair<Circuit, qubit_vector_t> Partition;
typedef std::vector<Partition> PartitionVec;

struct node {
  Circuit circuit;
  double cost_estimate;
  double distance;
  int cnot_count;
};

Circuit optimise(Circuit &circ, Architecture &arch, unsigned int k);

PartitionVec partition(Circuit &circ, Architecture &arch, unsigned int k);

std::vector<node_set_t> get_connected_subarch(Architecture &arch, unsigned int k);

bool expand(
  node_set_t &current, node_set_t &to_expand, node_set_t &to_ignore, 
  Architecture &arch, unsigned int k, std::vector<node_set_t> &result);

void get_all_predecessors(Circuit &circ, Vertex &vertex, VertexSet &result);

Subcircuit get_max_partition(Circuit &circ, qubit_vector_t &qubits);

Partition synthesise(Partition &partition);

std::vector<double> optimise_circuit(
  int indexA, int indexB, Eigen::MatrixXcd &U, Eigen::MatrixXcd &T);

std::vector<double> optimise_u3_gates(
  int indexA, int indexB, Eigen::MatrixXcd &U, Eigen::MatrixXcd &T, double parameters[6]);

inline Eigen::MatrixXcd place(Eigen::MatrixXcd u, int pos, int size) {
  int a = std::pow(2, pos);
  int b = std::pow(2, (size-(pos+1)));
  Eigen::MatrixXcd above = Eigen::MatrixXcd::Identity(a, a);
  Eigen::MatrixXcd below = Eigen::MatrixXcd::Identity(b, b);
  return kroneckerProduct(below, kroneckerProduct(u, above).eval()).eval();
}

inline Eigen::MatrixXcd evaluate_u3(double x, double y, double z, int pos, int size) {
  std::complex<double> i(0, 1);
  Eigen::MatrixXcd u3(2, 2);
  u3(0, 0) = std::cos(x/2);
  u3(0, 1) = -std::sin(x/2) * std::exp(i*z);
  u3(1, 0) = std::sin(x/2) * std::exp(i*y);
  u3(1, 1) = std::cos(x/2) * std::exp(i*(y+z));
  u3 = place(u3, pos, size);
  return u3;
}

inline double evaluate_distance(Eigen::MatrixXcd U, Eigen::MatrixXcd T) {
  return 1 - std::abs((U.adjoint() * T).trace()) / T.cols();
}

class CircuitCostFunction : public ceres::SizedCostFunction<1, 6> {
public:

  int size;
  int indexA;
  int indexB;
  Eigen::MatrixXcd U;
  Eigen::MatrixXcd T;

  CircuitCostFunction(int indexA, int indexB, 
    Eigen::MatrixXcd U, Eigen::MatrixXcd T);

  virtual bool Evaluate(double const* const* parameters,
      double* residuals, double** jacobians) const;

  std::vector<double> evaluate_jacs(const double* p) const;

  std::vector<Eigen::MatrixXcd> jac_matrices(double x, double y, double z) const;

};

} // namespace tket