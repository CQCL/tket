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

Circuit optimise(Circuit &circ, Architecture &arch, unsigned int k);

PartitionVec partition(Circuit &circ, Architecture &arch, unsigned int k);

std::vector<node_set_t> get_connected_subarch(Architecture &arch, unsigned int k);

bool expand(
  node_set_t &current, node_set_t &to_expand, node_set_t &to_ignore, 
  Architecture &arch, unsigned int k, std::vector<node_set_t> &result);

void get_all_predecessors(Circuit &circ, Vertex &vertex, VertexSet &result);

Subcircuit get_max_partition(Circuit &circ, qubit_vector_t &qubits);

Partition synthesise(Partition &partition);

std::vector<double> optimise_u3_gates(
  int indexA, int indexB, Eigen::MatrixXcd &U, Eigen::MatrixXcd &T);

class CircuitCostFunction : public ceres::SizedCostFunction<1, 6> {
public:
  // number of qubits in the circuit
  int size;
  // index of the first u3 gate's qubit
  int indexA;
  // index of the second u3 gate's qubit
  int indexB;
  // unitary of the circuit constructed so far
  Eigen::MatrixXcd U;
  // target unitary for the circuit
  Eigen::MatrixXcd T;

  CircuitCostFunction(int indexA, int indexB, 
    Eigen::MatrixXcd U, Eigen::MatrixXcd T);

  virtual bool Evaluate(double const* const* parameters,
      double* residuals, double** jacobians) const;

  std::vector<double> evaluate_distance(const double* p) const;

  std::vector<Eigen::MatrixXcd> u3_matrices(double x, double y, double z) const;

  void place(Eigen::MatrixXcd u, int pos) const;
};

} // namespace tket