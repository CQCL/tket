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
#include <algorithm>
#include <queue>
#include <random>
#include <string.h>

#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/sized_cost_function.h>
#include <ceres/solver.h>

#include "Circuit/Circuit.hpp"
#include "Architecture/Architecture.hpp"

namespace tket {
// number of parametes in gate to optimise
static const int num_param = 3;
// number of initial starting points to evaluate
static const int num_global_start = 100;
// number of start locations near the best initial start point
static const int num_local_start = 10;
// range by which local start points may vary
static const double local_adjustment = 0.1;

/**
 * @brief Finds optimum parameters for a U3 gate given an existing circuit,
 * a target circuit and the index of the qubit the U3 is being added to. 
 * This is achieved by evaluating random start points and running the ceres 
 * solver on the best ones.
 * 
 * @param pos the index of the qubit we want to add a U3 gate to
 * @param U the unitary representing the circuit we want to add a U3 gate to
 * @param T the target unitary to implement
 * @return a vector containing the optimum parameters (x, y, z)
 */
std::vector<double> optimise_u3(
  int pos, Eigen::MatrixXcd U, Eigen::MatrixXcd T);

/**
 * @brief Creates and runs the ceres solver to optimise the parameters of a
 * U3 gate with a specified start point,
 * 
 * @param pos the index of the qubit we want to add a U3 gate to
 * @param U the unitary representing the circuit we want to add a U3 gate to
 * @param T the target unitary to implement
 * @param parameters the start point of each parameter to solve
 * @return a vector containing the optimum parameters (x, y, z) and the new
 * distance of the circuit to the target
 */
std::vector<double> solve(int pos, Eigen::MatrixXcd U, 
  Eigen::MatrixXcd T, double parameters[num_param]);

/**
 * @brief Places the matrix representation of a gate on a given qubit of
 * the matrix representation of an empty circuit.
 * 
 * @param u the matrix representation of a gate
 * @param pos the index of the qubit to place the gate on
 * @param size the number of qubits in the circuit
 * @return the matrix representation of the gate placed on an empty circuit
 */
inline Eigen::MatrixXcd place(Eigen::MatrixXcd u, int pos, int size) {
  int a = std::pow(2, pos);
  int b = std::pow(2, (size-(pos+1)));
  Eigen::MatrixXcd above = Eigen::MatrixXcd::Identity(a, a);
  Eigen::MatrixXcd below = Eigen::MatrixXcd::Identity(b, b);
  return kroneckerProduct(below, kroneckerProduct(u, above).eval()).eval();
}

/**
 * @brief Computes the matrix representation of a U3 gate given its parameters
 * and places it on the specified qubit of an empty circuit. Note: this is
 * also done in evaluate matrices but this inline method only computes what
 * is necessary for a quick evaluation of the initial starting parameters
 * before actually running the solver.
 * 
 * @param x theta
 * @param y phi
 * @param z lambda
 * @param pos the index of the qubit to place the gate on
 * @param size the number of qubits in the cirucit it will be placed on.
 * @return the matrix representation of the U3 gate placed on an empty circuit
 */
inline Eigen::MatrixXcd evaluate_u3(
  double x, double y, double z, int pos, int size) {
  std::complex<double> i(0, 1);
  Eigen::MatrixXcd u3(2, 2);
  u3(0, 0) = std::cos(x/2);
  u3(0, 1) = -std::sin(x/2) * std::exp(i*z);
  u3(1, 0) = std::sin(x/2) * std::exp(i*y);
  u3(1, 1) = std::cos(x/2) * std::exp(i*(y+z));
  u3 = place(u3, pos, size);
  return u3;
}
/**
 * @brief Computes the distance between two unitaries as a variation of the
 * Hilbert-Schmidt norm. Has useful properties of being 0 when compilation
 * is exact and being fast to compute.
 * 
 * @param U First unitary to compare
 * @param T Second unitary to compare
 * @return distance between the two unitaries
 */
inline double evaluate_distance(Eigen::MatrixXcd U, const Eigen::MatrixXcd T) {
  return 1 - std::abs((U.adjoint() * T).trace()) / T.cols();
}

/** Cost function defining the optimisation problem of finding parameters
 * for a U3 gate.
 */
class CircuitCostFunction : public ceres::SizedCostFunction<1, num_param> {
public:
  /** number of qubits in the circuit */
  int size;
  /** index of the qubit the u3 gate is being added to */
  int pos;
  /** the unitary representing the circuit we want to add a U3 gate to */
  Eigen::MatrixXcd U;
  /** the target unitary to implement */
  Eigen::MatrixXcd T;

  CircuitCostFunction(int pos, Eigen::MatrixXcd U, Eigen::MatrixXcd T);

  /**
   * @brief Called by ceres to evaluate the cost and jacobians at a specific
   * point on the cost function. Overriden to store the correct cost in 
   * residuals[0] and the jacobians in jacobians[0].
   * 
   * @param parameters the point at which to evaluate the cost function
   * @param residuals an array of residuals to set. Our problem has a 
   * single residual so we only need to set residuals[0]
   * @param jacobians an array of jacobians to set. Our problem has 3
   * parameters so set jacobians[0][0], jacobians[0][1], and jacobians[0][2]
   * @return true: successful evaluation
   * @return false: unsucessful evaluation
   */
  virtual bool Evaluate(double const* const* parameters,
      double* residuals, double** jacobians) const;

  /**
   * @brief Evaluates the cost and jacobians at a specified point
   * 
   * @param p array of parameters defining the point to evaluate
   * @return array containing [jacobian x, jacobian y, jacobian z, cost]
   */
  std::vector<double> evaluate_costs(const double* p) const;

  /**
   * @brief Computes the matrix of a u3 gate given its parameters and the
   * matrices of the partial derivatives with respect to each parameter.
   * 
   * @param x theta
   * @param y phi
   * @param z lambda
   * @return array containing [d cost/dx, d cost/dy, d cost/dz, cost]
   */
  std::vector<Eigen::MatrixXcd> evaluate_matrices(double x, double y, double z) const;
};
} // namespace tket