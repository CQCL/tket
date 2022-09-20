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

#include <catch2/catch_test_macros.hpp>

#include "NumericalOptimiser/CeresSolver.hpp"

namespace tket {
namespace test_CeresSolver {

// Some methods to check approximate equality
// Set d for desired number of decimal places to check
bool approx_equality(std::complex<double> a, std::complex<double> b) {
  int d = 5;
  double r = pow(10, d);
  if (round(a.real() * r)/r != round(b.real() * r)/r) {
    return false;
  }
  if (round(a.imag() * r)/r != round(b.imag() * r)/r) {
    return false;
  }
  return true; 
}

bool approx_matrix_equality(Eigen::MatrixXcd a, Eigen::MatrixXcd b) {
  if(a.cols() != b.cols()) { return false; }
  if(a.rows() != b.rows()) { return false; }
  
  for(int i=0; i<a.rows(); i++) {
    for(int j=0; j<a.cols(); j++) {
      if (!approx_equality(a(i, j), b(i, j))) {
        return false;
      }
    }
  }
  return true; 
}

SCENARIO("Testing optimise_u3") {
  GIVEN("Finding parameters for Hadamard gate") {
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);
    Eigen::MatrixXcd H(2, 2);
    H << 1/sqrt(2),  1/sqrt(2),
          1/sqrt(2), -1/sqrt(2);

    // TODO: stdout is being polluted by ceres info
    // uncomment this test when this has been fixed
    // std::vector<double> p = optimise_u3(0, ID, H);

    // REQUIRE(evaluate_distance(evaluate_u3(p[0], p[1], p[2], 0, 1), H) < 0.1);
  }
}

SCENARIO("Testing solve") {
  GIVEN("A few random parameters to improve") {
    std::uniform_real_distribution<double> unif(0, 2*M_PI);
    std::default_random_engine re;
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);
    double p[num_param], o[num_param];
    Eigen::MatrixXcd H(2, 2);
    H << 1/sqrt(2),  1/sqrt(2),
          1/sqrt(2), -1/sqrt(2);

    int num_tests = 3;
    for (int i=0; i<num_tests; i++) {
      for (int j=0; j<num_param; j++) {
        p[j] = unif(re);
      }
      std::copy(p, p+num_param, o); 
      // same as above comment
      // solve(0, ID, H, p);
      
      // REQUIRE(
      //   evaluate_distance(evaluate_u3(p[0], p[1], p[2], 0, 1), H) <
      //   evaluate_distance(evaluate_u3(o[0], o[1], o[2], 0, 1), H));
    }
  }
}

SCENARIO("Testing place") {
  GIVEN("an identity matrix in a single qubit circuit") {
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);

    Eigen::MatrixXcd placed = place(ID, 0, 1);

    REQUIRE(placed == ID);
  }
  GIVEN("a hadamard on the second qubit of a 2 qubit circuit") {
    Eigen::MatrixXcd H(2, 2), placed_H(4, 4);
    
    H << 1/sqrt(2),  1/sqrt(2),
         1/sqrt(2), -1/sqrt(2);
    
    placed_H << 1/sqrt(2), 0,  1/sqrt(2), 0,
                0, 1/sqrt(2),  0, 1/sqrt(2),
                1/sqrt(2), 0, -1/sqrt(2), 0,
                0, 1/sqrt(2),  0,-1/sqrt(2);

    Eigen::MatrixXcd placed = place(H, 1, 2);

    REQUIRE(placed == placed_H);
  }
}

SCENARIO("Testing evaluate_u3") {
  GIVEN("Parameters for a single qubit identity gate") {
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);
    
    Eigen::MatrixXcd U3 = evaluate_u3(0, 0, 0, 0, 1);

    REQUIRE(approx_matrix_equality(U3, ID));
  }
  GIVEN("Parameters for a hadamard gate on qubit 1 of 3 a qubit circuit") {
    Eigen::MatrixXcd H(4, 4);
    H << 1/sqrt(2), 0,  1/sqrt(2), 0,
         0, 1/sqrt(2),  0, 1/sqrt(2),
         1/sqrt(2), 0, -1/sqrt(2), 0,
         0, 1/sqrt(2),  0,-1/sqrt(2);

    Eigen::MatrixXcd U3 = evaluate_u3(M_PI/2, 0, M_PI, 1, 2);

    REQUIRE(approx_matrix_equality(U3, H));
  }
}

SCENARIO("Testing evaluate_distance") {
  std::complex<double> i(0, 1);
  Eigen::MatrixXcd X(2, 2), Y(2, 2), ID(2, 2);

    X << 0, 1,
         1, 0;

    Y << 0, -i,
         i,  0;

    ID << 1, 0,
          0, 1;

  GIVEN("two X gates") {
    double result = evaluate_distance(X, X);
    REQUIRE(result == 0.0);
  }
  GIVEN("an X and a Y gate") {
    double result = evaluate_distance(X, Y);
    REQUIRE(result == 1.0);
  }
  GIVEN("an X and an ID gate") {
    double result = evaluate_distance(X, ID);
    REQUIRE(result == 1.0);
  }
}

SCENARIO("Testing CostFunction constructor") {
  GIVEN("3 qubit identity circuit and target, gate at pos 1") {
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(8, 8);

    CircuitCostFunction* cost_function = new CircuitCostFunction(1, ID, ID);

    REQUIRE(cost_function->size == 3);
    REQUIRE(cost_function->pos == 1);
    REQUIRE(cost_function->U == ID);
    REQUIRE(cost_function->T == ID);
  }
}

SCENARIO("Testing evaluate_costs") {
  GIVEN("all parameters set to 0, identity circuit, identity target") {
    double p[3] = {0.0, 0.0, 0.0};
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);

    CircuitCostFunction* cost_function = new CircuitCostFunction(0, ID, ID);
    std::vector<double> costs = cost_function->evaluate_costs(p);

    REQUIRE(costs[0] == 0.0);
    REQUIRE(costs[1] == 0.0);
    REQUIRE(costs[2] == 0.0);
    REQUIRE(costs[3] == 0.0);
  }
  GIVEN("parameters pi/2, 0, pi (hadamard), identity circuit, X target") {
    double p[3] = {M_PI/2, 0.0, M_PI};
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);
    Eigen::MatrixXcd X(2, 2);

    X << 0, 1,
         1, 0;

    CircuitCostFunction* cost_function = new CircuitCostFunction(0, ID, X);
    std::vector<double> costs = cost_function->evaluate_costs(p);

    REQUIRE(approx_equality(costs[0], -sqrt(2)));
    REQUIRE(approx_equality(costs[1], 0.0));
    REQUIRE(approx_equality(costs[2], 0.0));
    REQUIRE(approx_equality(costs[3], 1/(2+sqrt(2))));
  }
}

SCENARIO("Testing evaluate_matrices") {
  GIVEN("all parameters set to 0") {
    std::complex<double> i(0, 1);
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);
    Eigen::MatrixXcd J1(2, 2), J2(2, 2), J3(2, 2), U(2, 2);

    J1 << 0, -0.5,
          0.5,  0;

    J2 << 0, 0,
          0, i;

    J3 << 0, 0,
          0, i;

    U = ID;

    CircuitCostFunction* cost_function = new CircuitCostFunction(0, ID, ID);
    std::vector<Eigen::MatrixXcd> matrices = cost_function->evaluate_matrices(0, 0, 0);
    
    REQUIRE(matrices[0] == J1);
    REQUIRE(matrices[1] == J2);
    REQUIRE(matrices[2] == J3);
    REQUIRE(matrices[3] == U);
  }
  GIVEN("parameters pi/2, 0, pi (hadamard)") {
    std::complex<double> i(0, 1);
    Eigen::MatrixXcd ID = Eigen::MatrixXcd::Identity(2, 2);
    Eigen::MatrixXcd J1(2, 2), J2(2, 2), J3(2, 2), U(2, 2);

    J1 << -sqrt(2)/4.0, sqrt(2)/4.0,
           sqrt(2)/4.0, sqrt(2)/4.0;

    J2 << 0, 0,
          i*sqrt(2)/2.0, -i*sqrt(2)/2.0;

    J3 << 0,  i*sqrt(2)/2.0,
          0, -i*sqrt(2)/2.0;

    U << sqrt(2)/2.0,  sqrt(2)/2.0,
         sqrt(2)/2.0, -sqrt(2)/2.0;

    CircuitCostFunction* cost_function = new CircuitCostFunction(0, ID, ID);
    std::vector<Eigen::MatrixXcd> matrices = cost_function->evaluate_matrices(M_PI/2, 0, M_PI);
    
    REQUIRE(approx_matrix_equality(matrices[0], J1));
    REQUIRE(approx_matrix_equality(matrices[1], J2));
    REQUIRE(approx_matrix_equality(matrices[2], J3));
    REQUIRE(approx_matrix_equality(matrices[3], U));
  }
}
}  // namespace test_CeresSolver
}  // namespace tket