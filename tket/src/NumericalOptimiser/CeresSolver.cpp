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

#include "CeresSolver.hpp"

namespace tket {

std::vector<double> optimise_u3(
  int pos, Eigen::MatrixXcd U, Eigen::MatrixXcd T) {
  double lower_bound = 0;
  double upper_bound = 2*M_PI;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::uniform_real_distribution<double> adj(0, local_adjustment);
  std::default_random_engine re;
  double best_cost = 1.0;
  std::vector<double> best_params;

  auto compare = [](std::vector<double> lhs, std::vector<double> rhs) { 
    return lhs[num_param] > rhs[num_param];
  };
  std::priority_queue<std::vector<double>, std::vector<std::vector<double>>, 
    decltype(compare)> param_queue(compare);
    
  for(int i=0; i<num_global_start; i++) {
    std::vector<double> p;
    for(int j=0; j<num_param; j++) {
      p.push_back(unif(re));
    }
    Eigen::MatrixXcd U3 = evaluate_u3(p[0], p[1], p[2], pos, log2(T.cols()));
    p.push_back(evaluate_distance(U3*U, T));
    param_queue.push(p);
  }

  std::vector<double> start = param_queue.top();
  for(int j=0; j<num_local_start; j++) {
    double parameters[num_param];
    for(int k=0; k<num_param; k++) {
      parameters[k] = start[k] + (adj(re) - local_adjustment/2.0);
    }
    std::vector<double> result = solve(pos, U, T, parameters);
    if (result[num_param] < best_cost) { 
      best_cost = result[num_param];
      best_params = result;
    }
  }

  return best_params;
}

std::vector<double> solve(int pos, Eigen::MatrixXcd U, 
  Eigen::MatrixXcd T, double parameters[num_param]) {
  ceres::Problem problem;
  CircuitCostFunction* cost_function = new CircuitCostFunction(pos, U, T);
  std::vector<ceres::ResidualBlockId> residual_block_ids;
  ceres::ResidualBlockId block_id = problem.AddResidualBlock(cost_function, nullptr, parameters);
  residual_block_ids.push_back(block_id);
  
  for (int i=0; i<num_param; i++) {
    problem.SetParameterLowerBound(parameters, i, 0);
    problem.SetParameterUpperBound(parameters, i, 2*M_PI);
  }
  
  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = false;
  options.linear_solver_type = ceres::DENSE_QR;
  options.function_tolerance = 5e-16;
  options.gradient_tolerance = 1e-15;
  // options.num_threads = num_threads;
  // options.max_num_iterations = max_iters;

  ceres::Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::vector<double> result;
  for(int i=0; i<num_param; i++)
    result.push_back(parameters[i]);
  
  Eigen::MatrixXcd u3 = evaluate_u3(result[0], result[1], result[2], pos, log2(T.cols()));
  result.push_back(evaluate_distance(u3*U, T));

  return result;
}

CircuitCostFunction::CircuitCostFunction(int pos, Eigen::MatrixXcd U, Eigen::MatrixXcd T) {
  this->size = log2(T.cols());
  this->pos = pos;
  this->U = U;
  this->T = T;
}

bool CircuitCostFunction::Evaluate(double const* const* parameters,
  double* residuals, double** jacobians) const {
  std::vector<double> costs = evaluate_costs(parameters[0]);
  residuals[0] = costs[num_param];
  
  if (jacobians != nullptr && jacobians[0] != nullptr) {
    for (int i=0; i<num_param; i++) {
      jacobians[0][i] = costs[i];
    }
  }
  
  return true;
}

std::vector<double> CircuitCostFunction::evaluate_costs(const double* p) const {
  std::vector<Eigen::MatrixXcd> U3 = evaluate_matrices(p[0], p[1], p[2]);
  std::vector<Eigen::MatrixXcd> ju;
  std::vector<std::complex<double>> jus;
  std::vector<double> costs;
  
  for(int i=0; i<num_param+1; i++) {
    U3[i] = place(U3[i], pos, size);
  }

  Eigen::MatrixXcd C = U3[num_param]*U;
  std::complex<double> S = T.cwiseProduct(C.conjugate()).sum();
  double dsq = 1 - std::abs(S)/T.cols();

  for(int i=0; i<num_param; i++) {
    ju.push_back(T.cwiseProduct(U3[i].conjugate()));
    jus.push_back(ju[i].sum());
    costs.push_back(
      -(S.real()*jus[i].real() + S.imag()*jus[i].imag())*T.cols() / std::abs(S));
  }

  costs.push_back(dsq);
  
  return costs;
}

std::vector<Eigen::MatrixXcd> CircuitCostFunction::evaluate_matrices(
  double x, double y, double z) const {
  std::complex<double> i(0, 1);
  double cx = std::cos(x/2);
  double sx = std::sin(x/2);
  double cy = std::cos(y);
  double sy = std::sin(y);
  double cz = std::cos(z);
  double sz = std::sin(z);

  Eigen::MatrixXcd U(2,2), J1(2,2), J2(2,2), J3(2,2);
  
  U(0, 0) = cx;
  U(0, 1) = -sx * (cz + i * sz);
  U(1, 0) = sx * (cy + i * sy);
  U(1, 1) = cx * (cz * cy - sz * sy + i * cz * sy + i * sz * cy);

  J1(0, 0) = -0.5*sx;
  J1(0, 1) = -0.5*cx * (cz + i * sz);
  J1(1, 0) = 0.5*cx * (cy + i * sy);
  J1(1, 1) = -0.5*sx * (cz * cy - sz * sy + i * cz * sy + i * sz * cy);

  J2(0, 0) = 0;
  J2(0, 1) = 0;
  J2(1, 0) = sx *(-sy + i * cy);
  J2(1, 1) = cx *(cz * -sy - sz * cy + i * cz * cy + i * sz * -sy);

  J3(0, 0) = 0;
  J3(0, 1) = -sx *(-sz + i * cz);
  J3(1, 0) = 0;
  J3(1, 1) = cx *(-sz * cy - cz * sy + i * -sz * sy + i * cz * cy);

  std::vector<Eigen::MatrixXcd> v;
  v.push_back(J1);
  v.push_back(J2);
  v.push_back(J3);
  v.push_back(U);

  return v;
}
} // namespace tket