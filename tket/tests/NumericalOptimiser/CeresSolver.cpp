namespace tket {

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