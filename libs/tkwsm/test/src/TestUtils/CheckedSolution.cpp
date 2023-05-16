// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "CheckedSolution.hpp"

#include <catch2/catch_test_macros.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/EndToEndWrappers/MainSolver.hpp>

#include "../TestUtils/TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static void check_with_no_solution(
    CheckedSolution::ProblemInformation info,
    // const MainSolverParameters& solver_params,
    const SolutionData& soln_data, CheckedSolution::Statistics& stats,
    const OStreamWrapper& os) {
  if (soln_data.finished) {
    os << "; no soln.";
    switch (info.existence) {
      case CheckedSolution::ProblemInformation::SolutionsExistence::
          KNOWN_TO_BE_INSOLUBLE:
        // fallthrough.
      case CheckedSolution::ProblemInformation::SolutionsExistence::UNKNOWN:
        ++stats.success_count;
        return;
      case CheckedSolution::ProblemInformation::SolutionsExistence::
          KNOWN_TO_BE_SOLUBLE:
        ++stats.failure_count;
        return;
      default:
        CHECK(false);
    }
  } else {
    os << "; TIMED OUT";
    ++stats.timeout_count;
  }
}

static bool check_solution_scalar_product(
    CheckedSolution::ProblemInformation info,
    // const MainSolverParameters& solver_params,
    const SolutionData& solution_data,
    // CheckedSolution::Statistics& stats,
    // CheckedSolution& checked_solution,
    const SolutionWSM& best_solution, const OStreamWrapper& os) {
  // OK, maybe without realising the caller might just happen to pass in
  // an unweighted problem, but ignore that.
  const bool unweighted_problem =
      info.known_lower_bound && info.known_upper_bound &&
      info.known_lower_bound.value() == info.known_upper_bound.value();

  auto lower_bound = solution_data.trivial_weight_lower_bound;
  if (info.known_lower_bound) {
    lower_bound = std::max(lower_bound, info.known_lower_bound.value());
  }
  auto upper_bound = solution_data.trivial_weight_initial_upper_bound;
  if (info.known_upper_bound) {
    // Do we expect an OPTIMAL solution?
    if (solution_data.finished || unweighted_problem) {
      upper_bound = std::min(upper_bound, info.known_upper_bound.value());
    }
  }
  if (!(lower_bound <= best_solution.scalar_product &&
        best_solution.scalar_product <= upper_bound)) {
    os << ": CALC " << best_solution.scalar_product
       << " violates known soln bounds [" << lower_bound << ", " << upper_bound
       << "]";
    return false;
  }
  return true;
}

static void check_calculated_solution(
    CheckedSolution::ProblemInformation info,
    const MainSolverParameters& solver_params, const SolutionData& soln_data,
    CheckedSolution::Statistics& stats, CheckedSolution& checked_solution,
    const OStreamWrapper& os) {
  REQUIRE(soln_data.solutions.size() == 1);
  const auto& best_solution = soln_data.solutions[0];
  checked_solution.scalar_product = best_solution.scalar_product;

  checked_solution.assignments = best_solution.assignments;
  REQUIRE(is_sorted_and_unique(checked_solution.assignments));
  REQUIRE(!checked_solution.assignments.empty());
  if (!info.known_optimal_solution) {
    os << "; soln " << checked_solution.scalar_product;
  }
  const bool valid_scalar_product =
      check_solution_scalar_product(info, soln_data, best_solution, os);
  if (soln_data.finished || solver_params.terminate_with_first_full_solution) {
    if (valid_scalar_product) {
      ++stats.success_count;
    } else {
      ++stats.failure_count;
    }
  } else {
    os << "; TIMED OUT";
    ++stats.timeout_count;
  }
}

static void check_known_solution_information(
    CheckedSolution::ProblemInformation& info) {
  if (info.known_optimal_solution) {
    REQUIRE(
        info.existence != CheckedSolution::ProblemInformation::
                              SolutionsExistence::KNOWN_TO_BE_INSOLUBLE);
    info.existence = CheckedSolution::ProblemInformation::SolutionsExistence::
        KNOWN_TO_BE_SOLUBLE;
    if (info.known_lower_bound) {
      REQUIRE(
          info.known_lower_bound.value() <=
          info.known_optimal_solution.value());
    }
    if (info.known_upper_bound) {
      REQUIRE(
          info.known_upper_bound.value() >=
          info.known_optimal_solution.value());
    }
    info.known_lower_bound = info.known_optimal_solution.value();
    info.known_upper_bound = info.known_optimal_solution.value();
    return;
  }
  if (info.known_lower_bound && info.known_upper_bound) {
    REQUIRE(info.known_lower_bound.value() <= info.known_upper_bound.value());
    if (info.known_lower_bound.value() == info.known_upper_bound.value()) {
      info.known_optimal_solution = info.known_lower_bound.value();
    }
  }
}

static void check_for_impossible_weight_constraint(
    CheckedSolution::ProblemInformation& info,
    const MainSolverParameters& solver_params) {
  if (!solver_params.weight_upper_bound_constraint || !info.known_lower_bound) {
    return;
  }
  const auto imposed_upper_bound =
      solver_params.weight_upper_bound_constraint.value();
  const auto known_lower_bound = info.known_lower_bound.value();
  if (known_lower_bound <= imposed_upper_bound) {
    return;
  }
  // The problem has now become insoluble, with the extra weight constraint.
  info.known_optimal_solution = {};
  info.known_lower_bound = {};
  info.known_upper_bound = {};
  info.existence = CheckedSolution::ProblemInformation::SolutionsExistence::
      KNOWN_TO_BE_INSOLUBLE;
}

static void check_solver_object(
    const MainSolver& solver, const GraphEdgeWeights& pdata,
    const GraphEdgeWeights& tdata, CheckedSolution::ProblemInformation info,
    const MainSolverParameters& solver_params,
    CheckedSolution::Statistics& stats, CheckedSolution& checked_solution,
    const OStreamWrapper& os) {
  const auto& solution_data = solver.get_solution_data();
  checked_solution.iterations = solution_data.iterations;
  checked_solution.finished = solution_data.finished;
  CHECK(solution_data.solutions.size() <= 1);

  stats.total_init_time_ms += solution_data.initialisation_time_ms;
  stats.total_search_time_ms += solution_data.search_time_ms;
  stats.total_iterations += solution_data.iterations;

  if (TestSettings::get().print_solution_times) {
    os << "; time " << solution_data.initialisation_time_ms << "+"
       << solution_data.search_time_ms;
  }
  os << "; " << solution_data.iterations << " iters";
  if (info.known_optimal_solution) {
    os << "; known opt.val. " << info.known_optimal_solution.value();
  }

  for (const auto& solution : solution_data.solutions) {
    const auto errors = solution.get_errors(pdata, tdata);
    REQUIRE("" == errors);
    if (!errors.empty()) {
      ++stats.failure_count;
      return;
    }
  }

  if (solution_data.solutions.empty()) {
    check_with_no_solution(info, solution_data, stats, os);
  } else {
    check_calculated_solution(
        info, solver_params, solution_data, stats, checked_solution, os);
  }
}

static void solve_problem(
    const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
    CheckedSolution::ProblemInformation info,
    const MainSolverParameters& solver_params,
    CheckedSolution::Statistics& stats, CheckedSolution& checked_solution,
    const OStreamWrapper& os) {
  check_known_solution_information(info);
  check_for_impossible_weight_constraint(info, solver_params);
  const MainSolver solver(pdata, tdata, solver_params);
  check_solver_object(
      solver, pdata, tdata, info, solver_params, stats, checked_solution, os);
  const auto& extra_statistics = solver.get_solution_data().extra_statistics;
  if (!extra_statistics.impossible_target_vertices.empty()) {
    checked_solution.impossible_target_vertices =
        extra_statistics.impossible_target_vertices;
  }
  if (!TestSettings::get().print_verbose_solution_data) {
    return;
  }

  if (extra_statistics.number_of_pattern_vertices != 0 &&
      extra_statistics.number_of_target_vertices != 0) {
    TestSettings::get().os
        << "; tot.ass: "
        << extra_statistics.number_of_pattern_vertices *
               extra_statistics.number_of_target_vertices
        << "->(" << extra_statistics.initial_number_of_possible_assignments
        << "," << extra_statistics.total_number_of_assignments_tried << ","
        << extra_statistics.total_number_of_impossible_assignments << ")";
  }
  if (!extra_statistics.impossible_target_vertices.empty()) {
    TestSettings::get().os << "; imposs.tv: "
                           << str(extra_statistics.impossible_target_vertices);
  }
  if (extra_statistics.n_tv_initially_passed_to_weight_nogood_detector) {
    TestSettings::get().os
        << "; wngd.tv.: "
        << extra_statistics.n_tv_initially_passed_to_weight_nogood_detector
               .value()
        << "->";
    if (extra_statistics.n_tv_still_valid_in_weight_nogood_detector) {
      TestSettings::get().os
          << extra_statistics.n_tv_still_valid_in_weight_nogood_detector
                 .value();
    } else {
      TestSettings::get().os << "?";
    }
  }
}

CheckedSolution::Statistics::Statistics(const std::string& test_name)
    : success_count(0),
      failure_count(0),
      timeout_count(0),
      total_init_time_ms(0),
      total_search_time_ms(0),
      total_iterations(0) {
  TestSettings::get().os << "\n##### BEGIN test '" << test_name
                         << "' ########\n";
}

CheckedSolution::Statistics::Statistics(
    const std::string& test_name, std::size_t number_of_graphs)
    : Statistics(
          test_name + "; " + std::to_string(number_of_graphs) + " graphs") {}

void CheckedSolution::Statistics::finish(Expectation expectation) const {
  TestSettings::get().os << "\n##### END test: " << success_count
                         << " successes, " << failure_count << " failures, "
                         << timeout_count << " timeouts;\nTotal time "
                         << total_init_time_ms << "+" << total_search_time_ms
                         << " ms; total " << total_iterations
                         << " iterations\n";

  if (expectation == Expectation::ALL_SUCCESS) {
    CHECK(failure_count == 0);
    CHECK(timeout_count == 0);
  }
  if (expectation == Expectation::ALL_SUCCESS_OR_TIMEOUT) {
    CHECK(failure_count == 0);
  }
}

CheckedSolution::CheckedSolution(
    const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
    ProblemInformation info, const MainSolverParameters& solver_params,
    Statistics& stats) {
  REQUIRE(
      solver_params.for_multiple_full_solutions_the_max_number_to_obtain == 0);
  const auto& os = TestSettings::get().os;
  solve_problem(pdata, tdata, info, solver_params, stats, *this, os);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
