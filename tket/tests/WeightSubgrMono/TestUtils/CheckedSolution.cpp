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

#include "CheckedSolution.hpp"

#include <catch2/catch.hpp>

#include "../TestUtils/TestSettings.hpp"
#include "WeightSubgrMono/EndToEndWrappers/MainSolver.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static void check_finished_complete_solution(
    const CheckedSolution::ProblemInformation& info,
    const SolutionStatistics& soln_statistics, const SolutionWSM& solution,
    CheckedSolution::Statistics& stats, const OStreamWrapper& os) {
  // Finished and complete; this MUST be the optimal solution.
  if (!info.known_optimal_solution) {
    os << "; soln " << solution.total_scalar_product_weight;
  }
  os << "; " << soln_statistics.iterations << " iters.";
  bool satisfies_optimal_bounds = true;
  if (info.known_lower_bound) {
    const auto lb = info.known_lower_bound.value();
    if (lb > solution.total_scalar_product_weight) {
      satisfies_optimal_bounds = false;
    }
  }
  if (info.known_upper_bound) {
    const auto ub = info.known_upper_bound.value();
    if (ub < solution.total_scalar_product_weight) {
      satisfies_optimal_bounds = false;
    }
  }
  if (satisfies_optimal_bounds) {
    ++stats.success_count;
  } else {
    ++stats.failure_count;
    os << " violates known soln bounds!";
  }
}

static void check_unfinished_solution(
    const CheckedSolution::ProblemInformation& info,
    const SolutionStatistics& soln_statistics, const SolutionWSM& solution,
    const MainSolverParameters& solver_params,
    CheckedSolution::Statistics& stats, const OStreamWrapper& os) {
  if (soln_statistics.iterations >= solver_params.max_iterations) {
    os << " - hit iterations limit: " << soln_statistics.iterations;
  } else {
    os << " - TIMED OUT.";
  }

  ++stats.timeout_count;
  if (!solution.complete) {
    // Not really anything useful to test for an incomplete solution.
    return;
  }
  // Just check the lower bound, nothing else we can do.
  // (It was already checked to be VALID, with get_errors).
  if (info.known_lower_bound) {
    const auto lb = info.known_lower_bound.value();
    if (lb > solution.total_scalar_product_weight) {
      os << " - error - we know that the optimal solution value (if it exists) "
            "has lower bound "
         << lb << ". But we found a complete solution with value "
         << solution.total_scalar_product_weight;

      ++stats.failure_count;
    }
  }
  if (info.existence == CheckedSolution::ProblemInformation::
                            SolutionsExistence::KNOWN_TO_BE_INSOLUBLE) {
    os << " - error - found a complete solution with value "
       << solution.total_scalar_product_weight
       << ", even though we know no solution exists";

    ++stats.failure_count;
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

static void solve_problem(
    const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
    CheckedSolution::ProblemInformation info,
    const MainSolverParameters& solver_params,
    CheckedSolution::Statistics& stats, CheckedSolution& checked_solution,
    const OStreamWrapper& os,
    const std::vector<std::pair<VertexWSM, VertexWSM>>& suggested_assignments) {
  check_known_solution_information(info);
  check_for_impossible_weight_constraint(info, solver_params);
  MainSolver solver;
  if(suggested_assignments.empty()) {
    solver.solve(pdata, tdata, solver_params);
  } else {
    solver.initialise(pdata, tdata);
    solver.do_one_solve_iteration_with_suggestion(suggested_assignments);
    solver.solve(solver_params);
  }
  const auto& solution = solver.get_best_solution();
  const auto& soln_statistics = solver.get_solution_statistics();
  checked_solution.iterations = soln_statistics.iterations;
  checked_solution.finished = soln_statistics.finished;
  checked_solution.assignments = solution.assignments;
  if (solution.complete) {
    checked_solution.complete_solution_weight =
        solution.total_scalar_product_weight;
  }

  stats.total_init_time_ms += soln_statistics.initialisation_time_ms;
  stats.total_search_time_ms += soln_statistics.search_time_ms;

  os << " - time " << soln_statistics.initialisation_time_ms << "+"
     << soln_statistics.search_time_ms;
  if (info.known_optimal_solution) {
    os << " - known opt.val. " << info.known_optimal_solution.value();
  }

  {
    const auto errors = solution.get_errors(pdata, tdata);
    CHECK("" == errors);
    if (!errors.empty()) {
      ++stats.failure_count;
      return;
    }
    // The solution, if complete, is valid.
  }

  if (soln_statistics.finished) {
    if (solution.complete) {
      check_finished_complete_solution(
          info, soln_statistics, solution, stats, os);
      return;
    }
    // It is finished, but with no COMPLETE solution.
    os << " - no soln.";

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
    return;
  }
  if(solver_params.terminate_with_first_full_solution && solution.complete) {
    check_finished_complete_solution(
          info, soln_statistics, solution, stats, os);
    return;
  }
  check_unfinished_solution(
      info, soln_statistics, solution, solver_params, stats, os);
}


CheckedSolution::CheckedSolution(
    const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
    ProblemInformation info, const MainSolverParameters& solver_params,
    Statistics& stats,
    const std::vector<std::pair<VertexWSM, VertexWSM>>& suggested_assignments) {
  const auto& os = TestSettings::get().os;

  const auto orig_init_time = stats.total_init_time_ms;
  const auto orig_search_time = stats.total_search_time_ms;
  solve_problem(pdata, tdata, info, solver_params, stats, *this, os, suggested_assignments);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
