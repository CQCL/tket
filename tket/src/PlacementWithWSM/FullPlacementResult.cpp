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

#include "PlacementWithWSM/FullPlacementResult.hpp"

#include <algorithm>
#include <sstream>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SpecialExceptions.hpp"
#include "WeightSubgrMono/EndToEndWrappers/MainSolver.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

static GraphEdgeWeights get_complete_target_graph(
    const GraphEdgeWeights& enlarged_target_graph) {
  WeightWSM max_edge_weight = 0;
  for (const auto& entry : enlarged_target_graph) {
    max_edge_weight = std::max(max_edge_weight, entry.second);
  }
  max_edge_weight = get_product_or_throw(max_edge_weight, WeightWSM(3));
  const auto t_vertices = get_vertices(enlarged_target_graph);
  auto complete_target_graph = enlarged_target_graph;
  for (auto citer1 = t_vertices.cbegin(); citer1 != t_vertices.cend();
       ++citer1) {
    for (auto citer2 = citer1 + 1; citer2 != t_vertices.cend(); ++citer2) {
      const auto t_edge = get_edge(*citer1, *citer2);
      if (complete_target_graph.count(t_edge) == 0) {
        complete_target_graph[t_edge] = max_edge_weight;
      }
    }
  }
  TKET_ASSERT(
      complete_target_graph.size() ==
      (t_vertices.size() * (t_vertices.size() - 1)) / 2);
  return complete_target_graph;
}

static FullPlacementResult get_single_pass_result(
    const GraphEdgeWeights& pattern_graph,
    const GraphEdgeWeights& original_target_graph,
    // "enlarged_target_graph" should be original_target_graph
    // with extra edges and weights added.
    const GraphEdgeWeights& enlarged_target_graph,
    // The actual target graph to use; may be different from the above.
    const GraphEdgeWeights& target_graph_to_use,
    const std::vector<std::set<VertexWSM>>& gates, unsigned timeout1_ms,
    unsigned timeout2_ms, bool& final_solution_is_complete,
    std::optional<std::size_t> max_iterations_opt) {
  MainSolver solver;

  // Note that "statistics" and "best_solution" are stored within the MainSolver
  // and remain valid; although we cannot alter them, the MainSolver object can.
  const auto& statistics =
      solver.initialise(pattern_graph, target_graph_to_use);

  TKET_ASSERT(statistics.search_time_ms == 0);
  if (statistics.initialisation_time_ms >= timeout1_ms) {
    std::stringstream ss;
    ss << "Initialisation took " << statistics.initialisation_time_ms
       << " ms, already longer than timeout " << timeout1_ms << " ms.";
    throw InitialisationTimeout(ss.str());
  }

  // Continue with the solve.
  long long remaining_time = timeout1_ms - statistics.initialisation_time_ms;
  MainSolver::Parameters solver_params(remaining_time);
  if (max_iterations_opt) {
    solver_params.max_iterations = max_iterations_opt.value();
  }
  solver.solve(solver_params);
  const auto& best_solution = solver.get_best_solution();

  if (best_solution.complete && timeout2_ms > 0) {
    // We have at least one full solution,
    // so we use up the rest of the timeout time
    // to make it even better. (Note that if it's already finished,
    // it won't compute any further; so no time is wasted).

    // Note that if we DON'T have a full solution yet,
    // the caller is going to try something else.
    const auto total_time_so_far =
        statistics.initialisation_time_ms + statistics.search_time_ms;
    const auto allowed_total_time = timeout1_ms + timeout2_ms;
    if (total_time_so_far < allowed_total_time) {
      solver_params.timeout_ms = allowed_total_time - total_time_so_far;
      solver.solve(solver_params);
    }
  }

  FullPlacementResult full_result;
  full_result.result = PlacementAndStatistics(
      pattern_graph, original_target_graph, enlarged_target_graph, gates,
      best_solution);

  full_result.iterations_for_pass = statistics.iterations;
  full_result.total_init_time_ms = statistics.initialisation_time_ms;
  full_result.total_search_time_ms = statistics.search_time_ms;

  final_solution_is_complete =
      best_solution.complete && best_solution.assignments.size() ==
                                    full_result.result.valid_assignments.size();
  return full_result;
}

FullPlacementResult::FullPlacementResult() {}

FullPlacementResult::FullPlacementResult(
    const GraphEdgeWeights& pattern_graph,
    const GraphEdgeWeights& original_target_graph,
    // "enlarged_target_graph" should be original_target_graph
    // with extra edges and weights added.
    const GraphEdgeWeights& enlarged_target_graph,
    const std::vector<std::set<VertexWSM>>& gates,
    const Parameters& parameters) {
  // The caller should never be using such stupid short timeouts anyway.
  const unsigned min_timeout_ms = 4;
  const unsigned total_timeout_ms =
      std::max(parameters.timeout_ms, min_timeout_ms);
  bool final_solution_is_complete;

  if (parameters.pass_data_opt) {
    const auto pass_data = parameters.pass_data_opt.value();
    if (pass_data.first == Pass::INITIAL) {
      // We'll use all our time on the original solve.
      *this = get_single_pass_result(
          pattern_graph, original_target_graph, enlarged_target_graph,
          enlarged_target_graph, gates, total_timeout_ms, 0,
          final_solution_is_complete, pass_data.second);
    } else {
      TKET_ASSERT(pass_data.first == Pass::COMPLETE_TARGET_GRAPH);
      const auto complete_target_graph =
          get_complete_target_graph(enlarged_target_graph);
      *this = get_single_pass_result(
          pattern_graph, original_target_graph, enlarged_target_graph,
          complete_target_graph, gates, total_timeout_ms, 0,
          final_solution_is_complete, pass_data.second);
    }
    pass = pass_data.first;
    number_of_passes = 1;
    return;
  }

  // We're doing the complete collection of passes.
  const unsigned timeout1_ms = total_timeout_ms / 4;
  const unsigned timeout2_ms = total_timeout_ms - timeout1_ms;
  *this = get_single_pass_result(
      pattern_graph, original_target_graph, enlarged_target_graph,
      enlarged_target_graph, gates, timeout1_ms, timeout2_ms,
      final_solution_is_complete, parameters.max_iterations_opt);

  pass = Pass::INITIAL;
  number_of_passes = 1;

  if (!final_solution_is_complete) {
    const auto total_time_so_far = total_init_time_ms + total_search_time_ms;
    if (total_time_so_far + min_timeout_ms < total_timeout_ms) {
      // Solving with the original target graph failed to find a complete
      // solution, yet there is still some time budget left. Therefore, spend
      // all remaining time with the COMPLETE target graph, for which full
      // solutions are GUARANTEED, even if poor.
      const auto complete_target_graph =
          get_complete_target_graph(enlarged_target_graph);
      const auto new_solution = get_single_pass_result(
          pattern_graph, original_target_graph, enlarged_target_graph,
          complete_target_graph, gates, total_timeout_ms - total_time_so_far, 0,
          final_solution_is_complete, parameters.max_iterations_opt);
      if (result.prefer_other_solution(new_solution.result)) {
        // We're going to switch over!
        result = new_solution.result;
        pass = Pass::COMPLETE_TARGET_GRAPH;
        iterations_for_pass = new_solution.iterations_for_pass;
      }
      total_init_time_ms += new_solution.total_init_time_ms;
      total_search_time_ms += new_solution.total_search_time_ms;
      ++number_of_passes;
    }
  }
}

std::string FullPlacementResult::str(bool print_times) const {
  std::stringstream ss;
  ss << result.str() << "\nPasses: " << number_of_passes << "; best: ";
  switch (pass) {
    case Pass::INITIAL:
      ss << "INITIAL";
      break;
    case Pass::COMPLETE_TARGET_GRAPH:
      ss << "COMPLETE_TARGET_GRAPH";
      break;
    default:
      ss << "UNKNOWN! ERROR!";
  }
  ss << "; iterations: " << iterations_for_pass;
  if (print_times) {
    ss << "\nTotal time: " << total_init_time_ms << "+" << total_search_time_ms
       << " = " << total_init_time_ms + total_search_time_ms;
  }
  return ss.str();
}

}  // namespace tket
