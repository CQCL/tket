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
#include "MainSolverData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** The main class which takes a raw WSM problem and tries to solve it. */
class MainSolver {
 public:
  // NOTE: all the returned const SolutionStatistics& actually refer
  // to a single internally stored object, so the same reference
  // remains valid and can be reused throughout the lifetime of this object.

  const SolutionStatistics& get_solution_statistics() const;

  /** Do not actually solve, just set up the initial data etc. READY to solve.
   * Note that the GraphEdgeWeights objects are not needed after this,
   * as the data has been processed and converted.
   * If any previous problem was underway,
   * clears all the data and resets for this new problem.
   * @param pattern_edges The pattern graph, with edge weights
   * @param target_edges The target graph, with edge weights
   * @return Information about the solution so far (it's possible that the
   * problem is easy enough to solve in the initialisation stage, without any
   * search).
   */
  const SolutionStatistics& initialise(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges);

  /** Immediately after "initialise", before calling any "solve" function:
   * Take the suggested assignments and move down the search branch once,
   * trying to set as many assignments as possible.
   * Solve will behave in future as though this initial move
   * had arisen normally from the variable and value ordering heuristics.
   */
  void do_one_solve_iteration_with_suggestion(
      const std::vector<std::pair<VertexWSM, VertexWSM>>&
          suggested_assignments);

  /** After "initialise" has been called, actually run the solver.
   * This can be done repeatedly. The statistics are updated
   * with the cumulative results (time, etc.) over all runs.
   * If the solver has already finished (i.e., gone through the
   * complete search space, EITHER finding the optimal solution OR proving
   * that no solution exists), does nothing and wastes no time.
   */
  const SolutionStatistics& solve(const MainSolverParameters& parameters = {});

  /** Run the solver with a specified timeout in milliseconds. */
  const SolutionStatistics& solve(long long timeout_ms);

  /** All-in-one initialise and solve function. */
  const SolutionStatistics& solve(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges,
      const MainSolverParameters& parameters = {});

  /** All-in-one initialise and solve function. */
  const SolutionStatistics& solve(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges, long long timeout_ms);

  /** Return the best solution found so far. It may change as "solve" is
   * called more times, but the object is stored internally in this class
   * so the reference remains valid as long as this class does.
   */
  const SolutionWSM& get_best_solution() const;

 private:
  MainSolverData m_data;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
