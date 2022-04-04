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

  /** Upon construction, try to solve the problem.
   * @param pattern_edges The pattern graph, with edge weights
   * @param target_edges The target graph, with edge weights
   * @param parameters The parameters which configure the solving algorithm.
   */
  MainSolver(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges,
      const MainSolverParameters& parameters);

  /** Upon construction, do a single solve iteration, trying to use as many of the suggested assignments as possible.
   * Partial, incomplete, and even inconsistent suggestions will still be accepted;   
   * they will not invalidate the solver, but could lead to slowdowns (if the suggestion is very bad!)
   * The function "solve" will then behave in future exactly as though this initial move
   * had arisen normally from the variable and value ordering heuristics.
   * @param pattern_edges The pattern graph, with edge weights
   * @param target_edges The target graph, with edge weights
   * @param suggested_assignments A list of assignments to try, which need not be complete or valid.
   */
  MainSolver(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges,
      const std::vector<std::pair<VertexWSM, VertexWSM>>&
          suggested_assignments);

  /** One can continue solving, if the original solve terminated early
   * (which may occur before the end, due to timeouts, maximum iterations being hit, etc.)
   * Does nothing if it finished already (i.e., has gone through a complete search, so can guarantee that EITHER it has found
   * a jointly optimal solution, OR that no solutions exist).
   * Note that the timeout refers ONLY to the extra time taken by this function, not the total elapsed time so far which is stored and updated in the SolutionStatistics object within this class.
   * @param parameters The parameters which configure the solving algorithm.
   */
  void solve(const MainSolverParameters& parameters);

  const SolutionStatistics& get_solution_statistics() const;

  /** Return the best solution found so far. It may change as "solve" is
   * called more times, but the object is stored internally in this class
   * so the reference remains valid as long as this class does.
   */
  const SolutionWSM& get_best_solution() const;

  typedef std::vector<std::vector<std::pair<VertexWSM, VertexWSM>>>
      FullSolutionsList;

  /** A wrapper around the function in the internal SolutionStorage class.
   * See that class for more explanation. (This can be tricky. For example,
   * this might be empty, even if some full solutions were found;
   * it depends on the input parameters. In the unweighted case,
   * if it runs to completion and the upper bound on the number of solutions
   * is not hit, then this is guaranteed to contain ALL solutions; but in the
   * weighted case, it is NOT, even though it must contain at least one
   * jointly optimal solution).
   * @return The solutions returned by the function of the same name in the
   * internal SolutionStorage class.
   */
  const FullSolutionsList& get_some_full_solutions() const;

 private:
  MainSolverData m_data;

  /** Do not actually solve, just set up the initial data etc. READY to solve.
   * Note that the GraphEdgeWeights objects are not needed after this,
   * as the data has been processed and converted.
   * @param pattern_edges The pattern graph, with edge weights
   * @param target_edges The target graph, with edge weights
   */
  MainSolver(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges);

  /** Before calling the "solve" function,
   * take the suggested assignments and move down the search branch once,
   * trying to set as many assignments as possible.
   * Solve will behave in future as though this initial move
   * had arisen normally from the variable and value ordering heuristics.
   * @param suggested_assignments A list of assignments to try, which need not be complete or valid.
   */
  void do_one_solve_iteration_with_suggestion(
      const std::vector<std::pair<VertexWSM, VertexWSM>>&
          suggested_assignments);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
