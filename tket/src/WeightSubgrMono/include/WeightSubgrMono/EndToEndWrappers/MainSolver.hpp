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
#include <memory>
#include <optional>

#include "../Searching/SolutionWSM.hpp"
#include "SolutionStatistics.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** The main class which takes a raw WSM problem and tries to solve it. */
class MainSolver {
 public:
  MainSolver();

  ~MainSolver();

  // NOTE: all the returned const SolutionStatistics& actually refer
  // to a single internally stored object, so the same reference
  // remains valid and can be reused throughout the lifetime of this object.

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

  /** Extra parameters to configure the algorithm. */
  struct Parameters {
    /** The timeout only for the current call to "solve",
     * not the cumulative time.
     */
    long long timeout_ms;

    /** If set to TRUE, it will break off early when a complete solution
     * is found, no matter how poor (but still count as UNFINISHED,
     * in case the caller wants to continue looking for more solutions).
     */
    bool terminate_with_first_full_solution;

    /** For internal testing; set an upper bound on
     * the number of search iterations.
     */
    std::size_t max_iterations;

    /** If non-null, it means that we artificially impose
     * an upper bound on the total weight for a full solution,
     * i.e. an extra constraint on the solution.
     * This should speed up the solving, BUT may prevent us
     * from finding any solution.
     */
    std::optional<WeightWSM> weight_upper_bound_constraint;

    Parameters();

    /** Just set the timeout in milliseconds; the most common parameter. */
    explicit Parameters(long long timeout_ms);
  };

  /** After "initialise" has been called, actually run the solver.
   * This can be done repeatedly. The statistics are updated
   * with the cumulative results (time, etc.) over all runs.
   * If the solver has already finished (i.e., gone through the
   * complete search space, EITHER finding the optimal solution OR proving
   * that no solution exists), does nothing and wastes no time.
   */
  const SolutionStatistics& solve(const Parameters& parameters = {});

  /** Run the solver with a specified timeout in milliseconds. */
  const SolutionStatistics& solve(long long timeout_ms);

  /** All-in-one initialise and solve function. */
  const SolutionStatistics& solve(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges, const Parameters& parameters = {});

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
  // The pimpl idiom, to prevent big header files...
  struct Impl;
  std::unique_ptr<Impl> m_pimpl;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
