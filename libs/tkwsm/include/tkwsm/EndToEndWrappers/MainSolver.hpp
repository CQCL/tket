// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <chrono>
#include <memory>

#include "tkwsm/EndToEndWrappers/MainSolverParameters.hpp"
#include "tkwsm/EndToEndWrappers/SolutionData.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/GraphTheoretic/VertexRelabelling.hpp"
#include "tkwsm/Searching/SearchBranch.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct PreSearchComponents;
struct SearchComponents;

/** The main class which takes a raw WSM problem and tries to solve it. */
class MainSolver {
 public:
  /** Upon construction, try to solve the problem.
   * @param pattern_edges The pattern graph, with edge weights
   * @param target_edges The target graph, with edge weights
   * @param parameters Parameters which configure the solving algorithm.
   */
  MainSolver(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges,
      const MainSolverParameters& parameters);

  ~MainSolver();

  /** After construction, do further solving, if the original solve terminated
   * early (which may occur due to timeouts, maximum iterations being hit, etc.)
   * Does nothing if it finished already (i.e., has gone through a complete
   * search, so can guarantee that EITHER it has found a jointly optimal
   * solution, OR that no solutions exist). Note that the timeout refers ONLY to
   * the extra time taken by this function, not the total elapsed time so far
   * which is stored and updated in the SolutionStatistics object within this
   * class.
   * @param parameters The parameters which configure the solving algorithm.
   */
  void solve(const MainSolverParameters& parameters);

  /** Returns information about the solving and the best solution
   * found so far. */
  const SolutionData& get_solution_data() const;

 private:
  const VertexRelabelling m_pattern_vertex_relabelling;
  const VertexRelabelling m_target_vertex_relabelling;

  NeighboursData m_pattern_neighbours_data;
  NeighboursData m_target_neighbours_data;

  SolutionData m_solution_data;
  mutable SolutionData m_solution_data_original_vertices;

  // If the problem is trivially insoluble, no need to spend time constructing
  // these.
  std::unique_ptr<PreSearchComponents> m_pre_search_components_ptr;
  std::unique_ptr<SearchComponents> m_search_components_ptr;
  std::unique_ptr<SearchBranch> m_search_branch_ptr;

  /** Do NOT backtrack, just move down directly from the current node as far as
   * possible, i.e. a single solve iteration. Returns TRUE if we end with a full
   * solution, false otherwise.
   */
  bool move_down_from_reduced_node(
      const SearchBranch::ReductionParameters& reduction_parameters);

  /** Performs the solve.
   * We should NOT time things by timing each individual iteration and summing
   * them; instead, we set the END time and stop when we go over. This converts
   * the max EXTRA iterations and time into an absolute maximum time and
   * iterations.
   */
  void internal_solve(
      const MainSolverParameters& parameters, std::size_t max_iterations,
      const std::chrono::steady_clock::time_point& desired_end_time);

  void add_solution_from_final_node(
      const MainSolverParameters& parameters,
      const SearchBranch::ReductionParameters& reduction_parameters);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
