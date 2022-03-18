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
#include "MainSolverParameters.hpp"
#include "SolutionStatistics.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/EndToEndWrappers/CheckedWeightBounds.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchBranch.hpp"
#include "WeightSubgrMono/Searching/SharedData.hpp"
#include "WeightSubgrMono/Searching/ValueOrdering.hpp"
#include "WeightSubgrMono/Searching/VariableOrdering.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Internal data stored within the MainSolver object,
 * separate from the input parameters, to be used during the search. */
struct MainSolverData {
  bool initialised = false;
  FixedData fixed_data;
  std::unique_ptr<SharedData> shared_data_ptr;
  SearchBranch branch;
  VariableOrdering var_ordering;
  ValueOrdering val_ordering;
  SolutionStatistics statistics;
  std::optional<WeightWSM> previous_upper_bound_constraint;

  // Hack: in case we finish upon initialisation; an empty solution.
  // Otherwise, use the one stored in shared_data_ptr.
  const SolutionWSM empty_solution;

  /** Sets "fixed_data" and, if successful, shared_data_ptr.
   * Does some basic reduction also: some problems can be solved, or shown to be
   * impossible, without any searching.
   * @param pattern_edges The raw data for the pattern graph.
   * @param target_edges The raw data for the target graph.
   * @return The status: whether the problem is trivially soluble or insoluble
   * without deep searching, or still needs searching.
   */
  ReductionResult initialise(
      const GraphEdgeWeights& pattern_edges,
      const GraphEdgeWeights& target_edges);

  /** The main searching loop to be called, after it has been initialised. */
  void solve_loop_after_initialisation(const MainSolverParameters& parameters);

  void do_one_solve_iteration_with_suggestion(
      const std::vector<std::pair<VertexWSM, VertexWSM>>&
          suggested_assignments);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
