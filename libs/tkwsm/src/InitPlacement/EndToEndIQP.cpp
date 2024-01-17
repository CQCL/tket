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

#include "tkwsm/InitPlacement/EndToEndIQP.hpp"

#include <chrono>
#include <sstream>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/EndToEndWrappers/MainSolver.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/GraphTheoretic/VertexRelabelling.hpp"
#include "tkwsm/InitPlacement/InputStructs.hpp"
#include "tkwsm/InitPlacement/MonteCarloCompleteTargetSolution.hpp"
#include "tkwsm/InitPlacement/PrunedTargetEdges.hpp"
#include "tkwsm/InitPlacement/UtilsIQP.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

typedef std::chrono::steady_clock Clock;

namespace {

// Gets an initial placement using MCCT (Monte Carlo complete target graph),
// and also contains intermediate
// objects needed for other calculations.
struct MCCTData {
  MCCTData(
      const GraphEdgeWeights& pattern_graph_weights,
      const GraphEdgeWeights& target_architecture_with_error_weights,
      WeightWSM& scalar_product, unsigned& number_of_iterations)
      : pattern_relabelling(pattern_graph_weights),
        relabelled_pattern_ndata(pattern_relabelling.new_edges_and_weights),
        target_relabelling(target_architecture_with_error_weights),

        expanded_target_graph_data(target_architecture_with_error_weights),
        relabelled_explicit_target_ndata(get_relabelled_graph_data(
            expanded_target_graph_data.explicit_target_graph_weights,
            target_relabelling))

  {
    const MonteCarloCompleteTargetSolution mcct_solution(
        relabelled_pattern_ndata, relabelled_explicit_target_ndata,
        expanded_target_graph_data.implicit_weight);

    number_of_iterations = mcct_solution.iterations();
    scalar_product = mcct_solution.get_best_scalar_product();
    new_label_assignments = mcct_solution.get_best_assignments();
  }

  // The below data is needed later, for target pruning
  // (MCCT, of course, does NOT use target edge pruning).

  // We need to relabel the pattern vertices for neighbours data.
  // TODO refactor: MainSolver also does this internally,
  // which is inelegent wasted duplication
  // (although completely negligible impact).
  const VertexRelabelling pattern_relabelling;
  const NeighboursData relabelled_pattern_ndata;

  // We also need to relabel the target vertices for neighbours data.
  /// TODO Again, refactor: MainSolver also does this internally.
  // Note that this includes ALL TV (although relabelled).
  // This is the same relabelling for both the original and expanded
  // target graphs (as expansion only adds edges, not vertices).
  const VertexRelabelling target_relabelling;

  const TargetGraphData expanded_target_graph_data;

  // We've expanded the target; this contains all the explicit weights.
  const NeighboursData relabelled_explicit_target_ndata;

  // Element[PV] = TV, using the relabelled PV, TV.
  std::vector<unsigned> new_label_assignments;
};

}  // namespace

IQPResult::IQPResult(
    const GraphEdgeWeights& pattern_graph_weights,
    const GraphEdgeWeights& target_architecture_with_error_weights,
    unsigned timeout_ms, const IQPParameters& iqp_parameters) {
  const auto start = Clock::now();
  const MCCTData mcct_data(
      pattern_graph_weights, target_architecture_with_error_weights,
      mcct_scalar_product, mcct_iterations);

  mcct_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                     Clock::now() - start)
                     .count();

  TKET_ASSERT(
      mcct_data.new_label_assignments.size() ==
      mcct_data.pattern_relabelling.number_of_vertices);

  initial_qubit_placement.reserve(
      mcct_data.pattern_relabelling.number_of_vertices);
  for (unsigned new_pv = 0;
       new_pv < mcct_data.pattern_relabelling.number_of_vertices; ++new_pv) {
    initial_qubit_placement.emplace_back(
        mcct_data.pattern_relabelling.get_old_label(new_pv),
        mcct_data.target_relabelling.get_old_label(
            mcct_data.new_label_assignments[new_pv]));
  }

  if (mcct_time_ms >= timeout_ms || mcct_scalar_product == 0) {
    // (Unless we've got many zero weights, which is extremely unusual):
    // we're out of time - no WSM.
    total_time_ms = mcct_time_ms;
    std::sort(initial_qubit_placement.begin(), initial_qubit_placement.end());
    return;
  }

  // Prune the complete target edges.
  // This uses the same new TV labels as in MCCT.
  // But beware: some new TV might not appear in this data,
  // because we've deliberately cut it down to make WSM solving fast
  // (otherwise, we'd just retain the complete target graph).

  const GraphEdgeWeights relabelled_pruned_wsm_target_graph_data =
      get_new_target_graph_data(
          mcct_data.relabelled_pattern_ndata,
          mcct_data.relabelled_explicit_target_ndata,
          mcct_data.expanded_target_graph_data.implicit_weight,
          mcct_data.new_label_assignments);

  wsm_number_of_pruned_tv =
      get_number_of_vertices(relabelled_pruned_wsm_target_graph_data);
  wsm_number_of_pruned_t_edges = relabelled_pruned_wsm_target_graph_data.size();

  // Be paranoid and recalculate the scalar product...
  // this can be removed when it's been thoroughly tested
  // (although it does no harm - negligible extra time)
  {
    auto old_pv_new_tv_assignments = initial_qubit_placement;
    for (auto& entry : old_pv_new_tv_assignments) {
      entry.second = mcct_data.target_relabelling.get_new_label(entry.second);
    }
    TKET_ASSERT(
        mcct_scalar_product == get_checked_scalar_product(
                                   pattern_graph_weights,
                                   relabelled_pruned_wsm_target_graph_data,
                                   old_pv_new_tv_assignments));
  }
  MainSolverParameters wsm_solver_parameters;

  // We already checked that there is positive time remaining.
  wsm_solver_parameters.timeout_ms = timeout_ms - mcct_time_ms;
  wsm_solver_parameters.iterations_timeout = iqp_parameters.max_wsm_iterations;

  // We must strictly improve on our MCCT solution, otherwise no point.
  wsm_solver_parameters.weight_upper_bound_constraint = mcct_scalar_product - 1;

  const MainSolver wsm_solver(
      pattern_graph_weights, relabelled_pruned_wsm_target_graph_data,
      wsm_solver_parameters);
  const SolutionData& wsm_solution = wsm_solver.get_solution_data();

  wsm_init_time_ms = wsm_solution.initialisation_time_ms;
  wsm_solve_time_ms = wsm_solution.search_time_ms;
  total_time_ms = wsm_init_time_ms + wsm_solve_time_ms + mcct_time_ms;
  wsm_iterations = wsm_solution.iterations;
  if (wsm_solution.solutions.empty()) {
    return;
  }
  TKET_ASSERT(wsm_solution.solutions.size() == 1);
  TKET_ASSERT(wsm_solution.solutions[0].scalar_product <= mcct_scalar_product);

  wsm_scalar_product_opt = wsm_solution.solutions[0].scalar_product;
  initial_qubit_placement = wsm_solution.solutions[0].assignments;

  // Care: the WSM problem used new TV labels.
  for (auto& entry : initial_qubit_placement) {
    entry.second = mcct_data.target_relabelling.get_old_label(entry.second);
  }
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
