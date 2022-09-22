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

#include <tkwsm/EndToEndWrappers/MainSolver.hpp>

#include "Placement/Placement.hpp"
#include "RelabelledGraphWSM.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

using RelabelledPatternGraph =
    RelabelledGraphWSM<Qubit, QubitGraph::UndirectedConnGraph>;
using RelabelledTargetGraph =
    RelabelledGraphWSM<Node, Architecture::UndirectedConnGraph>;
using BimapValue = qubit_bimap_t::value_type;

// Where should isolated pattern vertices be assigned?
// They might NOT have been isolated originally; it may be
// that we deliberately erased some pattern edges.
// Thus, we still want them connected to useful target components,
// so assign to nonisolated target vertices first.
static void assign_isolated_pattern_vertices(
    qubit_bimap_t& map, const RelabelledPatternGraph& relabelled_pattern_graph,
    const RelabelledTargetGraph& relabelled_target_graph) {
  if (relabelled_pattern_graph.get_relabelled_isolated_vertices().empty()) {
    return;
  }
  // The PV assigned so far must be exactly the nonisolated ones.
  TKET_ASSERT(
      map.size() ==
      relabelled_pattern_graph.get_relabelled_nonisolated_vertices().size());

  // Also, all PV so far must have been assigned to nonisolated TV.
  TKET_ASSERT(
      map.size() <=
      relabelled_target_graph.get_relabelled_nonisolated_vertices().size());

  std::set<VertexWSM> unused_target_vertices;
  bool refilled_with_isolated_tv;
  if (map.size() <
      relabelled_target_graph.get_relabelled_nonisolated_vertices().size()) {
    refilled_with_isolated_tv = false;
    unused_target_vertices =
        relabelled_target_graph.get_relabelled_nonisolated_vertices();
    for (const auto& entry : map.left) {
      // We CANNOT have assigned a nonisolated pattern vertex
      // to an isolated target vertex.
      const VertexWSM new_pv =
          relabelled_pattern_graph.get_relabelled_vertex(entry.first);
      const VertexWSM new_tv =
          relabelled_target_graph.get_relabelled_vertex(entry.second);
      TKET_ASSERT(
          relabelled_pattern_graph.get_relabelled_nonisolated_vertices().count(
              new_pv) != 0);
      TKET_ASSERT(
          relabelled_target_graph.get_relabelled_nonisolated_vertices().count(
              new_tv) != 0);
      TKET_ASSERT(unused_target_vertices.erase(new_tv) == 1);
    }
  } else {
    refilled_with_isolated_tv = true;
    unused_target_vertices =
        relabelled_target_graph.get_relabelled_isolated_vertices();
  }
  for (VertexWSM isolated_pv :
       relabelled_pattern_graph.get_relabelled_isolated_vertices()) {
    if (unused_target_vertices.empty()) {
      // We cannot run out of vertices.
      TKET_ASSERT(!refilled_with_isolated_tv);
      refilled_with_isolated_tv = true;
      unused_target_vertices =
          relabelled_target_graph.get_relabelled_isolated_vertices();
      TKET_ASSERT(!unused_target_vertices.empty());
    }
    const VertexWSM next_tv = *unused_target_vertices.cbegin();
    unused_target_vertices.erase(next_tv);

    map.insert(BimapValue(
        relabelled_pattern_graph.get_original_vertices().at(isolated_pv),
        relabelled_target_graph.get_original_vertices().at(next_tv)));
  }
}

static void write_solver_solutions(
    std::vector<qubit_bimap_t>& all_maps,
    const std::vector<SolutionWSM>& solutions,
    const RelabelledPatternGraph& relabelled_pattern_graph,
    const RelabelledTargetGraph& relabelled_target_graph) {
  TKET_ASSERT(all_maps.empty());
  all_maps.resize(solutions.size());
  const WeightWSM expected_weight =
      relabelled_pattern_graph.get_relabelled_edges_and_weights().size();

  for (unsigned ii = 0; ii < solutions.size(); ++ii) {
    const auto& solution = solutions[ii];
    auto& map = all_maps[ii];
    TKET_ASSERT(solution.scalar_product == expected_weight);
    TKET_ASSERT(solution.total_p_edges_weight == expected_weight);
    for (const auto& relabelled_pv_tv : solution.assignments) {
      map.insert(BimapValue(
          relabelled_pattern_graph.get_original_vertices().at(
              relabelled_pv_tv.first),
          relabelled_target_graph.get_original_vertices().at(
              relabelled_pv_tv.second)));
    }
    assign_isolated_pattern_vertices(
        map, relabelled_pattern_graph, relabelled_target_graph);
    TKET_ASSERT(
        map.size() == relabelled_pattern_graph.get_original_vertices().size());
  }
}

/**
 * \cond Somehow doxygen 1.9.1 complains about this. Tell it to be quiet.
 */

std::vector<qubit_bimap_t> get_unweighted_subgraph_monomorphisms(
    const QubitGraph::UndirectedConnGraph& pattern_graph,
    const Architecture::UndirectedConnGraph& target_graph, unsigned max_matches,
    unsigned timeout_ms) {
  std::vector<qubit_bimap_t> all_maps;

  const RelabelledPatternGraph relabelled_pattern_graph(pattern_graph);
  const RelabelledTargetGraph relabelled_target_graph(target_graph);

  if (relabelled_pattern_graph.get_relabelled_edges_and_weights().size() >
          relabelled_target_graph.get_relabelled_edges_and_weights().size() ||
      relabelled_pattern_graph.get_relabelled_nonisolated_vertices().size() >
          relabelled_target_graph.get_relabelled_nonisolated_vertices()
              .size() ||

      relabelled_pattern_graph.get_relabelled_isolated_vertices().size() +
              relabelled_pattern_graph.get_relabelled_nonisolated_vertices()
                  .size() >
          relabelled_target_graph.get_relabelled_isolated_vertices().size() +
              relabelled_target_graph.get_relabelled_nonisolated_vertices()
                  .size()) {
    // The problem is trivially insoluble.
    return all_maps;
  }

  if (relabelled_pattern_graph.get_relabelled_nonisolated_vertices().empty()) {
    // Trivial special case: all pattern vertices are isolated!
    all_maps.resize(1);
    assign_isolated_pattern_vertices(
        all_maps[0], relabelled_pattern_graph, relabelled_target_graph);
    return all_maps;
  }

  MainSolverParameters solver_parameters;
  solver_parameters.terminate_with_first_full_solution = false;
  solver_parameters.for_multiple_full_solutions_the_max_number_to_obtain =
      max_matches;
  solver_parameters.timeout_ms = timeout_ms;

  const MainSolver main_solver(
      relabelled_pattern_graph.get_relabelled_edges_and_weights(),
      relabelled_target_graph.get_relabelled_edges_and_weights(),
      solver_parameters);

  const auto& solution_data = main_solver.get_solution_data();

  write_solver_solutions(
      all_maps, solution_data.solutions, relabelled_pattern_graph,
      relabelled_target_graph);

  return all_maps;
}

/**
 * \endcond
 */

}  // namespace tket
