// Copyright Quantinuum
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

#include "tket/Architecture/SubgraphMonomorphisms.hpp"

#include <set>
#include <tkassert/Assert.hpp>
#include <tkwsm/EndToEndWrappers/MainSolver.hpp>

namespace tket {

using namespace WeightedSubgraphMonomorphism;

static GraphEdgeWeights get_weight_one_edges(
    const std::vector<Swap>& edges, std::size_t number_of_vertices) {
  GraphEdgeWeights result;
  for (const auto& edge : edges) {
    TKET_ASSERT(edge.first < number_of_vertices);
    TKET_ASSERT(edge.second < number_of_vertices);
    result[get_edge(edge.first, edge.second)] = 1;
  }
  TKET_ASSERT(edges.size() == result.size());
  return result;
}

static void add_solver_solutions(
    const std::vector<SolutionWSM>& solutions, std::size_t pattern_n_vertices,
    std::size_t target_n_vertices, SubgraphMonomorphisms& monomorphisms) {
  monomorphisms.mappings.reserve(solutions.size());

  // The solver ignores isolated vertices, so we'll fill them in separately.
  std::set<std::size_t> used_pattern_vertices;
  std::set<std::size_t> used_pattern_vertices_copy;
  std::set<std::size_t> used_target_vertices;

  for (const auto& solution : solutions) {
    monomorphisms.mappings.emplace_back();
    monomorphisms.mappings.back().resize(pattern_n_vertices);
    used_pattern_vertices.clear();
    used_target_vertices.clear();
    for (const auto& pair : solution.assignments) {
      const auto& pv = pair.first;
      const auto& tv = pair.second;
      TKET_ASSERT(used_pattern_vertices.insert(pv).second);
      TKET_ASSERT(used_target_vertices.insert(tv).second);
      if (!used_pattern_vertices_copy.empty()) {
        // It's always the same nonisolated pattern vertices.
        TKET_ASSERT(used_pattern_vertices_copy.count(pv) != 0);
      }
      // The architecture mapping guarantees
      // contiguous vertex numbers {0,1,2,...,N}
      TKET_ASSERT(pv < pattern_n_vertices);
      TKET_ASSERT(tv < target_n_vertices);
      monomorphisms.mappings.back()[pv] = tv;
    }
    if (used_pattern_vertices_copy.empty()) {
      used_pattern_vertices_copy = used_pattern_vertices;
    } else {
      TKET_ASSERT(
          used_pattern_vertices_copy.size() == used_pattern_vertices.size());
    }

    // Now, fill in isolated p-vertex values.
    std::size_t next_tv = 0;
    for (std::size_t pv = 0; pv < pattern_n_vertices; ++pv) {
      // We could save time by saving the isolated PVs,
      // but surely not worth it
      if (used_pattern_vertices.count(pv) != 0) {
        continue;
      }
      while (used_target_vertices.count(next_tv) != 0) {
        ++next_tv;
      }
      monomorphisms.mappings.back()[pv] = next_tv;
      ++next_tv;
    }
  }
}

// Strictly speaking, we should go through all (T choose k) subsets
// of target vertices, and all permutations of k pattern vertices,
// but don't bother for now - just give one solution.
static void fill_with_all_isolated_pattern_vertices(
    SubgraphMonomorphisms& monomorphisms, std::size_t pattern_n_vertices) {
  monomorphisms.mappings.resize(1);
  monomorphisms.mappings[0].resize(pattern_n_vertices);
  for (std::size_t ii = 0; ii < pattern_n_vertices; ++ii) {
    monomorphisms.mappings[0][ii] = ii;
  }
}

SubgraphMonomorphisms::SubgraphMonomorphisms(
    const ArchitectureMapping& pattern_arch_mapping,
    const ArchitectureMapping& target_arch_mapping,
    const Parameters& parameters)
    : time_taken_ms(0) {
  if (parameters.max_number_of_mappings == 0) {
    return;
  }
  const auto pattern_n_vertices = pattern_arch_mapping.number_of_vertices();
  const auto target_n_vertices = target_arch_mapping.number_of_vertices();
  if (pattern_n_vertices > target_n_vertices) {
    return;
  }
  const auto pattern_edges = pattern_arch_mapping.get_edges();
  const auto target_edges = target_arch_mapping.get_edges();
  if (pattern_edges.size() > target_edges.size()) {
    return;
  }
  if (pattern_edges.empty()) {
    // A pointless special case: all pattern vertices are isolated!
    fill_with_all_isolated_pattern_vertices(*this, pattern_n_vertices);
    return;
  }
  const auto pattern_edges_and_weights =
      get_weight_one_edges(pattern_edges, pattern_n_vertices);
  const auto target_edges_and_weights =
      get_weight_one_edges(target_edges, target_n_vertices);

  MainSolverParameters solver_parameters;
  solver_parameters.terminate_with_first_full_solution = false;
  solver_parameters.for_multiple_full_solutions_the_max_number_to_obtain =
      parameters.max_number_of_mappings;
  solver_parameters.timeout_ms = parameters.timeout_ms;

  const MainSolver main_solver(
      pattern_edges_and_weights, target_edges_and_weights, solver_parameters);

  const auto& solution_data = main_solver.get_solution_data();

  time_taken_ms =
      solution_data.initialisation_time_ms + solution_data.search_time_ms;

  add_solver_solutions(
      solution_data.solutions, pattern_n_vertices, target_n_vertices, *this);
}

}  // namespace tket
