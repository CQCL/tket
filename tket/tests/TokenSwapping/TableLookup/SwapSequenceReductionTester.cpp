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

#include "SwapSequenceReductionTester.hpp"

#include <catch2/catch_test_macros.hpp>

#include "NeighboursFromEdges.hpp"
#include "TokenSwapping/SwapListSegmentOptimiser.hpp"
#include "TokenSwapping/VertexMapResizing.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"
#include "TokenSwapping/VertexSwapResult.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

static void reduce_sequence(
    const vector<Swap>& swaps, const VertexMapping& vertex_mapping,
    NeighboursFromEdges& neighbours, SwapList& raw_swap_list,
    SwapListOptimiser& general_optimiser,
    const SwapSequenceReductionTester::Options& options) {
  REQUIRE(!swaps.empty());

  VertexMapResizing map_resizing(neighbours);
  SwapListTableOptimiser table_optimiser;
  SwapListSegmentOptimiser& segment_optimiser =
      table_optimiser.get_segment_optimiser();
  raw_swap_list.clear();
  for (const auto& swap : swaps) {
    raw_swap_list.push_back(swap);
  }
  std::set<size_t> vertices_with_tokens;
  for (const auto& entry : vertex_mapping) {
    vertices_with_tokens.insert(entry.first);
  }

  if (options.optimise_initial_segment_only) {
    general_optimiser.optimise_pass_with_frontward_travel(raw_swap_list);
    if (!raw_swap_list.empty()) {
      table_optimiser.get_segment_optimiser().optimise_segment(
          raw_swap_list.front_id().value(), vertices_with_tokens, map_resizing,
          raw_swap_list);
    }
    return;
  }
  table_optimiser.optimise(
      vertices_with_tokens, map_resizing, raw_swap_list, general_optimiser);
}

static void check_solution(
    VertexMapping problem_vertex_mapping, const SwapList& raw_swap_list) {
  // Every vertex swap on a source->target mapping converts it to a new
  // source->target map, i.e. map[v] = (token currently at v).
  // So we BEGIN with every token equalling its target,
  // thus at the end every token must equal its vertex.
  for (auto id_opt = raw_swap_list.front_id(); id_opt;) {
    const auto id = id_opt.value();
    id_opt = raw_swap_list.next(id);
    const auto& swap = raw_swap_list.at(id);
    const VertexSwapResult vswap_result(swap, problem_vertex_mapping);
  }
  REQUIRE(all_tokens_home(problem_vertex_mapping));
}

static size_t get_reduced_swaps_size_with_checks(
    const vector<Swap>& swaps, const VertexMapping& problem_vertex_mapping,
    NeighboursFromEdges& neighbours_calculator,
    SwapListOptimiser& general_optimiser,
    const SwapSequenceReductionTester::Options& options) {
  SwapList raw_swap_list;
  reduce_sequence(
      swaps, problem_vertex_mapping, neighbours_calculator, raw_swap_list,
      general_optimiser, options);
  check_solution(problem_vertex_mapping, raw_swap_list);
  REQUIRE(raw_swap_list.size() <= swaps.size());
  return raw_swap_list.size();
}

size_t SwapSequenceReductionTester::get_checked_solution_size(
    const DecodedProblemData& problem_data,
    const SwapSequenceReductionTester::Options& options) {
  NeighboursFromEdges neighbours_calculator(problem_data.swaps);
  return get_reduced_swaps_size_with_checks(
      problem_data.swaps, problem_data.vertex_mapping, neighbours_calculator,
      m_general_optimiser, options);
}

// Reduces the sequence of swaps, checks it, and returns the size.
size_t SwapSequenceReductionTester::get_checked_solution_size(
    const DecodedProblemData& problem_data,
    const DecodedArchitectureData& architecture_data,
    const SwapSequenceReductionTester::Options& options) {
  NeighboursFromEdges neighbours_calculator(architecture_data.edges);
  return get_reduced_swaps_size_with_checks(
      problem_data.swaps, problem_data.vertex_mapping, neighbours_calculator,
      m_general_optimiser, options);
}

SequenceReductionStats::SequenceReductionStats()
    : problems(0),
      reduced_problems(0),
      total_original_swaps(0),
      total_original_swaps_for_reduced_problems(0),
      total_reduced_swaps(0) {}

void SequenceReductionStats::add_solution(
    size_t original_swaps, size_t reduced_swaps) {
  REQUIRE(reduced_swaps <= original_swaps);
  ++problems;
  if (reduced_swaps < original_swaps) {
    ++reduced_problems;
    total_original_swaps_for_reduced_problems += original_swaps;
  }
  total_reduced_swaps += reduced_swaps;
  total_original_swaps += original_swaps;
}

std::string SequenceReductionStats::str() const {
  std::stringstream ss;
  const size_t swaps_for_equal_probs =
      total_original_swaps - total_original_swaps_for_reduced_problems;
  const size_t reduced_swaps_for_reduced_probs =
      total_reduced_swaps - swaps_for_equal_probs;
  const size_t overall_decrease = total_original_swaps - total_reduced_swaps;
  ss << "[" << problems - reduced_problems << " equal probs ("
     << swaps_for_equal_probs << "); " << reduced_problems << " reduced probs ("
     << reduced_swaps_for_reduced_probs << " vs "
     << total_original_swaps_for_reduced_problems << ")]\n[Overall reduction "
     << total_reduced_swaps << " vs " << total_original_swaps << ": ";
  if (total_original_swaps == 0) {
    ss << "0%";
  } else {
    ss << (100 * overall_decrease) / total_original_swaps << "%";
  }
  ss << "]";
  return ss.str();
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
