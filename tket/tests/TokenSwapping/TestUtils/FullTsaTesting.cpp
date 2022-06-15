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

#include "FullTsaTesting.hpp"

#include <catch2/catch_test_macros.hpp>

#include "Architecture/ArchitectureMapping.hpp"
#include "Architecture/DistancesFromArchitecture.hpp"
#include "Architecture/NeighboursFromArchitecture.hpp"
#include "DebugFunctions.hpp"
#include "TokenSwapping/DistanceFunctions.hpp"
#include "TokenSwapping/RiverFlowPathFinder.hpp"
#include "TokenSwapping/VertexSwapResult.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

void FullTsaTesting::check_solution(
    size_t counts_list_index, VertexMapping vertex_mapping, size_t lower_bound,
    AllowEmptySwaps allow_empty_swaps) {
  bool empty_swap_occurred = false;
  REQUIRE(m_swap_list.size() >= lower_bound);
  for (auto swap : m_swap_list.to_vector()) {
    const VertexSwapResult swap_res(swap, vertex_mapping);
    if (swap_res.tokens_moved == 0) {
      empty_swap_occurred = true;
    }
  }
  if (empty_swap_occurred && allow_empty_swaps == AllowEmptySwaps::NO) {
    INFO(
        "index=" << counts_list_index << ", " << vertex_mapping.size()
                 << " toks; lb=" << lower_bound << "; " << m_swap_list.size()
                 << " swaps");
    CHECK(false);
  }
  REQUIRE(all_tokens_home(vertex_mapping));
  auto& swaps = m_counts_list[counts_list_index].sorted_swaps;
  swaps = m_swap_list.to_vector();
  std::sort(swaps.begin(), swaps.end());
}

void FullTsaTesting::check_equivalent_good_solution(
    size_t existing_index, VertexMapping vertex_mapping,
    AllowEmptySwaps allow_empty_swaps) {
  check_solution(
      m_counts_list.size() - 1, vertex_mapping, 0, allow_empty_swaps);
  INFO("existing_index=" << existing_index);
  CHECK(
      m_counts_list[existing_index].sorted_swaps ==
      m_counts_list.back().sorted_swaps);
}

void FullTsaTesting::test_order(size_t index1, size_t index2) const {
  INFO("i1=" << index1 << ", i2=" << index2);
  CHECK(
      m_counts_list[index1].sorted_swaps.size() <=
      m_counts_list[index2].sorted_swaps.size());
}

void FullTsaTesting::complete_counts_list_for_single_problem() {
  size_t smallest_number = m_counts_list[0].sorted_swaps.size();

  // Ignore the last index, which is a dummy.
  for (size_t index = 0; index + 1 < m_counts_list.size(); ++index) {
    auto& counts = m_counts_list[index];
    counts.total_swaps += counts.sorted_swaps.size();
    smallest_number = std::min(smallest_number, counts.sorted_swaps.size());
  }
  // Now, we've got the (joint) winner.
  size_t best_index = m_counts_list.size();
  size_t num_winners = 0;

  for (size_t index = 0; index + 1 < m_counts_list.size(); ++index) {
    auto& counts = m_counts_list[index];
    REQUIRE(counts.sorted_swaps.size() >= smallest_number);
    if (counts.sorted_swaps.size() == smallest_number) {
      ++counts.problems_where_this_was_the_joint_winner;
      ++num_winners;
      best_index = index;
    }
  }
  REQUIRE(num_winners >= 1);
  if (num_winners == 1) {
    ++m_counts_list[best_index].problems_where_this_was_the_clear_winner;
  }
}

FullTsaTesting::FullTsaTesting() {
  m_counts_list.resize(7);
  for (auto& entry : m_counts_list) {
    entry.total_swaps = 0;
  }
}

void FullTsaTesting::add_problems(
    const ArchitectureMapping& arch_mapping,
    const vector<VertexMapping>& problems, const std::string& new_name,
    RNG& rng, PartialTsaInterface& full_tsa) {
  m_number_of_problems += problems.size();
  const std::string name_for_this = new_name + ":" + full_tsa.name();
  if (m_name.empty()) {
    m_name = name_for_this;
  } else {
    if (m_name != name_for_this) {
      m_name = m_name + ":" + name_for_this;
    }
  }
  DistancesFromArchitecture distances(arch_mapping);
  NeighboursFromArchitecture neighbours(arch_mapping);
  RiverFlowPathFinder path_finder(distances, neighbours, rng);
  vector<Swap> raw_calc_swaps;
  VertexMapping problem_copy_to_destroy;

  for (size_t prob_index = 0; prob_index < problems.size(); ++prob_index) {
    const auto& problem = problems[prob_index];
    const auto lower_bound = get_swaps_lower_bound(problem, distances);
    m_number_of_tokens += problem.size();
    m_total_lower_bounds += lower_bound;
    problem_copy_to_destroy = problem;
    m_swap_list.clear();
    rng.set_seed();
    full_tsa.append_partial_solution(
        m_swap_list, problem_copy_to_destroy, distances, neighbours,
        path_finder);
    raw_calc_swaps = m_swap_list.to_vector();

    // Now, let's check the calculated swaps.
    check_solution(0, problem, lower_bound, AllowEmptySwaps::NO);

    // Minimal travel optimising
    m_optimiser.optimise_pass_with_zero_travel(m_swap_list);
    check_solution(1, problem, lower_bound, AllowEmptySwaps::NO);
    test_order(1, 0);

    //...add artificial token tracking...(remembering that empty swaps
    // can be introduced, since it knows nothing about our tokens).
    m_optimiser.optimise_pass_with_token_tracking(m_swap_list);
    check_solution(2, problem, lower_bound, AllowEmptySwaps::YES);
    test_order(2, 1);

    m_optimiser.optimise_pass_remove_empty_swaps(m_swap_list, problem);
    check_solution(3, problem, lower_bound, AllowEmptySwaps::NO);
    test_order(3, 2);

    m_optimiser.full_optimise(m_swap_list, problem);
    check_solution(4, problem, lower_bound, AllowEmptySwaps::NO);
    test_order(4, 3);

    // Now, test various equalities.

    // The token tracking pass, by itself, is the same whether or not
    // we zero travel optimise first (which just makes things faster,
    // not better).
    m_swap_list.clear();
    for (auto swap : raw_calc_swaps) {
      m_swap_list.push_back(swap);
    }
    m_optimiser.optimise_pass_with_token_tracking(m_swap_list);
    m_optimiser.optimise_pass_with_frontward_travel(m_swap_list);
    // Is 5 the same as 2? No! Usually the same, but NOT always;
    // e.g. a test with random trees found a small difference.
    check_solution(5, problem, lower_bound, AllowEmptySwaps::YES);

    // Swap travels permute the swaps, but otherwise reduce them
    // no more than zero travel.
    m_swap_list.clear();
    for (auto swap : raw_calc_swaps) {
      m_swap_list.push_back(swap);
    }
    m_optimiser.optimise_pass_with_frontward_travel(m_swap_list);
    check_equivalent_good_solution(1, problem, AllowEmptySwaps::NO);

    // full optimise is no better when combined
    // with other passes.
    m_swap_list.clear();
    for (auto swap : raw_calc_swaps) {
      m_swap_list.push_back(swap);
    }
    m_optimiser.full_optimise(m_swap_list);
    check_equivalent_good_solution(2, problem, AllowEmptySwaps::YES);
    m_optimiser.optimise_pass_with_token_tracking(m_swap_list);
    check_equivalent_good_solution(2, problem, AllowEmptySwaps::YES);

    m_swap_list.clear();
    for (auto swap : raw_calc_swaps) {
      m_swap_list.push_back(swap);
    }
    m_optimiser.full_optimise(m_swap_list, problem);
    check_equivalent_good_solution(4, problem, AllowEmptySwaps::NO);

    complete_counts_list_for_single_problem();
  }
}

std::string FullTsaTesting::str() const {
  std::stringstream ss;
  ss << "[" << m_name << ": " << m_number_of_problems << " probs; "
     << m_number_of_tokens << " toks; " << m_total_lower_bounds
     << " tot.lb]\n[Total swaps:";
  // The last entry is a "dummy".
  for (size_t index = 0; index + 1 < m_counts_list.size(); ++index) {
    ss << " " << m_counts_list[index].total_swaps;
  }
  ss << "]\n[Winners: joint:";
  for (size_t index = 0; index + 1 < m_counts_list.size(); ++index) {
    ss << " " << m_counts_list[index].problems_where_this_was_the_joint_winner;
  }
  ss << "  undisputed:";
  for (size_t index = 0; index + 1 < m_counts_list.size(); ++index) {
    ss << " " << m_counts_list[index].problems_where_this_was_the_clear_winner;
  }
  ss << "]";
  return ss.str();
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
