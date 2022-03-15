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

#include "WeightSubgrMono/Searching/SharedData.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchBranch.hpp"
#include "WeightSubgrMono/Searching/ValueOrdering.hpp"
#include "WeightSubgrMono/Searching/VariableOrdering.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SharedData::SharedData(const FixedData& fd) : fixed_data(fd) {}

ReductionResult SharedData::initialise(SearchBranch& branch) {
  const auto result = branch.initialise(*this);
  if (result == ReductionResult::FINISHED) {
    solution_storage.add_full_solution(branch);
  }
  return result;
}


ReductionResult SharedData::search(
    SearchBranch& branch, VariableOrdering& var_ordering,
    ValueOrdering& val_ordering) {

  if (branch.move_down_has_been_called()) {
    // We have to move up.
    std::size_t level;
    branch.get_data(level);
    if(level == 0) {
      return ReductionResult::FINISHED;
    }
    const auto& this_node_new_assignments = branch.get_current_node_wrapper().get().chosen_assignments;
    TKET_ASSERT(!this_node_new_assignments.empty());

    // NOTE: when we move DOWN, we erase EXACTLY one possible PV->TV choice.
    // Therefore, when moving back up to here, Dom(PV) is the ONLY domain
    // with a potential problem; we need not check any others.
    const auto choice = this_node_new_assignments[0];
    if (!branch.backtrack()) {
      return ReductionResult::FINISHED;
    }
    // We've now moved up; but what happened to the domain of our
    // first chosen PV when we moved down from here last?
    auto& domain_for_chosen_pv = branch.get_current_node_wrapper().get()
        .pattern_v_to_possible_target_v.at(choice.first);
    if(domain_for_chosen_pv.empty()) {
      // Each node represents possibilities; we've exhausted all choices
      // from this node, BUT the nodes further up might still have
      // valid possibilities. So move up again! Recurse!
      return search(branch, var_ordering, val_ordering);
    }
    if(domain_for_chosen_pv.size() == 1) {
      // Only one possible choice, BUT treat it as a free choice.
      // This is needed to make nogoods correct,
      // when constructed only from the first chosen_assignments entries.
      const auto new_tv = *domain_for_chosen_pv.cbegin();
      TKET_ASSERT(choice.second != new_tv);
      branch.move_down(choice.first, new_tv);
    }
  }

  // Now, we've backtracked; we can start to move down as far as possible.
  WeightWSM max_weight;
  {
    const auto max_weight_opt =
        solution_storage.get_acceptable_scalar_product();
    if (max_weight_opt) {
      max_weight = max_weight_opt.value();
      if (max_weight == 0) {
        return ReductionResult::FINISHED;
      }
    } else {
      set_maximum(max_weight);
    }
  }

  for (;;) {
    const auto reduction_result = branch.reduce_current_node(*this, max_weight);

    if (reduction_result == ReductionResult::FAILURE) {
      solution_storage.add_partial_solution(branch);
      return reduction_result;
    }
    if (reduction_result == ReductionResult::SUCCESS) {
      // We can move down. So choose a variable and value...
      const auto& node = branch.get_current_node_wrapper().get();
      const auto next_pv = var_ordering.choose_next_variable(
          node, branch.get_assignments(), *this);

      const auto next_tv = val_ordering.get_target_value(
          node.pattern_v_to_possible_target_v.at(next_pv), *this);

      branch.move_down(next_pv, next_tv);
      continue;
    }
    TKET_ASSERT(reduction_result == ReductionResult::FINISHED);
    break;
  }
  solution_storage.add_full_solution(branch);
  // Reducing the node completely doesn't mean the SEARCH is finished,
  // it only means that moving down this particular branch has
  // found a full solution (and there may be more).
  return ReductionResult::SUCCESS;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
