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

bool SharedData::search(
    SearchBranch& branch, VariableOrdering& var_ordering,
    ValueOrdering& val_ordering) {
  if (branch.move_down_has_been_called()) {
    // We have to move up.
    if (!branch.move_up()) {
      return false;
    }
    // We've now backtracked.
  }
  // Move down as far as possible.
  WeightWSM max_weight;
  {
    const auto max_weight_opt =
        solution_storage.get_acceptable_scalar_product();
    if (max_weight_opt) {
      max_weight = max_weight_opt.value();
      if (max_weight == 0) {
        return false;
      }
    } else {
      set_maximum(max_weight);
    }
  }

  for (;;) {
    const auto reduction_result = branch.reduce_current_node(*this, max_weight);

    if (reduction_result == ReductionResult::FAILURE) {
      // We're at a nogood.
      solution_storage.add_partial_solution(branch);
      return true;
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
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
