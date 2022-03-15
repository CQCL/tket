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

#include "WeightSubgrMono/Reducing/AllDiffPropagator.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool AllDiffPropagator::reduce(
    Assignments& assignments, SearchNodeWrapper& node_wrapper,
    std::size_t& number_of_assignments_previously_processed_in_this_node)
    const {
  auto& node = node_wrapper.get_mutable();
  auto& pv_to_domains_map = node.pattern_v_to_possible_target_v;

  while (number_of_assignments_previously_processed_in_this_node <
         node.chosen_assignments.size()) {
    const auto& new_assignment =
        node.chosen_assignments
            [number_of_assignments_previously_processed_in_this_node];
    ++number_of_assignments_previously_processed_in_this_node;
    const auto& p_vertex = new_assignment.first;
    const auto& t_vertex = new_assignment.second;
    const auto tv_opt = get_optional_value(assignments, p_vertex);
    if (tv_opt) {
      if (tv_opt.value() != t_vertex) {
        return false;
      }
    } else {
      assignments[p_vertex] = t_vertex;
    }
    {
      const auto domain_citer = pv_to_domains_map.find(p_vertex);
      if (domain_citer != pv_to_domains_map.cend()) {
        const auto& domain = domain_citer->second;
        if (!domain.empty()) {
          TKET_ASSERT(domain.size() == 1);
          TKET_ASSERT(*domain.cbegin() == t_vertex);
        }
      }
    }
    pv_to_domains_map.erase(p_vertex);

    for (const auto& entry : pv_to_domains_map) {
      const auto reduced_size = node_wrapper.remove_element_from_domain(
          entry.first, t_vertex, assignments);
      if (reduced_size == 0) {
        return false;
      }
    }
  }

  // Now, all domains should have size >= 2...
  TKET_ASSERT(
      number_of_assignments_previously_processed_in_this_node ==
      node.chosen_assignments.size());

  for (const auto& entry : pv_to_domains_map) {
    TKET_ASSERT(entry.second.size() > 1);
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
