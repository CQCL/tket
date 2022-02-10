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

#include "WeightSubgrMono/Searching/SearchNode.hpp"

#include <cmath>
#include <sstream>

namespace tket {
namespace WeightedSubgraphMonomorphism {

double SearchNode::get_log10_search_space_size() const {
  double total = 0;
  for (const auto& entry : pattern_v_to_possible_target_v) {
    total += std::log10(entry.second.size());
  }
  return total;
}

void SearchNode::initialise_from_assignment(
    VertexWSM pattern_v, VertexWSM target_v, const SearchNode& previous_node) {
  current_scalar_product = previous_node.current_scalar_product;
  total_p_edge_weights = previous_node.total_p_edge_weights;

  chosen_assignments.resize(1);
  chosen_assignments[0].first = pattern_v;
  chosen_assignments[0].second = target_v;

  pattern_v_to_possible_target_v = previous_node.pattern_v_to_possible_target_v;
  pattern_v_to_possible_target_v.erase(pattern_v);
}

std::string SearchNode::str() const {
  std::stringstream ss;
  ss << "\n### Node: weight " << current_scalar_product
     << "; total p edge weights " << total_p_edge_weights << ". Made "
     << chosen_assignments.size() << " assignments: [";
  for (const auto& entry : chosen_assignments) {
    ss << entry.first << ":" << entry.second << "  ";
  }
  ss << "]\nStill " << pattern_v_to_possible_target_v.size()
     << " unassigned vars:";

  for (const auto& entry : pattern_v_to_possible_target_v) {
    ss << "\nDom(" << entry.first << ") = {";
    for (auto vv : entry.second) {
      ss << vv << " ";
    }
    ss << "}";
  }
  ss << "\n";
  return ss.str();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
