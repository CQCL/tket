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

#include "WeightSubgrMono/Searching/WeightUpdater.hpp"

#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchNodeWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static bool add_edge_weights(
    const FixedData& fixed_data, VertexWSM pv, VertexWSM tv, VertexWSM other_pv,
    WeightWSM p_edge_weight, const Assignments& assignments,
    std::set<std::pair<VertexWSM, VertexWSM>>& p_edges_processed,
    SearchNodeWrapper& node_wrapper, WeightWSM max_weight) {
  const auto other_tv_citer = assignments.find(other_pv);
  if (other_tv_citer == assignments.cend()) {
    // The neighbouring p-vertex is not yet assigned.
    return true;
  }
  const auto p_edge = get_edge(pv, other_pv);
  if (p_edges_processed.count(p_edge) != 0) {
    // This edge was already processed.
    return true;
  }
  const VertexWSM& other_tv = other_tv_citer->second;
  const auto target_weight_opt =
      fixed_data.target_neighbours_data.get_edge_weight_opt(tv, other_tv);

  if (!target_weight_opt) {
    return false;
  }
  p_edges_processed.insert(p_edge);

  // We have a p-edge AND a t-edge.
  node_wrapper.add_p_edge_weights(p_edge_weight)
      .add_scalar_product(p_edge_weight * target_weight_opt.value());

  if (node_wrapper.get().current_scalar_product > max_weight) {
    return false;
  }
  return true;
}

bool WeightUpdater::operator()(
    const FixedData& fixed_data, const Assignments& assignments,
    SearchNodeWrapper& node_wrapper,
    std::size_t number_of_assignments_previously_processed_in_this_node,
    WeightWSM max_weight) const {
  const auto& chosen_assignments = node_wrapper.get().chosen_assignments;
  m_p_edges_processed.clear();
  for (auto ii = number_of_assignments_previously_processed_in_this_node;
       ii < chosen_assignments.size(); ++ii) {
    const auto& p_vertex = chosen_assignments[ii].first;
    const auto& t_vertex = chosen_assignments[ii].second;

    const auto& edges_and_weights =
        fixed_data.pattern_neighbours_data.get_neighbours_and_weights(p_vertex);

    for (const auto& entry : edges_and_weights) {
      if (!add_edge_weights(
              fixed_data, p_vertex, t_vertex, entry.first, entry.second,
              assignments, m_p_edges_processed, node_wrapper, max_weight)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
