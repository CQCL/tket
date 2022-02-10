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

#include "WeightSubgrMono/Searching/FixedData.hpp"

#include <algorithm>
#include <stdexcept>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool FixedData::initialise(
    const GraphEdgeWeights& p_data, const GraphEdgeWeights& t_data,
    DomainInitialiser::Parameters parameters) {
  initial_node.current_scalar_product = 0;
  initial_node.chosen_assignments.clear();
  initial_node.total_p_edge_weights = 0;
  total_p_edge_weights = 0;

  if (p_data.size() > t_data.size()) {
    // Who knows? Who cares? Doesn't actually matter!
    target_is_complete = false;
    return false;
  }
  pattern_neighbours_data.initialise(p_data);
  target_neighbours_data.initialise(t_data);

  const auto p_vertices =
      pattern_neighbours_data.get_nonisolated_vertices_expensive();
  const auto t_vertices =
      target_neighbours_data.get_nonisolated_vertices_expensive();

  const auto number_of_possible_t_edges =
      (t_vertices.size() * (t_vertices.size() - 1)) / 2;
  if (!(t_data.size() <= number_of_possible_t_edges)) {
    throw std::runtime_error("Invalid target graph input data");
  }
  target_is_complete = t_data.size() == number_of_possible_t_edges;

  if (p_vertices.size() > t_vertices.size()) {
    return false;
  }
  if (p_vertices.size() <= 1) {
    return true;
  }

  auto& node_domains_map = initial_node.pattern_v_to_possible_target_v;
  node_domains_map.clear();

  if (target_is_complete) {
    const std::set<VertexWSM> full_domain{
        t_vertices.cbegin(), t_vertices.cend()};
    for (auto pv : p_vertices) {
      node_domains_map[pv] = full_domain;
    }
    return true;
  }

  DomainInitialiser initialiser;
  const bool success = initialiser.full_initialisation(
      initial_node.pattern_v_to_possible_target_v, p_vertices,
      pattern_neighbours_data, t_vertices, target_neighbours_data, parameters);

  if (!success) {
    return false;
  }
  const auto& assigned_pv_list = initialiser.get_assigned_vertices();
  initial_node.chosen_assignments.clear();
  for (auto pv : assigned_pv_list) {
    const auto iter = initial_node.pattern_v_to_possible_target_v.find(pv);
    TKET_ASSERT(iter != initial_node.pattern_v_to_possible_target_v.end());
    TKET_ASSERT(iter->second.size() == 1);

    initial_node.chosen_assignments.emplace_back(pv, *iter->second.cbegin());

    initial_node.pattern_v_to_possible_target_v.erase(iter);
  }
  for (const auto& entry : p_data) {
    total_p_edge_weights += entry.second;
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
