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

#include "WeightSubgrMono/Searching/CompleteTargetManager.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SearchNode.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static void fill_pattern_edge_sums(
    std::map<VertexWSM, WeightWSM>& edge_sums, const NeighboursData& ndata) {
  // We only care about the KEYS, we do NOT rely on the internal implementation.
  const auto& map = ndata.get_map();
  for (const auto& entry : map) {
    const auto& vertex = entry.first;

    // Automatically set to 0.
    auto& sum = edge_sums[vertex];
    // const auto& neighbours_and_weights =
    // ndata.get_neighbours_and_weights(vertex);
    for (const auto& pair : ndata.get_neighbours_and_weights(vertex)) {
      sum += pair.second;
    }
  }
}

static void fill_target_cumulative_edge_sums(
    std::map<VertexWSM, std::vector<WeightWSM>>& edge_sums,
    const NeighboursData& ndata) {
  const auto& map = ndata.get_map();
  TKET_ASSERT(map.size() > 0);
  std::vector<unsigned> indices;
  indices.resize(map.size() - 1);
  for (unsigned ii = 0; ii < indices.size(); ++ii) {
    indices[ii] = ii;
  }
  for (const auto& entry : map) {
    const auto& tv = entry.first;
    const auto& neighbours_and_weights = entry.second;
    TKET_ASSERT(neighbours_and_weights.size() == indices.size());

    // Sort by weights.
    // Nonstable sorts don't matter as we discard the vertex numbers
    // (i.e., the INDICES are platform-dependent, but the final results
    // are not!)
    std::sort(
        indices.begin(), indices.end(),
        [&neighbours_and_weights](unsigned ii, unsigned jj) {
          return neighbours_and_weights[ii].second <
                 neighbours_and_weights[jj].second;
        });

    WeightWSM cumulative_sum = 0;
    auto& sums_vector = edge_sums[tv];
    sums_vector.resize(indices.size());
    for (unsigned ii = 0; ii < indices.size(); ++ii) {
      cumulative_sum += neighbours_and_weights[indices[ii]].second;
      sums_vector[ii] = cumulative_sum;
    }
  }
}

CompleteTargetManager::CompleteTargetManager(const FixedData& fixed_data)
    : m_fixed_data(fixed_data) {
  TKET_ASSERT(fixed_data.target_is_complete);
  fill_pattern_edge_sums(
      m_pattern_edge_sums, fixed_data.pattern_neighbours_data);
  fill_target_cumulative_edge_sums(
      m_target_partial_edge_sums, fixed_data.target_neighbours_data);
}

VertexWSM CompleteTargetManager::choose_variable(
    const SearchNode& node, const Assignments& assignments) const {
  TKET_ASSERT(node.pattern_v_to_possible_target_v.size() > 0);
  VertexWSM pv_to_return = node.pattern_v_to_possible_target_v.cbegin()->first;
  WeightWSM highest_weight = 0;

  // Because it's a complete target graph, we don't worry about building up
  // the assignments by growing a connected subgraph.
  for (const auto& entry : node.pattern_v_to_possible_target_v) {
    const auto& new_pv = entry.first;
    TKET_ASSERT(assignments.count(new_pv) == 0);
    const auto& new_weight = m_pattern_edge_sums.at(new_pv);
    if (new_weight > highest_weight) {
      pv_to_return = new_pv;
      highest_weight = new_weight;
    }
  }
  return pv_to_return;
}

std::pair<VertexWSM, VertexWSM> CompleteTargetManager::choose_next_assignment(
    const SearchNode& node, const Assignments& assignments) const {
  const auto pv = choose_variable(node, assignments);
  const auto number_of_pv_neighbours =
      m_fixed_data.pattern_neighbours_data.get_neighbours_and_weights(pv)
          .size();
  TKET_ASSERT(number_of_pv_neighbours > 0);
  const auto& domain = node.pattern_v_to_possible_target_v.at(pv);
  TKET_ASSERT(domain.size() > 1);
  VertexWSM best_tv = *domain.cbegin();
  WeightWSM best_tv_sum;
  set_maximum(best_tv_sum);
  for (auto tv : domain) {
    const auto& tv_sum =
        m_target_partial_edge_sums.at(tv).at(number_of_pv_neighbours - 1);
    if (tv_sum < best_tv_sum) {
      best_tv = tv;
      best_tv_sum = tv_sum;
    }
  }
  return std::make_pair(pv, best_tv);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
