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

#include "WeightSubgrMono/EndToEndWrappers/CheckedWeightBounds.hpp"

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SpecialExceptions.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// KEY: a target vertex TV.
// VALUE: the max weight of any edge containing TV.
static std::map<VertexWSM, WeightWSM> get_max_t_edge_weight_map(
    const GraphEdgeWeights& target_edges) {
  std::map<VertexWSM, WeightWSM> max_t_edge_weight_map;
  for (const auto& entry : target_edges) {
    // If a new key appears in max_t_edge_weight_map,
    // the new value is zero.
    // Works because we want the maximum.
    {
      auto& weight = max_t_edge_weight_map[entry.first.first];
      weight = std::max(weight, entry.second);
    }
    {
      auto& weight = max_t_edge_weight_map[entry.first.second];
      weight = std::max(weight, entry.second);
    }
  }
  return max_t_edge_weight_map;
}

static WeightWSM get_maximum_weight(const GraphEdgeWeights& data) {
  WeightWSM max_weight = 0;
  for (const auto& entry : data) {
    max_weight = std::max(max_weight, entry.second);
  }
  return max_weight;
}

CheckedWeightBounds::CheckedWeightBounds(
    const FixedData& fixed_data, const GraphEdgeWeights& pattern_edges,
    const GraphEdgeWeights& target_edges, WeightWSM extra_safety_factor) {
  other_inconsistency_occurred = false;

  // First, most crude: largest p-weight * largest t-weight * number.
  {
    const auto weight_product_opt = get_checked_product(
        get_maximum_weight(pattern_edges), get_maximum_weight(target_edges));

    if (weight_product_opt) {
      const auto total_product_opt = get_checked_product(
          weight_product_opt.value(), WeightWSM(pattern_edges.size()));
      if (total_product_opt) {
        // The crude upper bound hasn't overflowed;
        // but remember the extra factor.
        const auto upper_bound_to_use = total_product_opt.value();
        const auto final_product_opt =
            get_checked_product(upper_bound_to_use, extra_safety_factor);
        if (final_product_opt) {
          // The crudest lower bound imaginable!
          lower_bound = 0;
          upper_bound = upper_bound_to_use;
          return;
        }
      }
    }
  }

  // The simplest check overflowed; try a bit more refined.
  // For each p-edge weight, multiply ONLY by the max t-edge weight
  // with a t-vertex in the domain.

  // For a TV, the maximum t-edge weight containing TV.
  const auto max_t_edge_weight_map = get_max_t_edge_weight_map(target_edges);

  // KEY: a PV
  // VALUE: the maximum t-edge weight that any p-edge containing PV
  //      could be assigned to.
  std::map<VertexWSM, WeightWSM> maximum_t_edge_weights_from_pattern_v;
  for (const auto& entry : fixed_data.initial_node.chosen_assignments) {
    maximum_t_edge_weights_from_pattern_v[entry.first] =
        max_t_edge_weight_map.at(entry.second);
  }
  // Now, fill values for the domains.
  for (const auto& entry :
       fixed_data.initial_node.pattern_v_to_possible_target_v) {
    const auto& pv = entry.first;
    const auto& domain = entry.second;
    // The value is automatically set to zero on first use.
    auto& weight = maximum_t_edge_weights_from_pattern_v[pv];
    for (auto tv : domain) {
      weight = std::max(weight, max_t_edge_weight_map.at(tv));
    }
  }

  // Finally, go through the p-edges one-by-one.
  WeightWSM current_upper_bound = 0;

  for (const auto& entry : pattern_edges) {
    const auto& pv1 = entry.first.first;
    auto t_weight_estimate = maximum_t_edge_weights_from_pattern_v.at(pv1);

    const auto& pv2 = entry.first.second;
    t_weight_estimate = std::max(
        t_weight_estimate, maximum_t_edge_weights_from_pattern_v.at(pv2));
    const auto scalar_product =
        get_product_or_throw(entry.second, t_weight_estimate);
    current_upper_bound = get_sum_or_throw(current_upper_bound, scalar_product);
  }
  get_product_or_throw(current_upper_bound, extra_safety_factor);
  // We've got an estimate, which doesn't overflow.
  lower_bound = 0;
  upper_bound = current_upper_bound;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
