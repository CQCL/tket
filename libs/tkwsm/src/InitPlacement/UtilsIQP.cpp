// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/InitPlacement/UtilsIQP.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/GraphTheoretic/VertexRelabelling.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

WeightWSM get_scalar_product_with_complete_target(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    WeightWSM implicit_target_weight,
    const std::vector<unsigned>& assignments) {
  WeightWSM scalar_product = 0;
  TKET_ASSERT(
      assignments.size() == pattern_ndata.get_number_of_nonisolated_vertices());
  for (VertexWSM pv = 0; pv < assignments.size(); ++pv) {
    const unsigned& tv = assignments[pv];
    TKET_ASSERT(tv < target_ndata.get_number_of_nonisolated_vertices());

    const auto& neighbours_and_weights =
        pattern_ndata.get_neighbours_and_weights(pv);
    // Only use edges (v1,v2) with v1<v2.
    for (auto citer = std::lower_bound(
             neighbours_and_weights.cbegin(), neighbours_and_weights.cend(),
             std::make_pair(pv, WeightWSM(0)));
         citer != neighbours_and_weights.cend(); ++citer) {
      const VertexWSM& other_pv = citer->first;
      const unsigned& other_tv = assignments.at(other_pv);
      TKET_ASSERT(other_tv < target_ndata.get_number_of_nonisolated_vertices());

      const WeightWSM& p_edge_weight = citer->second;
      WeightWSM t_edge_weight = implicit_target_weight;
      const auto explicit_t_edge_weight_opt =
          target_ndata.get_edge_weight_opt(tv, other_tv);
      if (explicit_t_edge_weight_opt) {
        t_edge_weight = explicit_t_edge_weight_opt.value();
      }
      scalar_product += p_edge_weight * t_edge_weight;
    }
  }
  return scalar_product;
}

GraphEdgeWeights get_relabelled_graph_data(
    const GraphEdgeWeights& graph_data, const VertexRelabelling& relabelling) {
  GraphEdgeWeights new_edges_and_weights;
  for (const auto& entry : graph_data) {
    const VertexWSM& old_v1 = entry.first.first;
    const VertexWSM& old_v2 = entry.first.second;

    new_edges_and_weights[get_edge(
        relabelling.get_new_label(old_v1), relabelling.get_new_label(old_v2))] =
        entry.second;
  }
  return new_edges_and_weights;
}

WeightWSM get_scalar_product_upper_bound_for_complete_target_graph(
    const NeighboursData& pattern_ndata,
    const NeighboursData& explicit_target_ndata,
    WeightWSM implicit_target_weight, WeightWSM extra_safety_factor) {
  // A very crude faster check, first.
  const std::size_t number_of_tv =
      explicit_target_ndata.get_number_of_nonisolated_vertices();
  const std::size_t total_number_of_target_edges =
      get_product_or_throw(number_of_tv, std::size_t(number_of_tv - 1)) / 2;

  if (pattern_ndata.get_number_of_nonisolated_vertices() > number_of_tv) {
    throw std::runtime_error("not enough target vertices");
  }
  TKET_ASSERT(
      total_number_of_target_edges >= pattern_ndata.get_number_of_edges());

  std::vector<WeightWSM> pattern_weights =
      pattern_ndata.get_weights_expensive();
  std::sort(pattern_weights.begin(), pattern_weights.end());

  std::vector<WeightWSM> explicit_target_weights =
      explicit_target_ndata.get_weights_expensive();
  std::sort(explicit_target_weights.begin(), explicit_target_weights.end());

  TKET_ASSERT(pattern_weights.size() <= total_number_of_target_edges);
  TKET_ASSERT(explicit_target_weights.size() <= total_number_of_target_edges);
  TKET_ASSERT(!pattern_weights.empty());
  TKET_ASSERT(!explicit_target_weights.empty());
  const WeightWSM highest_t_weight =
      std::max(explicit_target_weights.back(), implicit_target_weight);

  // Do cheaper checks first.
  const auto crude_single_product_upper_bound_opt =
      get_checked_product(highest_t_weight, pattern_weights.back());
  if (crude_single_product_upper_bound_opt) {
    auto crude_scalar_product_upper_bound_opt = get_checked_product(
        crude_single_product_upper_bound_opt.value(),
        WeightWSM(pattern_weights.size()));

    if (crude_scalar_product_upper_bound_opt) {
      // Demand a bit more, for good measure
      // (conversion to signed int types, etc.)
      if (get_checked_product(
              crude_scalar_product_upper_bound_opt.value(),
              extra_safety_factor)) {
        return crude_single_product_upper_bound_opt.value();
      }
    }
  }

  // If we reached here, the crudest check failed.
  // So, do a more expensive, but thorough check.
  // To maximise  sum a(i).b(p(i)) over all permutations p,
  // where (a(i)) is increasing, p must be such that (b(p(i))) is INCREASING.
  // Let's be very crude and treat the weight vectors as stacks;
  // pop the largest weights off the back.
  WeightWSM total_scalar_product = 0;
  std::size_t remaining_implicit_target_edges =
      total_number_of_target_edges - explicit_target_weights.size();
  while (!pattern_weights.empty()) {
    bool use_implicit_weight = remaining_implicit_target_edges > 0;
    // Check if there are any larger explicit target weights.
    // (Although, usually the implicit weight is larger).
    if (use_implicit_weight && !explicit_target_weights.empty() &&
        implicit_target_weight < explicit_target_weights.back()) {
      use_implicit_weight = false;
    }
    WeightWSM largest_t_weight;
    if (use_implicit_weight) {
      TKET_ASSERT(remaining_implicit_target_edges > 0);
      --remaining_implicit_target_edges;
      largest_t_weight = implicit_target_weight;
    } else {
      TKET_ASSERT(!explicit_target_weights.empty());
      largest_t_weight = explicit_target_weights.back();
      explicit_target_weights.pop_back();
    }
    const WeightWSM single_product =
        get_product_or_throw(pattern_weights.back(), largest_t_weight);
    pattern_weights.pop_back();
    total_scalar_product =
        get_sum_or_throw(total_scalar_product, single_product);
  }
  get_product_or_throw(total_scalar_product, extra_safety_factor);
  return total_scalar_product;
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
