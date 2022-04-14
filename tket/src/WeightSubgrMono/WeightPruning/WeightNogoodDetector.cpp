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

#include "WeightSubgrMono/WeightPruning/WeightNogoodDetector.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

WeightNogoodDetector::WeightNogoodDetector(
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data,
    const std::set<VertexWSM>& initial_used_target_vertices)
    : m_pattern_neighbours_data(pattern_neighbours_data),
      m_target_neighbours_data(target_neighbours_data),
      m_valid_target_vertices(
          initial_used_target_vertices.cbegin(),
          initial_used_target_vertices.cend()) {}

std::optional<WeightWSM> WeightNogoodDetector::get_min_weight_for_tv(
    VertexWSM tv) const {
  if (m_valid_target_vertices.count(tv) == 0) {
    return {};
  }
  const auto weight_opt = get_optional_value(m_minimum_t_weights_from_tv, tv);
  if (weight_opt) {
    return weight_opt;
  }
  // We must find the minimum weight, by looking at all neighbours.
  WeightWSM min_weight;
  set_maximum(min_weight);
  const auto& data = m_target_neighbours_data.get_neighbours_and_weights(tv);
  for (const auto& entry : data) {
    const auto& neighbour_tv = entry.first;
    if (m_valid_target_vertices.count(neighbour_tv) == 0) {
      continue;
    }
    min_weight = std::min(min_weight, entry.second);
  }
  // Really, to do this properly, we should update
  // target_neighbours_data dynamically,
  // erasing invalid target vertices and updating neighbour lists.
  if (is_maximum(min_weight)) {
    m_valid_target_vertices.erase(tv);
    return {};
  }
  m_minimum_t_weights_from_tv[tv] = min_weight;
  return min_weight;
}

WeightWSM WeightNogoodDetector::get_t_weight_lower_bound(VertexWSM pv) const {
  // Of course, std::lower_bound is unfortunate terminology;
  // completely unrelated to our weight lower bounds!
  const auto citer = std::lower_bound(
      m_t_weight_lower_bounds_for_p_edges_containing_pv.cbegin(),
      m_t_weight_lower_bounds_for_p_edges_containing_pv.cend(),
      std::make_pair(pv, WeightWSM(0)));
  TKET_ASSERT(
      citer != m_t_weight_lower_bounds_for_p_edges_containing_pv.cend());
  TKET_ASSERT(citer->first == pv);
  return citer->second;
}

bool WeightNogoodDetector::fill_t_weight_lower_bounds_for_p_edges_containing_pv(
    const PossibleAssignments& possible_assignments) const {
  // Every unassigned edge must connect to an unassigned p-vertex.
  m_t_weight_lower_bounds_for_p_edges_containing_pv.clear();

  for (const auto& entry : possible_assignments) {
    WeightWSM weight;
    set_maximum(weight);
    const auto& domain = entry.second;
    TKET_ASSERT(!domain.empty());
    for (auto tv : domain) {
      const auto weight_opt_for_tv = get_min_weight_for_tv(tv);
      if (weight_opt_for_tv) {
        weight = std::min(weight, weight_opt_for_tv.value());
      }
    }
    if (is_maximum(weight)) {
      // A nogood found already!
      return false;
    }
    m_t_weight_lower_bounds_for_p_edges_containing_pv.emplace_back(
        entry.first, weight);
  }
  return true;
}

WeightNogoodDetector::Result WeightNogoodDetector::operator()(
    const PossibleAssignments& possible_assignments,
    WeightWSM max_extra_scalar_product) const {
  Result result;
  if (!fill_t_weight_lower_bounds_for_p_edges_containing_pv(
          possible_assignments)) {
    // A nogood!
    return result;
  }
  WeightWSM weight_lower_bound = 0;

  // Now, for each unassigned p-vertex, look at its neighbours
  // and deduce the p-edge weights.
  // There's a problem: if BOTH p-edge endpoints are unassigned,
  // the edge will be counted twice.
  // To solve this: only add the data when pv1 < pv2.
  for (const auto& entry : possible_assignments) {
    if (entry.second.size() == 1) {
      // It's assigned.
      continue;
    }
    const auto& pv1 = entry.first;

    // Note: we don't care about the domain of pv1,
    // we've already gone through it before (when we called
    // fill_t_weight_lower_bounds_for_p_edges_containing_pv).

    // A lower bound for all t-edges
    // which could possibly contain f(pv).
    const auto minimum_t_weight = get_t_weight_lower_bound(pv1);

    const auto& p_neighbours_and_weights =
        m_pattern_neighbours_data.get_neighbours_and_weights(pv1);

    for (const auto& pv2_weight_pair : p_neighbours_and_weights) {
      const VertexWSM& pv2 = pv2_weight_pair.first;
      const WeightWSM& p_weight = pv2_weight_pair.second;
      WeightWSM t_weight_estimate = minimum_t_weight;
      const auto& domain2 = possible_assignments.at(pv2);
      if (domain2.size() == 1) {
        // This other p-vertex PV2 is assigned already.
        // So let's check the other t-weight estimate, it may be better.
        const auto& tv2 = *domain2.cbegin();

        // We already know that this edge contains pv1,
        // so DEFINITELY has t-weight >= this current estimate.
        // But we also know that it will be assigned to a target edge
        // containing tv2, which has t_weight >= the other value.
        const auto other_tv_weight_bound_opt = get_min_weight_for_tv(tv2);

        if (!other_tv_weight_bound_opt) {
          // We're at a nogood!
          // (This COULD actually happen. It means that TV2 is invalid,
          // and in fact always was;
          // but we didn't realise this at the time we made the assignment).
          // Definitely worth the caller trying to make use of this new
          // information, although it is algorithmically complicated.
          result.invalid_t_vertex = tv2;
          return result;
        }

        // We take the MAX to get a valid LOWER bound
        // as LARGE as possible.
        t_weight_estimate =
            std::max(t_weight_estimate, other_tv_weight_bound_opt.value());
      } else {
        // pv2 is ALSO unassigned. Beware of double counting!
        if (pv1 > pv2) {
          continue;
        }
        // We'll do BOTH pv1--pv2 and pv2--pv1 NOW.
        // If pv2 produces a larger value, we can use that instead.
        t_weight_estimate =
            std::max(t_weight_estimate, get_t_weight_lower_bound(pv2));
      }
      weight_lower_bound += p_weight * t_weight_estimate;
      if (weight_lower_bound > max_extra_scalar_product) {
        return result;
      }
    }
  }
  result.extra_scalar_product_lower_bound = weight_lower_bound;
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
