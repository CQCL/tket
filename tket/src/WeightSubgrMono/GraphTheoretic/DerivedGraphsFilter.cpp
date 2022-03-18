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

#include "WeightSubgrMono/GraphTheoretic/DerivedGraphsFilter.hpp"
#include "WeightSubgrMono/GraphTheoretic/DerivedGraphs.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"

#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphsFilter::DerivedGraphsFilter(const FixedData& fixed_data) {
  m_updater_pair_ptr = std::make_unique<DerivedGraphsUpdaterPair>(
      fixed_data.pattern_neighbours_data,
      fixed_data.target_neighbours_data,
      m_calculator,
      m_storage);
  TKET_ASSERT(m_updater_pair_ptr);
}

DerivedGraphsFilter::~DerivedGraphsFilter() {}

DerivedGraphs& DerivedGraphsFilter::get_derived_pattern_graphs() {
  return m_updater_pair_ptr->patterns_updater.get_derived_graphs();
}

DerivedGraphs& DerivedGraphsFilter::get_derived_target_graphs() {
  return m_updater_pair_ptr->targets_updater.get_derived_graphs();
}

static void fill_and_sort_weights(
      const DerivedGraphStructs::NeighboursAndCounts& raw_data,
      std::vector<DerivedGraphStructs::Count>& weights) {
  if(!weights.empty()) {
    TKET_ASSERT(weights.size() == raw_data.size());
    return;
  }
  weights.reserve(raw_data.size());
  for(const auto& entry : raw_data) {
    weights.push_back(entry.second);
  }
  std::sort(weights.begin(), weights.end());
}


template<class DGraph>
static bool compatible_using_weight_maps(VertexWSM pv, VertexWSM tv,
      DGraph& pattern_graph,
      std::map<VertexWSM, std::vector<DerivedGraphStructs::Count>>& pattern_weights_map,
      DGraph& target_graph,
      std::map<VertexWSM, std::vector<DerivedGraphStructs::Count>>& target_weights_map) {
  const auto& pattern_neighbours = pattern_graph.get_neighbours(pv);
  const auto& target_neighbours = target_graph.get_neighbours(tv);
  if(pattern_neighbours.size() > target_neighbours.size()) {
    return false;
  }
  auto& pattern_weights = pattern_weights_map[pv];
  fill_and_sort_weights(pattern_neighbours, pattern_weights);

  auto& target_weights = target_weights_map[tv];
  fill_and_sort_weights(target_neighbours, target_weights);
  return FilterUtils::compatible_sorted_degree_sequences(pattern_weights, target_weights);
}



bool DerivedGraphsFilter::is_compatible(VertexWSM pv, VertexWSM tv, const FixedData& fixed_data) {
  {
    // Do we already know that it's compatible?
    const auto domain_citer = m_compatible_assignments.find(pv);
    if(domain_citer != m_compatible_assignments.cend() &&
          domain_citer->second.count(tv) != 0) {
      return true;
    }
  }
  {
    // Do we already know that it's impossible?
    const auto domain_complement_citer = m_impossible_assignments.find(pv);
    if(domain_complement_citer != m_impossible_assignments.cend() &&
          domain_complement_citer->second.count(tv) != 0) {
      return false;
    }
  }

  // We don't know, we must calculate!
  // Break out of the loop if we detect an impossibility.
  for(;;) {
    if(get_derived_pattern_graphs().triangle_counts.get_count(pv) >
          get_derived_target_graphs().triangle_counts.get_count(tv)) {
      break;
    }
    if(!compatible_using_weight_maps(pv, tv, 
            get_derived_pattern_graphs().d2_graph,
            m_d2_pattern_weights,
            get_derived_target_graphs().d2_graph,
            m_d2_target_weights)) {
      break;
    }
    if(!compatible_using_weight_maps(pv, tv, 
            get_derived_pattern_graphs().d3_graph,
            m_d3_pattern_weights,
            get_derived_target_graphs().d3_graph,
            m_d3_target_weights)) {
      break;
    }
    m_compatible_assignments[pv].insert(tv);
    return true;
  }
  m_impossible_assignments[pv].insert(tv);
  return false;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
