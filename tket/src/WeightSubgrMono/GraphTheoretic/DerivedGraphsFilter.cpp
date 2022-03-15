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
#include "WeightSubgrMono/Searching/FixedData.hpp"

#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphsContainer& DerivedGraphsFilter::get_container() {
  return m_container;
}


static void fill_sorted_counts(const DerivedGraphsCalculator::NeighboursAndCounts& neighbours,
    std::vector<std::size_t>& counts) {
  if(counts.size() == neighbours.size()) {
    return;
  }
  counts.resize(neighbours.size());
  for(unsigned ii=0; ii<neighbours.size(); ++ii) {
    counts[ii] = neighbours[ii].second;
  }
  std::sort(counts.begin(), counts.end());
}

template<class GetDegreeFunc>
static void fill_degree_sequence(
      const DerivedGraphsCalculator::NeighboursAndCounts& neighbours,
      std::vector<std::size_t>& degree_sequence,
      const GetDegreeFunc& func) {
  if(degree_sequence.size() == neighbours.size()) {
    return;
  }
  degree_sequence.resize(neighbours.size());
  for(unsigned ii=0; ii<neighbours.size(); ++ii) {
    degree_sequence[ii] = func(neighbours[ii].first);
  }
  std::sort(degree_sequence.begin(), degree_sequence.end());
}


bool DerivedGraphsFilter::is_compatible(VertexWSM pv, VertexWSM tv, const FixedData& fixed_data) {
  {
    const auto domain_citer = m_compatible_assignments.find(pv);
    if(domain_citer != m_compatible_assignments.cend() &&
          domain_citer->second.count(tv) != 0) {
      return true;
    }
  }
  {
    const auto domain_complement_citer = m_impossible_assignments.find(pv);
    if(domain_complement_citer != m_impossible_assignments.cend() &&
          domain_complement_citer->second.count(tv) != 0) {
      return false;
    }
  }
  // It must be calculated.
  auto& pv_data = m_container.get_pattern_v_data_permanent_reference(pv, fixed_data.pattern_neighbours_data);
  auto& tv_data = m_container.get_target_v_data_permanent_reference(tv, fixed_data.target_neighbours_data);

  // Do increaingly expensive checks.
  // Break out of the loop if invalid.
  for(;;) {
    if(pv_data.triangle_count > tv_data.triangle_count ||
          pv_data.depth_2_neighbours.size() > tv_data.depth_2_neighbours.size() ||
          pv_data.depth_3_neighbours.size() > tv_data.depth_3_neighbours.size()) {
            
      break;
    }
    // Match the counts, i.e. edge weights in the derived graph.
    
    fill_sorted_counts(pv_data.depth_2_neighbours, pv_data.depth_2_counts);
    fill_sorted_counts(tv_data.depth_2_neighbours, tv_data.depth_2_counts);
    if(!FilterUtils::compatible_sorted_degree_sequences(pv_data.depth_2_counts, tv_data.depth_2_counts)) {
      break;
    }
    fill_sorted_counts(pv_data.depth_3_neighbours, pv_data.depth_3_counts);
    fill_sorted_counts(tv_data.depth_3_neighbours, tv_data.depth_3_counts);
    if(!FilterUtils::compatible_sorted_degree_sequences(pv_data.depth_3_counts, tv_data.depth_3_counts)) {
      break;
    }
    // Now, check the degree sequences in the derived graphs themselves.
    auto& container = m_container;
    fill_degree_sequence(pv_data.depth_2_neighbours, pv_data.depth_2_degree_sequence,
        [pv, &container, &fixed_data](VertexWSM other_v) {
          return container.get_pattern_v_data_permanent_reference(
              other_v, fixed_data.pattern_neighbours_data).depth_2_neighbours.size();
        });

    fill_degree_sequence(tv_data.depth_2_neighbours, tv_data.depth_2_degree_sequence,
        [tv, &container, &fixed_data](VertexWSM other_v) {
          return container.get_target_v_data_permanent_reference(
              other_v, fixed_data.target_neighbours_data).depth_2_neighbours.size();
        });
    
    if(!FilterUtils::compatible_sorted_degree_sequences(pv_data.depth_2_degree_sequence, tv_data.depth_2_degree_sequence)) {
      break;
    }

    // Finally, depth 3. A bit ugly to repeat similar code,
    // but not worth fancy stuff to deduplicate (we are NOT doing
    // deeper levels of iterated derived graphs; experiments showed that
    // it took longer to compute than the time saved).
    fill_degree_sequence(pv_data.depth_3_neighbours, pv_data.depth_3_degree_sequence,
        [pv, &container, &fixed_data](VertexWSM other_v) {
          return container.get_pattern_v_data_permanent_reference(
              other_v, fixed_data.pattern_neighbours_data).depth_3_neighbours.size();
        });

    fill_degree_sequence(tv_data.depth_3_neighbours, tv_data.depth_3_degree_sequence,
        [tv, &container, &fixed_data](VertexWSM other_v) {
          return container.get_target_v_data_permanent_reference(
              other_v, fixed_data.target_neighbours_data).depth_3_neighbours.size();
        });
    
    if(!FilterUtils::compatible_sorted_degree_sequences(pv_data.depth_3_degree_sequence, tv_data.depth_3_degree_sequence)) {
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
