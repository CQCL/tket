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

#include "WeightSubgrMono/WeightPruning/WeightChecker.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"
#include "WeightSubgrMono/Searching/DomainsAccessor.hpp"
#include "WeightSubgrMono/Searching/SearchBranch.hpp"
#include "WeightSubgrMono/WeightPruning/WeightNogoodDetector.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

WeightChecker::WeightChecker(
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data,
    const SearchBranch& search_branch, WeightWSM total_p_edge_weights)
    : m_pattern_neighbours_data(pattern_neighbours_data),
      m_target_neighbours_data(target_neighbours_data),
      m_search_branch(search_branch),
      m_manager(total_p_edge_weights) {}

WeightChecker::~WeightChecker() {}

WeightChecker::Result WeightChecker::operator()(
    const DomainsAccessor& accessor, WeightWSM max_extra_scalar_product) {
  Result result;
  std::size_t current_number_of_unassigned_vertices = 0;

  for (VertexWSM pv : accessor.get_unassigned_pattern_vertices_superset()) {
    switch (accessor.get_domain(pv).size()) {
      case 0:
        result.nogood = true;
        return result;
      case 1:
        break;
      default:
        ++current_number_of_unassigned_vertices;
    }
  }
  const std::size_t current_number_of_assigned_vertices =
      accessor.get_pattern_vertices().size() -
      current_number_of_unassigned_vertices;

  const auto current_scalar_product = accessor.get_scalar_product();
  const auto max_scalar_product =
      current_scalar_product + max_extra_scalar_product;

  if (!m_manager.should_activate_detector(
          current_scalar_product, max_scalar_product,
          accessor.get_total_p_edge_weights(),
          current_number_of_assigned_vertices,
          current_number_of_unassigned_vertices)) {
    result.nogood = false;
    return result;
  }

  // Finally, we use the detector; check if it's initialised.
  if (!m_detector_ptr) {
    m_detector_ptr = std::make_unique<WeightNogoodDetector>(
        m_pattern_neighbours_data, m_target_neighbours_data,
        m_search_branch.get_used_target_vertices());
    TKET_ASSERT(m_detector_ptr);
  }
  const auto detector_result =
      m_detector_ptr->operator()(accessor, max_extra_scalar_product);
  result.invalid_t_vertex = detector_result.invalid_t_vertex;

  if (detector_result.extra_scalar_product_lower_bound) {
    result.nogood = false;
    m_manager.register_lower_bound_failure(
        current_scalar_product, max_scalar_product,
        detector_result.extra_scalar_product_lower_bound.value());
    return result;
  }
  result.nogood = true;
  m_manager.register_success();
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
