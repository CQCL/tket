// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tkwsm/WeightPruning/WeightChecker.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Common/TemporaryRefactorCode.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"
#include "tkwsm/Searching/SearchBranch.hpp"
#include "tkwsm/WeightPruning/WeightNogoodDetector.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

WeightChecker::WeightChecker(
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data,
    const SearchBranch& search_branch, WeightWSM total_p_edge_weights,
    std::set<VertexWSM>& impossible_target_vertices)
    : m_pattern_neighbours_data(pattern_neighbours_data),
      m_target_neighbours_data(target_neighbours_data),
      m_search_branch(search_branch),
      m_manager(total_p_edge_weights),
      m_impossible_target_vertices(impossible_target_vertices) {}

WeightChecker::~WeightChecker() {}

std::optional<WeightChecker::TVData> WeightChecker::get_tv_data_opt() {
  if (!m_detector_ptr) {
    return {};
  }
  m_tv_data.final_number_of_tv = m_detector_ptr->get_number_of_possible_tv();
  return m_tv_data;
}

bool WeightChecker::check(
    const DomainsAccessor& accessor, WeightWSM max_extra_scalar_product) {
  std::size_t current_number_of_unassigned_vertices = 0;

  for (VertexWSM pv : accessor.get_unassigned_pattern_vertices_superset()) {
    switch (accessor.get_domain_size(pv)) {
      // TODO: make a test case for this; but it's fiddly.
      // GCOVR_EXCL_START
      case 0:
        return false;
      // GCOVR_EXCL_STOP
      case 1:
        break;
      default:
        ++current_number_of_unassigned_vertices;
    }
  }
  const std::size_t current_number_of_assigned_vertices =
      accessor.get_number_of_pattern_vertices() -
      current_number_of_unassigned_vertices;

  const auto current_scalar_product = accessor.get_scalar_product();
  const auto max_scalar_product =
      current_scalar_product + max_extra_scalar_product;

  if (!m_manager.should_activate_detector(
          current_scalar_product, max_scalar_product,
          accessor.get_total_p_edge_weights(),
          current_number_of_assigned_vertices,
          current_number_of_unassigned_vertices)) {
    return true;
  }

  // Finally, we use the detector; check if it's initialised.
  if (!m_detector_ptr) {
    std::set<VertexWSM> used_tv;
    TemporaryRefactorCode::set_domain_from_bitset(
        used_tv, m_search_branch.get_used_target_vertices());
    m_tv_data.initial_number_of_tv = used_tv.size();

    m_detector_ptr = std::make_unique<WeightNogoodDetector>(
        m_pattern_neighbours_data, m_target_neighbours_data, used_tv,
        m_impossible_target_vertices);
    TKET_ASSERT(m_detector_ptr);
  }
  const auto extra_scalar_product_lower_bound_opt =
      m_detector_ptr->get_extra_scalar_product_lower_bound(
          accessor, max_extra_scalar_product);

  if (extra_scalar_product_lower_bound_opt) {
    // We have a lower bound, but it's not big enough to force the
    // scalar product beyond the limit; thus we haven't found a weight nogood,
    // so this counts as a failure.
    m_manager.register_lower_bound_failure(
        current_scalar_product, max_scalar_product,
        extra_scalar_product_lower_bound_opt.value());
    return true;
  }
  // We HAVE found a nogood, which is a SUCCESS
  // (since it's the whole purpose of the detector).
  m_manager.register_success();
  return false;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
