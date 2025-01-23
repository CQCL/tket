// Copyright Quantinuum
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

#include "tkwsm/InitPlacement/SolutionJumper.hpp"

#include <sstream>
#include <stdexcept>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"
#include "tkwsm/InitPlacement/UtilsIQP.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

SolutionJumper::SolutionJumper(
    const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
    WeightWSM implicit_target_weight)
    : m_pattern_ndata(pattern_ndata),
      m_target_ndata(target_ndata),
      m_implicit_target_weight(implicit_target_weight) {
  m_assigned_target_vertices.resize(
      m_pattern_ndata.get_number_of_nonisolated_vertices());

  get_scalar_product_upper_bound_for_complete_target_graph(
      m_pattern_ndata, m_target_ndata, m_implicit_target_weight);
}

const std::vector<unsigned>& SolutionJumper::get_assignments() const {
  return m_assigned_target_vertices;
}

std::vector<unsigned>& SolutionJumper::get_assignments_to_overwrite() {
  // Invalidate the weight contributions; will be recalculated lazily.
  m_scalar_product_contributions.resize(m_assigned_target_vertices.size());
  for (unsigned pv = 0; pv < m_scalar_product_contributions.size(); ++pv) {
    m_scalar_product_contributions[pv] = 0;
  }
  return m_assigned_target_vertices;
}

const NeighboursData& SolutionJumper::get_pattern_ndata() const {
  return m_pattern_ndata;
}

const NeighboursData& SolutionJumper::get_target_ndata() const {
  return m_target_ndata;
}

WeightWSM SolutionJumper::reset_and_get_new_scalar_product() {
  TKET_ASSERT(
      m_assigned_target_vertices.size() ==
      m_pattern_ndata.get_number_of_nonisolated_vertices());
  m_source_pattern_vertices.resize(
      m_target_ndata.get_number_of_nonisolated_vertices());

  // Don't need sets, just fill with dummy values (a lot quicker!)
  for (auto& pv : m_source_pattern_vertices) {
    set_maximum(pv);
  }

  for (unsigned pv = 0; pv < m_assigned_target_vertices.size(); ++pv) {
    const auto tv = m_assigned_target_vertices[pv];
    if (tv >= m_source_pattern_vertices.size()) {
      throw std::runtime_error("PV assigned to invalid TV");
    }
    auto& pv_to_set = m_source_pattern_vertices[tv];
    if (pv_to_set < m_assigned_target_vertices.size()) {
      throw std::runtime_error("PV assigned to multiple TV");
    }
    pv_to_set = pv;
  }

  // We have to recalculate all the contributions
  // (which also gives us the total scalar product).
  WeightWSM total_contribution = 0;
  m_scalar_product_contributions.resize(m_assigned_target_vertices.size());
  for (unsigned pv = 0; pv < m_scalar_product_contributions.size(); ++pv) {
    m_scalar_product_contributions[pv] =
        get_hypothetical_scalar_product_contribution_disallowing_case_c(
            pv, m_assigned_target_vertices.at(pv));
    total_contribution += m_scalar_product_contributions[pv];
  }
  TKET_ASSERT(
      total_contribution / 2 == get_scalar_product_with_complete_target(
                                    m_pattern_ndata, m_target_ndata,
                                    m_implicit_target_weight,
                                    get_assignments()));
  // Every edge was counted twice.
  return total_contribution / 2;
}

WeightWSM SolutionJumper::get_current_scalar_product_contribution(
    unsigned pv) const {
  auto& contribution = m_scalar_product_contributions.at(pv);
  if (contribution == 0) {
    contribution =
        get_hypothetical_scalar_product_contribution_disallowing_case_c(
            pv, m_assigned_target_vertices.at(pv));
  }
  return contribution;
}

std::optional<WeightWSM>
SolutionJumper::perform_move_and_get_scalar_product_decrease(
    unsigned pv1, unsigned tv2, WeightWSM minimum_decrease) {
  check_validity();

  const unsigned pv2 = m_source_pattern_vertices.at(tv2);
  if (pv1 == pv2) {
    // No actual change!
    return {};
  }

  const unsigned tv1 = m_assigned_target_vertices.at(pv1);
  TKET_ASSERT(tv1 != tv2);

  WeightWSM existing_contribution_to_erase =
      get_current_scalar_product_contribution(pv1);

  const auto pv1_new_contribution_data =
      get_hypothetical_scalar_product_contribution(pv1, tv2);

  WeightWSM new_contribution_to_add = pv1_new_contribution_data.contribution;

  if (pv2 < m_assigned_target_vertices.size()) {
    // We currently have PV1 -> TV1, and PV2 -> TV2.
    existing_contribution_to_erase +=
        get_current_scalar_product_contribution(pv2);

    const auto pv2_new_contribution_data =
        get_hypothetical_scalar_product_contribution(pv2, tv1);
    new_contribution_to_add += pv2_new_contribution_data.contribution;

    if (existing_contribution_to_erase <
        new_contribution_to_add + minimum_decrease) {
      return {};
    }
    invalidate_neighbour_contributions(pv1);
    invalidate_neighbour_contributions(pv2);

    // Make the new assignment PV2->TV1.
    m_scalar_product_contributions.at(pv2) =
        pv2_new_contribution_data.contribution;
    m_assigned_target_vertices.at(pv2) = tv1;

    // Check for case C. Doesn't actually matter, just a validity check.
    if (pv1_new_contribution_data.case_c_other_pv_opt) {
      TKET_ASSERT(pv2_new_contribution_data.case_c_other_pv_opt);
      TKET_ASSERT(pv1_new_contribution_data.case_c_other_pv_opt.value() == pv2);
      TKET_ASSERT(pv2_new_contribution_data.case_c_other_pv_opt.value() == pv1);
    } else {
      TKET_ASSERT(!pv2_new_contribution_data.case_c_other_pv_opt);
    }
  } else {
    // TV2 is unoccupied; there is no PV2 mapping to it!
    TKET_ASSERT(!pv1_new_contribution_data.case_c_other_pv_opt);
    if (existing_contribution_to_erase <
        new_contribution_to_add + minimum_decrease) {
      return {};
    }
    invalidate_neighbour_contributions(pv1);
  }
  m_source_pattern_vertices.at(tv1) = pv2;

  // Make the new assignment PV1->TV2.
  m_scalar_product_contributions.at(pv1) =
      pv1_new_contribution_data.contribution;
  m_assigned_target_vertices.at(pv1) = tv2;
  m_source_pattern_vertices.at(tv2) = pv1;
  return existing_contribution_to_erase - new_contribution_to_add;
}

WeightWSM SolutionJumper::get_target_edge_weight(
    unsigned tv1, unsigned tv2) const {
  const auto t_edge_weight_opt = m_target_ndata.get_edge_weight_opt(tv1, tv2);
  if (t_edge_weight_opt) {
    return t_edge_weight_opt.value();
  }
  return m_implicit_target_weight;
}

void SolutionJumper::check_validity() const {
  if (m_assigned_target_vertices.empty() || m_source_pattern_vertices.empty()) {
    return;
  }
  for (unsigned pv = 0; pv < m_assigned_target_vertices.size(); ++pv) {
    TKET_ASSERT(
        m_source_pattern_vertices.at(m_assigned_target_vertices[pv]) == pv);
  }
  for (unsigned tv = 0; tv < m_source_pattern_vertices.size(); ++tv) {
    const auto pv = m_source_pattern_vertices[tv];
    TKET_ASSERT(
        pv >= m_assigned_target_vertices.size() ||
        m_assigned_target_vertices[pv] == tv);
  }
}

WeightWSM
SolutionJumper::get_hypothetical_scalar_product_contribution_disallowing_case_c(
    unsigned pv, unsigned tv) const {
  const auto result = get_hypothetical_scalar_product_contribution(pv, tv);
  TKET_ASSERT(!result.case_c_other_pv_opt);
  return result.contribution;
}

SolutionJumper::HypotheticalScalarProductContribution
SolutionJumper::get_hypothetical_scalar_product_contribution(
    unsigned pv, unsigned tv) const {
  HypotheticalScalarProductContribution result;
  result.contribution = 0;

  for (const auto& pv_weight_pair :
       m_pattern_ndata.get_neighbours_and_weights(pv)) {
    const unsigned& other_pv = pv_weight_pair.first;
    const unsigned& other_tv = m_assigned_target_vertices.at(other_pv);

    const WeightWSM& p_edge_weight = pv_weight_pair.second;
    WeightWSM t_edge_weight;
    if (tv == other_tv) {
      t_edge_weight =
          get_target_edge_weight(tv, m_assigned_target_vertices.at(pv));
      TKET_ASSERT(!result.case_c_other_pv_opt);
      result.case_c_other_pv_opt = other_pv;
    } else {
      t_edge_weight = get_target_edge_weight(tv, other_tv);
    }

    result.contribution = get_sum_or_throw(
        result.contribution,
        get_product_or_throw(p_edge_weight, t_edge_weight));
  }
  return result;
}

void SolutionJumper::invalidate_neighbour_contributions(unsigned pv) {
  for (const auto& p_neighbour_weight_pair :
       m_pattern_ndata.get_neighbours_and_weights(pv)) {
    m_scalar_product_contributions[p_neighbour_weight_pair.first] = 0;
  }
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
