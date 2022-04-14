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

#include "WeightSubgrMono/Reducing/HallSetDetector.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool HallSetDetector::Data::operator<(const Data& other) const {
  return domain_size > other.domain_size ||
         (domain_size == other.domain_size && pv < other.pv);
}

bool HallSetDetector::fill_new_domain_sizes_data(
    const PossibleAssignments& possible_assignments) {
  m_domain_sizes_data.clear();
  Data data;

  for (const auto& entry : possible_assignments) {
    data.domain_size = entry.second.size();
    if (data.domain_size == 0) {
      return false;
    }
    if (data.domain_size == 1) {
      continue;
    }
    data.pv = entry.first;
    m_domain_sizes_data.emplace_back(data);
  }
  return true;
}

bool HallSetDetector::fill_existing_domain_sizes_data(
    const PossibleAssignments& possible_assignments) {
  // DON'T increment ii, swap with the back to erase
  for (unsigned ii = 0; ii < m_domain_sizes_data.size();) {
    auto& element = m_domain_sizes_data[ii];
    element.domain_size = possible_assignments.at(element.pv).size();
    if (element.domain_size == 0) {
      return false;
    }
    if (element.domain_size > 1) {
      ++ii;
      continue;
    }
    // The element is to be "erased". We do this by swapping with the back,
    // then popping off the back element.
    if (ii + 1 < m_domain_sizes_data.size()) {
      element.pv = m_domain_sizes_data.back().pv;
    }
    m_domain_sizes_data.pop_back();
  }
  return true;
}

const HallSetDetector::Result& HallSetDetector::get_hall_set(
    const PossibleAssignments& possible_assignments, Action action) {
  m_result.pattern_vertices.clear();
  m_result.union_of_domains.clear();
  m_result.status = Status::UNINTERESTING;

  bool success;
  if (action == Action::CLEAR_DATA) {
    success = fill_new_domain_sizes_data(possible_assignments);
  } else {
    success = fill_existing_domain_sizes_data(possible_assignments);
  }
  if (!success) {
    m_result.status = Status::NOGOOD;
    return m_result;
  }
  if (m_domain_sizes_data.empty()) {
    return m_result;
  }
  // Now we try to detect a nogood.
  std::sort(m_domain_sizes_data.begin(), m_domain_sizes_data.end());
  fill_result_using_nonempty_domain_sizes_data(possible_assignments);
  return m_result;
}

void HallSetDetector::fill_result_using_nonempty_domain_sizes_data(
    const PossibleAssignments& possible_assignments) {
  // The number of PV, starting from back(),
  // which were used to create the combined domain (the union).
  size_t current_number_of_pv = 0;
  TKET_ASSERT(m_result.pattern_vertices.empty());
  TKET_ASSERT(m_result.union_of_domains.empty());
  TKET_ASSERT(m_result.status == Status::UNINTERESTING);

  // Break only upon an interesting result; return if uninteresting
  // (when we are SURE that it will be uninteresting).
  for (;;) {
    size_t current_union_size = m_result.union_of_domains.size();

    if (current_number_of_pv > 0 &&
        current_union_size <= current_number_of_pv) {
      break;
    }
    if (current_number_of_pv >= m_domain_sizes_data.size()) {
      return;
    }
    // The combined domain size must be <= the number of pv in the data,
    // or it is not interesting.
    // Currently the size is > the number of pv, so we must keep adding more pv
    // in the hope that the number of pv overtakes the union size.

    // Check if it's POSSIBLE that adding more PVs could allow
    // the number of PVs to overtake the union size.
    bool continue_adding_pvs = false;
    auto union_size_lower_bound = current_union_size;

    for (auto number_of_pv = current_number_of_pv + 1;
         number_of_pv < m_domain_sizes_data.size(); ++number_of_pv) {
      // If we took the union of the current combined domain with
      // "number_of_pv", the size would definitely be >= this.
      union_size_lower_bound = std::max(
          current_union_size,
          m_domain_sizes_data[(m_domain_sizes_data.size() - number_of_pv) - 1]
              .domain_size);

      if (union_size_lower_bound <= number_of_pv) {
        // It's POSSIBLE that adding up to this many PVs would give
        // a small enough union.
        continue_adding_pvs = true;
        break;
      }
      // If we reached here, we can be SURE that the union would be too large.
    }
    if (!continue_adding_pvs) {
      // We are sure that adding ANY number of extra PVs
      // (at least, in the order we're considering - smallest domains first)
      // would give a union too large. So give up.
      return;
    }
    // Based upon our size estimate, there's the POSSIBILITY that
    // the union might be small enough if we add a certain additional number of
    // PVs. So just add one PV and re-evaluate.
    const VertexWSM new_pv =
        m_domain_sizes_data
            [(m_domain_sizes_data.size() - current_number_of_pv) - 1]
                .pv;
    const auto& new_domain = possible_assignments.at(new_pv);
    for (VertexWSM tv : new_domain) {
      m_result.union_of_domains.insert(tv);
      if (m_result.union_of_domains.size() > m_domain_sizes_data.size()) {
        return;
      }
    }
    ++current_number_of_pv;
  }
  fill_interesting_result(current_number_of_pv);
}

void HallSetDetector::fill_interesting_result(std::size_t number_of_pv) {
  if (m_result.union_of_domains.size() < number_of_pv) {
    m_result.status = Status::NOGOOD;
    return;
  }
  TKET_ASSERT(m_result.union_of_domains.size() == number_of_pv);
  TKET_ASSERT(m_result.pattern_vertices.empty());
  TKET_ASSERT(number_of_pv <= m_domain_sizes_data.size());
  m_result.pattern_vertices.reserve(number_of_pv);

  m_result.status = Status::HALL_SET;
  // We actually need the PV, so fill them in.
  for (std::size_t ii = 0; ii < number_of_pv; ++ii) {
    m_result.pattern_vertices.push_back(
        m_domain_sizes_data[m_domain_sizes_data.size() - ii - 1].pv);
  }
  std::sort(m_result.pattern_vertices.begin(), m_result.pattern_vertices.end());

  // Clear out the vertices used, which are all at the back,
  // ready for the next possible Hall set.
  m_domain_sizes_data.resize(m_domain_sizes_data.size() - number_of_pv);
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
