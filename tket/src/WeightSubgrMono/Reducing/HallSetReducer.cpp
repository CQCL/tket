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

#include "WeightSubgrMono/Reducing/HallSetReducer.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/SetIntersection.hpp"
#include "WeightSubgrMono/Searching/NodeWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// NOTE:
// If a Hall set is detected in node[i],
// it might seem like a good idea to check node[i-1], node[i-2], ...
// also in the search branch, if the check is cheap.
// We need to experiment and see. Maybe it seems a little unlikely
// to be useful, since domain sizes are (mostly) monotonic,
// and it might seem that it would have been detected sooner.
// On the other hand, the Hall set detector might fail to detect
// many sets, purely down to a slight reordering of variables,
// so maybe our chances of finding the Hall set in previous nodes is higher.


HallSetReducer::HallSetReducer() : m_awaiting_first_reduction(true) {}

void HallSetReducer::clear() { m_awaiting_first_reduction = true; }

HallSetReducer::Result HallSetReducer::operator()(NodeWSM& node) {
  HallSetDetector::Action action;
  if (m_awaiting_first_reduction) {
    action = HallSetDetector::Action::CLEAR_DATA;
    m_awaiting_first_reduction = false;
  } else {
    action = HallSetDetector::Action::USE_EXISTING_DATA;
  }
  const auto& current_domains = node.get_possible_assignments();
  const auto& detector_result =
      m_detector.get_hall_set(current_domains, action);
  switch (detector_result.status) {
    case HallSetDetector::Status::UNINTERESTING:
      return Result::FINISHED;
    case HallSetDetector::Status::NOGOOD:
      return Result::NOGOOD;
    case HallSetDetector::Status::HALL_SET:
      break;
  }
  // We've got a Hall set.
  const auto fill_result =
      fill_new_domains_data_from_hall_set(current_domains, detector_result);

  if (fill_result.nogood) {
    return Result::NOGOOD;
  }
  if (!fill_result.changed) {
    return Result::FINISHED;
  }
  TKET_ASSERT(fill_result.number_of_new_domains <= m_new_domains_data.size());
  for (unsigned ii = 0; ii < fill_result.number_of_new_domains; ++ii) {
    node.overwrite_domain(
        m_new_domains_data[ii].first, m_new_domains_data[ii].second);
  }
  if (fill_result.new_assignments) {
    return Result::NEW_ASSIGNMENTS;
  }
  // We've reduced, but not got new assignments; so recurse.
  return operator()(node);
}

HallSetReducer::FillResult HallSetReducer::fill_new_domains_data_from_hall_set(
    const PossibleAssignments& current_domains,
    const HallSetDetector::Result& detector_result) {
  TKET_ASSERT(
      detector_result.pattern_vertices.size() ==
      detector_result.union_of_domains.size());

  FillResult result;
  result.number_of_new_domains = 0;
  result.new_assignments = false;
  result.nogood = false;
  result.changed = false;

  for (const auto& entry : current_domains) {
    const VertexWSM& pv = entry.first;
    if (std::binary_search(
            detector_result.pattern_vertices.cbegin(),
            detector_result.pattern_vertices.cend(), pv)) {
      continue;
    }
    const auto& current_domain = entry.second;
    if (disjoint(current_domain, detector_result.union_of_domains)) {
      continue;
    }
    // Now we find the new domain; it is definitely different.
    result.changed = true;
    if (result.number_of_new_domains >= m_new_domains_data.size()) {
      m_new_domains_data.resize(result.number_of_new_domains + 1);
    }
    auto& new_domain_data = m_new_domains_data[result.number_of_new_domains];
    auto& new_domain = new_domain_data.second;
    ++result.number_of_new_domains;

    new_domain.clear();
    for (VertexWSM tv : current_domain) {
      if (detector_result.union_of_domains.count(tv) == 0) {
        new_domain.push_back(tv);
      }
    }
    if (new_domain.empty()) {
      result.nogood = true;
      return result;
    }
    new_domain_data.first = pv;
    if (new_domain.size() == 1) {
      result.new_assignments = true;
    }
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
