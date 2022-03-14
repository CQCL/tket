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

#include "TokenSwapping/CanonicalRelabelling.hpp"

#include <algorithm>
#include <numeric>

#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

CanonicalRelabelling::CanonicalRelabelling() {
  // no more than 6 vertices, so no more than 6 cycles ever needed.
  m_cycles.resize(6);
}

const CanonicalRelabelling::Result& CanonicalRelabelling::operator()(
    const VertexMapping& desired_mapping) {
  m_result.too_many_vertices = false;
  m_result.permutation_hash = 0;
  m_result.new_to_old_vertices.clear();
  m_result.old_to_new_vertices.clear();
  m_result.identity = all_tokens_home(desired_mapping);
  if (m_result.identity) {
    return m_result;
  }
  check_mapping(desired_mapping, m_work_mapping);
  if (desired_mapping.size() > 6) {
    m_result.too_many_vertices = true;
    return m_result;
  }
  // If not the identity, at least 2 vertices moved.
  TKET_ASSERT(desired_mapping.size() >= 2);
  TKET_ASSERT(desired_mapping.size() <= 6);

  m_desired_mapping = desired_mapping;
  unsigned next_cyc_index = 0;

  while (!m_desired_mapping.empty()) {
    // New cycle starts
    auto& this_cycle = m_cycles[next_cyc_index];
    ++next_cyc_index;
    this_cycle.clear();
    this_cycle.push_back(m_desired_mapping.cbegin()->first);
    bool terminated_correctly = false;
    for (unsigned infinite_loop_guard = 1 + m_desired_mapping.size();
         infinite_loop_guard != 0; --infinite_loop_guard) {
      const auto curr_v = this_cycle.back();
      const auto target_v = m_desired_mapping.at(curr_v);
      TKET_ASSERT(m_desired_mapping.erase(curr_v) == 1);
      if (target_v == this_cycle[0]) {
        terminated_correctly = true;
        break;
      }
      this_cycle.push_back(target_v);
    }
    TKET_ASSERT(terminated_correctly);
  }
  // Sort by cycle length, LONGEST cycles first.
  // But, also want a "stable-like" sort:
  // make a consistent choice across all platforms,
  // if cycle lengths are equal,
  // based only upon the vertex numbers.
  m_sorted_cycles_indices.resize(next_cyc_index);
  std::iota(m_sorted_cycles_indices.begin(), m_sorted_cycles_indices.end(), 0);
  const auto& cycles = m_cycles;

  std::sort(
      m_sorted_cycles_indices.begin(), m_sorted_cycles_indices.end(),
      [cycles](unsigned ii, unsigned jj) {
        const auto& cyc1 = cycles[ii];
        const auto& cyc2 = cycles[jj];
        return (cyc1.size() > cyc2.size()) ||
               // Using the raw vertex numbers is, of course, non-canonical,
               // but necessary if we are to have stable results
               // across ALL nonstable sorting algorithms
               // on different platforms/compilers.
               ((cyc1.size() == cyc2.size()) && cyc1[0] < cyc2[0]);
      });

  // Now we can set up the mapping.
  m_result.new_to_old_vertices.clear();
  for (auto ii : m_sorted_cycles_indices) {
    const auto& cyc = m_cycles[ii];
    TKET_ASSERT(!cyc.empty());
    TKET_ASSERT(cyc.size() <= 6);
    for (size_t old_v : cyc) {
      m_result.new_to_old_vertices.push_back(old_v);
    }
  }
  TKET_ASSERT(m_result.new_to_old_vertices.size() <= 6);
  m_result.old_to_new_vertices.clear();
  for (unsigned ii = 0; ii < m_result.new_to_old_vertices.size(); ++ii) {
    m_result.old_to_new_vertices[m_result.new_to_old_vertices[ii]] = ii;
  }
  // GCOVR_EXCL_START
  TKET_ASSERT(
      m_result.new_to_old_vertices.size() ==
      m_result.old_to_new_vertices.size());
  // GCOVR_EXCL_STOP
  // And finally, the permutation hash.
  m_result.permutation_hash = 0;
  for (auto ii : m_sorted_cycles_indices) {
    const auto& cyc = m_cycles[ii];
    if (cyc.size() == 1) {
      break;
    }
    m_result.permutation_hash *= 10;
    m_result.permutation_hash += cyc.size();
  }
  return m_result;
}

}  // namespace tsa_internal
}  // namespace tket
