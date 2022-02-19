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

#include "CyclesCandidateManager.hpp"

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <stdexcept>

#include "Utils/Assert.hpp"
#include "VertexSwapResult.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

size_t CyclesCandidateManager::fill_initial_cycle_ids(const Cycles& cycles) {
  m_cycle_with_vertex_hash.clear();
  m_cycles_to_keep.clear();
  size_t cycle_length = 0;
  for (auto id_opt = cycles.front_id(); id_opt;
       id_opt = cycles.next(id_opt.value())) {
    const auto& cycle = cycles.at(id_opt.value());
    const auto& vertices = cycle.vertices;

    if (cycle_length == 0) {
      cycle_length = vertices.size();
      TKET_ASSERT(cycle_length >= 2);
    } else {
      TKET_ASSERT(cycle_length == vertices.size());
    }
    TKET_ASSERT(cycle.decrease > 0);

    // We want 50*(decrease)/(num swaps) >=  min_candidate_power_percentage.
    // (We multiply by 50 because a swap can change L by 2, not 1).
    if (50 * static_cast<unsigned>(cycle.decrease) <
        (m_options.min_candidate_power_percentage * cycle_length)) {
      continue;
    }

    CycleData cycle_data;
    cycle_data.id = id_opt.value();
    cycle_data.first_vertex_index = 0;
    for (size_t ii = 1; ii < vertices.size(); ++ii) {
      if (vertices[ii] < vertices[cycle_data.first_vertex_index]) {
        cycle_data.first_vertex_index = ii;
      }
    }
    size_t hash = static_cast<size_t>(cycle.decrease);
    for (size_t ii = 0; ii < cycle_length; ++ii) {
      boost::hash_combine(
          hash, vertices[(ii + cycle_data.first_vertex_index) % cycle_length]);
    }
    const auto prev_cycle_citer = m_cycle_with_vertex_hash.find(hash);
    if (prev_cycle_citer == m_cycle_with_vertex_hash.cend()) {
      m_cycle_with_vertex_hash[hash] = cycle_data;
    } else {
      // A previous cycle with this hash; but is it equal?
      const auto& previous_cycle_data = prev_cycle_citer->second;
      const auto& previous_cycle = cycles.at(previous_cycle_data.id);
      if (previous_cycle.decrease == cycle.decrease) {
        bool equal_vertices = true;
        for (size_t ii = 0; ii < cycle_length; ++ii) {
          if (previous_cycle.vertices
                  [(ii + previous_cycle_data.first_vertex_index) %
                   cycle_length] !=
              cycle.vertices
                  [(ii + cycle_data.first_vertex_index) % cycle_length]) {
            equal_vertices = false;
            break;
          }
        }
        if (equal_vertices) {
          // This new cycle is just the previous cycle repeated,
          // but starting from a different vertex
          continue;
        }
      }
    }
    m_cycles_to_keep.push_back(cycle_data.id);
  }
  return cycle_length;
}

void CyclesCandidateManager::discard_lower_power_solutions(
    const Cycles& cycles) {
  int highest_decrease = 0;
  for (auto id : m_cycles_to_keep) {
    highest_decrease = std::max(highest_decrease, cycles.at(id).decrease);
  }
  TKET_ASSERT(highest_decrease > 0);

  for (size_t ii = 0; ii < m_cycles_to_keep.size();) {
    if (cycles.at(m_cycles_to_keep[ii]).decrease < highest_decrease) {
      // This cycle is not good enough.
      // Erase this ID, by swapping with the back
      m_cycles_to_keep[ii] = m_cycles_to_keep.back();
      m_cycles_to_keep.pop_back();
      continue;
    }
    // Keep this ID. Onto the next!
    ++ii;
  }
}

void CyclesCandidateManager::sort_candidates(const Cycles& cycles) {
  // Greedy heuristic: we want the maximal number of disjoint cycles.
  // So, choose those which touch few others first.
  // Experimentation is needed with other algorithms!
  m_touching_data.clear();
  for (size_t ii = 0; ii < m_cycles_to_keep.size(); ++ii) {
    // Automatically set to zero on first use.
    m_touching_data[m_cycles_to_keep[ii]];

    for (size_t jj = ii + 1; jj < m_cycles_to_keep.size(); ++jj) {
      bool touches = false;
      // For short cycles, not much slower than using sets
      // or sorted vectors.
      for (auto v1 : cycles.at(m_cycles_to_keep[ii]).vertices) {
        if (touches) {
          break;
        }
        for (auto v2 : cycles.at(m_cycles_to_keep[jj]).vertices) {
          if (v1 == v2) {
            touches = true;
            break;
          }
        }
      }
      if (touches) {
        ++m_touching_data[m_cycles_to_keep[ii]];
        ++m_touching_data[m_cycles_to_keep[jj]];
      }
    }
  }
  // Now, sort...
  auto& touching_data = m_touching_data;
  std::sort(
      m_cycles_to_keep.begin(), m_cycles_to_keep.end(),
      [&touching_data](Cycles::ID lhs, Cycles::ID rhs) {
        const auto lhs_touch_number = touching_data.at(lhs);
        const auto rhs_touch_number = touching_data.at(rhs);

        // Don't JUST sort on the touch number, because then the order
        // of equal-touch-number elements would be implementation dependent
        // (i.e., not a "stable" sort across all platforms/compilers).
        return (lhs_touch_number < rhs_touch_number) ||
               (lhs_touch_number == rhs_touch_number && lhs < rhs);
      });
}

bool CyclesCandidateManager::should_add_swaps_for_candidate(
    const Cycles& cycles, Cycles::ID id) {
  const auto& cycle = cycles.at(id);
  const auto& vertices = cycle.vertices;
  for (auto v : vertices) {
    if (m_vertices_used.count(v) != 0) {
      return false;
    }
  }
  for (auto v : vertices) {
    m_vertices_used.insert(v);
  }
  return true;
}

void CyclesCandidateManager::append_partial_solution(
    const CyclesGrowthManager& growth_manager, SwapList& swaps,
    VertexMapping& vertex_mapping) {
  const auto& cycles = growth_manager.get_cycles();
  const size_t cycle_size = fill_initial_cycle_ids(cycles);

  if (m_cycles_to_keep.empty()) {
    return;
  }
  const bool keep_lower_power_solutions =
      (cycle_size == 2)
          ? m_options.return_all_good_single_swaps
          : m_options.return_lower_power_solutions_for_multiswap_candidates;

  if (!keep_lower_power_solutions) {
    discard_lower_power_solutions(cycles);
  }
  sort_candidates(cycles);
  m_vertices_used.clear();

  // It's the final function, so don't bother erasing
  // elements in m_cycles_to_keep.
  for (auto id : m_cycles_to_keep) {
    if (!should_add_swaps_for_candidate(cycles, id)) {
      continue;
    }
    const auto& vertices = cycles.at(id).vertices;
    for (size_t ii = vertices.size() - 1; ii > 0; --ii) {
      VertexSwapResult(vertices[ii], vertices[ii - 1], vertex_mapping, swaps);
    }
  }
}

}  // namespace tsa_internal
}  // namespace tket
