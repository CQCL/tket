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

#include "CyclesGrowthManager.hpp"

#include <stdexcept>

#include "TokenSwapping/DistanceFunctions.hpp"
#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

bool Cycle::contains(size_t vertex) const {
  for (auto vv : vertices) {
    if (vertex == vv) {
      return true;
    }
  }
  return false;
}

CyclesGrowthManager::Options& CyclesGrowthManager::get_options() {
  return m_options;
}

const Cycles& CyclesGrowthManager::get_cycles(
    bool throw_if_cycles_are_not_candidates) const {
  // GCOVR_EXCL_START
  TKET_ASSERT(
      !(throw_if_cycles_are_not_candidates && !m_cycles_are_candidates));
  // GCOVR_EXCL_STOP
  return m_cycles;
}

bool CyclesGrowthManager::reset(
    const VertexMapping& vertex_mapping, DistancesInterface& distances,
    NeighboursInterface& neighbours) {
  m_cycles.clear();
  m_cycles_are_candidates = false;

  // OK, a bit inefficient, every really good swap (decreasing L by 2)
  // will appear twice, but not a disaster.
  // If no such swap exists, then this stored data will be necessary, because
  // direction matters in longer cycles: v0->v1->v2->v0 is very different
  // from v2->v1->v0->v2.
  // It's simplest just to treat swaps as special cases of cycles, on 2
  // vertices.
  for (const auto& entry : vertex_mapping) {
    const auto source = entry.first;
    const auto target = entry.second;
    const auto source_distance_to_target = distances(source, target);
    if (source_distance_to_target == 0) {
      continue;
    }
    const auto& adj_vertices = neighbours(source);
    for (auto adj_v : adj_vertices) {
      const auto other_v_distance_to_target = distances(adj_v, target);
      if (other_v_distance_to_target < source_distance_to_target) {
        const auto new_id = m_cycles.emplace_back();
        auto& cycle = m_cycles.at(new_id);
        cycle.decrease = 1;
        cycle.vertices.resize(2);
        cycle.vertices[0] = source;
        cycle.vertices[1] = adj_v;
        if (m_cycles.size() >= m_options.max_number_of_cycles) {
          return true;
        }
      }
    }
  }
  return !m_cycles.empty();
}

bool CyclesGrowthManager::attempt_to_close_cycles(
    const VertexMapping& vertex_mapping, DistancesInterface& distances) {
  TKET_ASSERT(!m_cycles_are_candidates);
  for (auto id_opt = m_cycles.front_id(); id_opt;) {
    const auto id = id_opt.value();
    id_opt = m_cycles.next(id);
    auto& cycle = m_cycles.at(id);
    const int decrease = get_move_decrease(
        vertex_mapping, cycle.vertices.back(), cycle.vertices[0], distances);
    const int new_decrease = cycle.decrease + decrease;
    if (new_decrease > 0) {
      cycle.decrease = new_decrease;
      if (!m_cycles_are_candidates) {
        // It's the first good one, so delete all previous.
        for (auto prev_id_opt = m_cycles.previous(id); prev_id_opt;) {
          const auto id_to_be_deleted = prev_id_opt.value();
          prev_id_opt = m_cycles.previous(id_to_be_deleted);
          m_cycles.erase(id_to_be_deleted);
        }
      }
      m_cycles_are_candidates = true;
    } else {
      // Not a good closed cycle; do we delete it?
      if (m_cycles_are_candidates) {
        m_cycles.erase(id);
      }
    }
  }
  return m_cycles_are_candidates;
}

CyclesGrowthManager::GrowthResult CyclesGrowthManager::attempt_to_grow(
    const VertexMapping& vertex_mapping, DistancesInterface& distances,
    NeighboursInterface& neighbours) {
  GrowthResult result;

  TKET_ASSERT(!m_cycles.empty());

  if (m_cycles.front().vertices.size() >= m_options.max_cycle_size) {
    m_cycles.clear();
    result.hit_cycle_length_limit = true;
    result.empty = true;
    return result;
  }
  for (auto id_opt = m_cycles.front_id(); id_opt;) {
    const auto id = id_opt.value();
    id_opt = m_cycles.next(id);

    // Add an arrow onto the back.
    const auto back_vertex = m_cycles.at(id).vertices.back();
    const auto& adj_vertices = neighbours(back_vertex);
    for (auto adj_v : adj_vertices) {
      int new_decr;
      {
        // Important not to reuse this once cycles are added,
        // as it may be invalidated
        const auto& cycle = m_cycles.at(id);
        if (cycle.contains(adj_v)) {
          continue;
        }
        new_decr =
            cycle.decrease +
            get_move_decrease(vertex_mapping, back_vertex, adj_v, distances);

        // If there are N moves, each move can only decrease L by at most one,
        // so it's unfair to demand a huge L-decrease, because shorter cycles
        // would be killed immediately.
        // With N vertices there are N-1 moves, but we are about to add
        // the new vertex adj_v to this partial cycle (unless we discard it),
        // taking it back up to N.
        const int num_moves = cycle.vertices.size();
        int min_decrease = num_moves;
        min_decrease =
            std::min(min_decrease, m_options.min_decrease_for_partial_path);

        // We want 100*(L-decr)/(num.moves) >=
        // min_power_percentage_for_partial_path. But we need the ceiling
        // because of interger division.
        min_decrease = std::max(
            min_decrease,
            (99 + m_options.min_power_percentage_for_partial_path * num_moves) /
                100);

        if (new_decr < min_decrease) {
          continue;
        }
      }
      // A new cycle to be added. Add it before the current position,
      // so we won't pass through it again in the main loop.
      const auto new_id = m_cycles.insert_before(id);
      auto& new_cycle = m_cycles.at(new_id);
      new_cycle.decrease = new_decr;
      new_cycle.vertices = m_cycles.at(id).vertices;
      new_cycle.vertices.push_back(adj_v);
      if (m_cycles.size() >= m_options.max_number_of_cycles) {
        // Break out of the INNER loop, i.e. neighbours for this
        // cycle endpoint. However, this cycle is about to be deleted,
        // creating space, so continue with further cycles.
        break;
      }
    }
    m_cycles.erase(id);
  }
  result.empty = m_cycles.empty();
  return result;
}

}  // namespace tsa_internal
}  // namespace tket
