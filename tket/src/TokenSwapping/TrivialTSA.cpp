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

#include "TrivialTSA.hpp"

#include <sstream>
#include <stdexcept>

#include "CyclicShiftCostEstimate.hpp"
#include "TokenSwapping/DistanceFunctions.hpp"
#include "TokenSwapping/GeneralFunctions.hpp"
#include "TokenSwapping/VertexSwapResult.hpp"
#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

// Make an arrow from each nonempty vertex to its target;
// what are the connected components of the resulting directed graph?
// Two different arrows cannot point INTO the same vertex.
// So, EITHER a cycle (so, a, abstract cyclic shift on tokens is performed),
// OR a path, with all except the final vertex being nonempty.
// In either case, we enact a cyclic shift.

// To find a component, we might have to go backwards along arrows
// as well as forwards.

TrivialTSA::TrivialTSA(Options options) : m_options(options) {
  m_name = "Trivial";
}

void TrivialTSA::set(Options options) { m_options = options; }

bool TrivialTSA::grow_cycle_forwards(
    const VertexMapping& vertex_mapping, Endpoints& endpoints) {
  auto current_id = endpoints.first;
  const auto start_vertex = m_abstract_cycles_vertices.at(current_id);

  // If valid, a single cycle contains at most one empty vertex.
  // Thus there are at most N+1 vertices.
  for (size_t infin_loop_guard = vertex_mapping.size() + 1;
       infin_loop_guard != 0; --infin_loop_guard) {
    const auto v1 = m_abstract_cycles_vertices.at(current_id);
    const auto citer = vertex_mapping.find(v1);
    if (citer == vertex_mapping.cend()) {
      // We end at an empty vertex.
      endpoints.second = current_id;
      return false;
    }
    if (citer->second == start_vertex) {
      // We've hit the start.
      endpoints.second = current_id;
      return true;
    }
    current_id = m_abstract_cycles_vertices.insert_after(current_id);
    m_abstract_cycles_vertices.at(current_id) = citer->second;
  }
  TKET_ASSERT(!"TrivialTSA::grow_cycle_forwards: "
      "hit vertex count limit; invalid vertex mapping");
  return false;
}

void TrivialTSA::grow_cycle_backwards(Endpoints& endpoints) {
  auto current_id = endpoints.first;

  // In a valid cycle, every vertex but one (the empty vertex)
  // is the target of something, and therefore there are <= N+1 vertices.
  for (size_t infin_loop_guard = m_reversed_vertex_mapping.size() + 1;
       infin_loop_guard != 0; --infin_loop_guard) {
    const auto v1 = m_abstract_cycles_vertices.at(current_id);
    const auto citer = m_reversed_vertex_mapping.find(v1);
    if (citer == m_reversed_vertex_mapping.cend()) {
      // Our vertex is not the target of anything.
      // So, it's the START.
      endpoints.first = current_id;
      return;
    }
    // Remember the reverse order!
    current_id = m_abstract_cycles_vertices.insert_before(current_id);
    m_abstract_cycles_vertices.at(current_id) = citer->second;
  }
  TKET_ASSERT(!"TrivialTSA::grow_cycle_backwards: "
      "hit vertex count limit; invalid vertex mapping");
}

void TrivialTSA::do_final_checks() const {
  m_vertices_seen.clear();
  for (const auto& entry : m_reversed_vertex_mapping) {
    m_vertices_seen.insert(entry.first);
    m_vertices_seen.insert(entry.second);
  }
  TKET_ASSERT(m_vertices_seen.size() == m_abstract_cycles_vertices.size());

  // Erase them again...!
  for (const auto& endpoints : m_cycle_endpoints) {
    for (auto id = endpoints.first;;
         id = m_abstract_cycles_vertices.next(id).value()) {
      // GCOVR_EXCL_START
      TKET_ASSERT(
          m_vertices_seen.erase(m_abstract_cycles_vertices.at(id)) == 1);
      // GCOVR_EXCL_STOP
      if (id == endpoints.second) {
        break;
      }
    }
  }
  TKET_ASSERT(m_vertices_seen.empty());
}

void TrivialTSA::fill_disjoint_abstract_cycles(
    const VertexMapping& vertex_mapping) {
  m_vertices_seen.clear();
  m_abstract_cycles_vertices.clear();
  m_cycle_endpoints.clear();
  m_reversed_vertex_mapping = get_reversed_map(vertex_mapping);
  Endpoints endpoints;

  // Get the disjoint abstract cycles.
  for (const auto& entry : vertex_mapping) {
    if (m_vertices_seen.count(entry.first) != 0) {
      continue;
    }
    m_abstract_cycles_vertices.push_back(entry.first);
    endpoints.first = m_abstract_cycles_vertices.back_id().value();
    if (!grow_cycle_forwards(vertex_mapping, endpoints)) {
      grow_cycle_backwards(endpoints);
    }
    m_cycle_endpoints.push_back(endpoints);

    // Now, add the vertices to vertices seen...
    for (auto id = endpoints.first;;
         id = m_abstract_cycles_vertices.next(id).value()) {
      // GCOVR_EXCL_START
      TKET_ASSERT(
          m_vertices_seen.insert(m_abstract_cycles_vertices.at(id)).second);
      // GCOVR_EXCL_STOP
      if (id == endpoints.second) {
        break;
      }
    }
  }
}

void TrivialTSA::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, NeighboursInterface& /*not needed*/,
    RiverFlowPathFinder& path_finder) {
  append_partial_solution(swaps, vertex_mapping, distances, path_finder);
}

void TrivialTSA::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, RiverFlowPathFinder& path_finder) {
  if (all_tokens_home(vertex_mapping)) {
    return;
  }
  fill_disjoint_abstract_cycles(vertex_mapping);
  do_final_checks();

  if (m_options == Options::FULL_TSA) {
    // OK, below, for a single cycle, we use CyclicShiftCostEstimate
    // to estimate, not ONLY the cheapest single cycle, but ALSO
    // the start vertex to enact it most cheaply.
    // We could do that here also and it might save a bit,
    // BUT the full Trivial TSA is really only used for testing now
    // so don't bother.
    append_partial_solution_with_all_cycles(swaps, vertex_mapping, path_finder);
    return;
  }
  TKET_ASSERT(m_options == Options::BREAK_AFTER_PROGRESS);
  // We're only going to do ONE cycle; so find which cycle
  // has the shortest estimated number of swaps
  size_t best_estimated_concrete_swaps = std::numeric_limits<size_t>::max();
  Endpoints best_endpoints;
  size_t start_v_index = std::numeric_limits<size_t>::max();

  for (const auto& endpoints : m_cycle_endpoints) {
    copy_vertices_to_work_vector(endpoints);
    if (m_vertices_work_vector.size() < 2) {
      TKET_ASSERT(m_vertices_work_vector.size() == 1);
      continue;
    }
    const CyclicShiftCostEstimate estimate(m_vertices_work_vector, distances);
    // GCOVR_EXCL_START
    TKET_ASSERT(
        estimate.estimated_concrete_swaps < std::numeric_limits<size_t>::max());
    TKET_ASSERT(estimate.start_v_index < m_vertices_work_vector.size());
    // GCOVR_EXCL_STOP
    if (estimate.estimated_concrete_swaps < best_estimated_concrete_swaps) {
      best_estimated_concrete_swaps = estimate.estimated_concrete_swaps;
      start_v_index = estimate.start_v_index;
      best_endpoints = endpoints;
    }
  }
  // GCOVR_EXCL_START
  TKET_ASSERT(
      best_estimated_concrete_swaps < std::numeric_limits<size_t>::max());
  // GCOVR_EXCL_STOP
  const auto swap_size_before = swaps.size();
  const auto decrease = append_partial_solution_with_single_cycle(
      best_endpoints, start_v_index, swaps, vertex_mapping, distances,
      path_finder);
  TKET_ASSERT(swap_size_before < swaps.size());
  TKET_ASSERT(decrease > 0);
}

void TrivialTSA::copy_vertices_to_work_vector(const Endpoints& endpoints) {
  m_vertices_work_vector.clear();
  for (auto id = endpoints.first;;
       id = m_abstract_cycles_vertices.next(id).value()) {
    m_vertices_work_vector.push_back(m_abstract_cycles_vertices.at(id));
    if (id == endpoints.second) {
      break;
    }
  }
}

void TrivialTSA::append_partial_solution_with_all_cycles(
    SwapList& swaps, VertexMapping& vertex_mapping,
    RiverFlowPathFinder& path_finder) {
  for (const auto& endpoints : m_cycle_endpoints) {
    copy_vertices_to_work_vector(endpoints);
    if (m_vertices_work_vector.size() < 2) {
      continue;
    }
    // Break the abstract cycle into abstract swaps...
    // To shift: [a,b,c,d] -> [d,a,b,c], we do abstract swaps in
    // opposite order of the shift direction, i.e.   cd bc ab
    for (size_t ii = m_vertices_work_vector.size() - 1; ii > 0; --ii) {
      // Abstract swap(v1, v2).
      const auto v1 = m_vertices_work_vector[ii];
      const auto v2 = m_vertices_work_vector[ii - 1];
      TKET_ASSERT(v1 != v2);
      const auto& path = path_finder(v1, v2);
      TKET_ASSERT(path.size() >= 2);
      append_swaps_to_interchange_path_ends(path, vertex_mapping, swaps);
    }
  }
}

size_t TrivialTSA::append_partial_solution_with_single_cycle(
    const Endpoints& endpoints, size_t start_v_index, SwapList& swaps,
    VertexMapping& vertex_mapping, DistancesInterface& distances,
    RiverFlowPathFinder& path_finder) {
  copy_vertices_to_work_vector(endpoints);
  TKET_ASSERT(m_vertices_work_vector.size() >= 2);
  TKET_ASSERT(start_v_index < m_vertices_work_vector.size());

  // Can go negative! But MUST be >= 1 at the end
  // (otherwise this cycle was useless and should never have occurred).
  int current_L_decrease = 0;

  // To shift: [a,b,c,d] -> [d,a,b,c], we do abstract swaps in the opposite
  // order to the shift direction, i.e.   cd bc ab
  for (size_t ii = m_vertices_work_vector.size() - 1; ii > 0; --ii) {
    // Abstract swap(v1, v2).
    const auto v1 = m_vertices_work_vector
        [(ii + start_v_index) % m_vertices_work_vector.size()];

    const auto v2 = m_vertices_work_vector
        [((ii - 1) + start_v_index) % m_vertices_work_vector.size()];

    TKET_ASSERT(v1 != v2);
    const auto& path = path_finder(v1, v2);
    TKET_ASSERT(path.size() >= 2);

    // e.g., to swap endpoints:  [x,a,b,c,y] -> [y,a,b,c,x],
    // do concrete swaps xa ab bc cy bc ab xa.

    //  xa ab bc cy ...(ascending)
    for (size_t jj = 1; jj < path.size(); ++jj) {
      current_L_decrease +=
          get_swap_decrease(vertex_mapping, path[jj], path[jj - 1], distances);

      VertexSwapResult(path[jj], path[jj - 1], vertex_mapping, swaps);
      if (current_L_decrease > 0) {
        return static_cast<size_t>(current_L_decrease);
      }
    }
    // Now the reverse: bc ab xa
    for (size_t kk = path.size() - 2; kk > 0; --kk) {
      current_L_decrease +=
          get_swap_decrease(vertex_mapping, path[kk], path[kk - 1], distances);

      VertexSwapResult(path[kk], path[kk - 1], vertex_mapping, swaps);
      if (current_L_decrease > 0) {
        return static_cast<size_t>(current_L_decrease);
      }
    }
  }
  // The cycle MUST have decreased L overall,
  // otherwise we shouldn't have done it.
  TKET_ASSERT(!"TrivialTSA::append_partial_solution_with_single_cycle");
  return 0;
}

}  // namespace tsa_internal
}  // namespace tket
