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

#include "BestTsaTester.hpp"

#include <catch2/catch_test_macros.hpp>

#include "Architecture/BestTsaWithArch.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"
#include "TokenSwapping/VertexSwapResult.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {

// We are going to treat the raw data in FixedSwapSequences etc. as
// the "correct" data, which we don't want to relabel or process further.
//
// But when an Architecture object is created with a vector of edges,
// given by <unsigned, unsigned> pairs ("raw" vertices),
// vertex relabelling takes place.
// Thus we need an extra layer of conversion to get back what we want.
struct VertexRelabellingManager {
  std::map<unsigned, size_t> raw_to_internal_map;
  // The internal indices are, of course, 0,1,2,...,N for some N,
  // and therefore we can use a vector instead of a map.
  vector<unsigned> internal_to_raw_map;

  // The exact same edges that were used to construct the Architecture object
  // (in the same order!) must be passed in.
  explicit VertexRelabellingManager(
      const vector<std::pair<unsigned, unsigned>>& raw_edges) {
    for (auto edge : raw_edges) {
      size_t next_index = raw_to_internal_map.size();
      if (raw_to_internal_map.count(edge.first) == 0) {
        raw_to_internal_map[edge.first] = next_index;
      }
      next_index = raw_to_internal_map.size();
      if (raw_to_internal_map.count(edge.second) == 0) {
        raw_to_internal_map[edge.second] = next_index;
      }
    }
    internal_to_raw_map.resize(raw_to_internal_map.size());
    for (const auto& entry : raw_to_internal_map) {
      internal_to_raw_map[entry.second] = entry.first;
    }
  }
  Swap get_raw_swap(Swap internal_swap) const {
    return get_swap(
        internal_to_raw_map.at(internal_swap.first),
        internal_to_raw_map.at(internal_swap.second));
  }

  // To be used as input to the TSA.
  // Gives the source->target mappings for INTERNAL vertices.
  VertexMapping get_internal_mapping_for_tsa_input(
      const VertexMapping& raw_mapping) const {
    VertexMapping mapping;
    for (const auto& entry : raw_mapping) {
      mapping[raw_to_internal_map.at(entry.first)] =
          raw_to_internal_map.at(entry.second);
    }
    return mapping;
  }
};
}  // namespace

size_t BestTsaTester::get_checked_solution_size(
    const DecodedProblemData& problem_data) {
  m_architecture_work_data.edges.clear();
  for (const auto& swap : problem_data.swaps) {
    m_architecture_work_data.edges.insert(swap);
  }
  m_architecture_work_data.number_of_vertices = 0;
  return get_checked_solution_size(problem_data, m_architecture_work_data);
}

size_t BestTsaTester::get_checked_solution_size(
    const DecodedProblemData& problem_data,
    const DecodedArchitectureData& architecture_data) {
  CHECK(problem_data.number_of_vertices >= 4);
  if (architecture_data.number_of_vertices > 0) {
    CHECK(
        architecture_data.number_of_vertices >=
        problem_data.number_of_vertices);
  }
  // problem_data.number_of_vertices only includes the vertices mentioned in the
  // solution swaps.
  // architecture_data.number_of_vertices is EITHER set to zero,
  // OR is calculated from the EDGES in the architecture, and hence is correct.
  const auto number_of_vertices = std::max(
      architecture_data.number_of_vertices, problem_data.number_of_vertices);

  check_mapping(problem_data.vertex_mapping);
  for (const auto& swap : problem_data.swaps) {
    REQUIRE(architecture_data.edges.count(swap) != 0);
  }
  for (const auto& edge : architecture_data.edges) {
    REQUIRE(edge.first < number_of_vertices);
    REQUIRE(edge.second < number_of_vertices);
  }
  m_edges_vect = vector<std::pair<unsigned, unsigned>>{
      architecture_data.edges.cbegin(), architecture_data.edges.cend()};

  REQUIRE(problem_data.vertex_mapping.size() >= 1);
  REQUIRE(problem_data.vertex_mapping.size() <= number_of_vertices);
  REQUIRE(problem_data.vertex_mapping.crbegin()->first < number_of_vertices);

  const bool full_tokens =
      problem_data.vertex_mapping.size() == number_of_vertices;

  const Architecture arch(m_edges_vect);
  const ArchitectureMapping arch_mapping(arch, m_edges_vect);
  const VertexRelabellingManager relabelling_manager(m_edges_vect);
  m_raw_swap_list.clear();
  m_vertex_mapping_copy =
      relabelling_manager.get_internal_mapping_for_tsa_input(
          problem_data.vertex_mapping);

  BestTsaWithArch::append_solution(
      m_raw_swap_list, m_vertex_mapping_copy, arch_mapping);

  // Now check the calculated solution.
  // Set it back to the raw, i.e. "proper" mapping.
  m_vertex_mapping_copy = problem_data.vertex_mapping;

  for (auto id_opt = m_raw_swap_list.front_id(); id_opt;) {
    const auto id = id_opt.value();
    id_opt = m_raw_swap_list.next(id);
    auto& swap = m_raw_swap_list.at(id);
    // This is an "internal" swap, so needs conversion back to "raw".
    swap = relabelling_manager.get_raw_swap(swap);

    const VertexSwapResult vswap_result(swap, m_vertex_mapping_copy);
    if (full_tokens) {
      REQUIRE(vswap_result.tokens_moved == 2);
    } else {
      // We require our best TSA to avoid empty swaps.
      REQUIRE(vswap_result.tokens_moved >= 1);
      REQUIRE(vswap_result.tokens_moved <= 2);
    }
    REQUIRE(architecture_data.edges.count(swap) != 0);
  }
  REQUIRE(all_tokens_home(m_vertex_mapping_copy));
  return m_raw_swap_list.size();
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
