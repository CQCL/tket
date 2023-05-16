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

#include "PlacementCostModelInterface.hpp"

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <iostream>
#include <sstream>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

typedef std::chrono::steady_clock Clock;

PlacementCostModelInterface::~PlacementCostModelInterface() {}

void PlacementCostModelInterface::initialise_with_qubit_placement(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& placement) {
  REQUIRE(placement.size() > 2);
  m_current_placement.clear();
  m_current_tokens.clear();
  for (const auto& pair : placement) {
    m_current_placement[pair.first] = pair.second;
    m_current_tokens[pair.second] = pair.first;
  }
  REQUIRE(m_current_placement.size() == placement.size());
  REQUIRE(m_current_tokens.size() == placement.size());
  REQUIRE(placement.crbegin()->first < m_invalid_token / 2);
  require_valid();
}

WeightWSM PlacementCostModelInterface::get_cost(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates) {
  WeightWSM cost = 0;
  for (const auto& gate : gates) {
    cost += do_token_swapping_and_apply_gate(gate.first, gate.second);
  }
  return cost;
}

WeightWSM PlacementCostModelInterface::do_token_swapping_and_apply_gate(
    VertexWSM pv1, VertexWSM pv2) {
  require_valid();
  REQUIRE(pv1 < m_invalid_token);
  REQUIRE(pv2 < m_invalid_token);
  REQUIRE(pv1 != pv2);
  const auto vertex1 = m_current_placement.at(pv1);
  const auto vertex2 = m_current_placement.at(pv2);
  const Path& path = get_path_to_use(vertex1, vertex2);
  REQUIRE(path.size() >= 2);
  REQUIRE(path[0].first == vertex1);
  REQUIRE(path[0].second == 0);
  REQUIRE(path.back().first == vertex2);

  const WeightsData weights_data(path);
  enact_swaps(path, weights_data);
  return weights_data.highest_edge_weight +
         weights_data.sum_of_swap_edge_weights *
             m_number_of_primitive_gates_in_swap;
}

PlacementCostModelInterface::WeightsData::WeightsData(
    const PlacementCostModelInterface::Path& path)
    : sum_of_swap_edge_weights(0),
      highest_edge_weight(0),
      index_of_largest_edge_weight(1) {
  REQUIRE(path.size() >= 2);

  // Dummy weight.
  REQUIRE(path[0].second == 0);

  for (unsigned ii = 1; ii < path.size(); ++ii) {
    sum_of_swap_edge_weights += path[ii].second;
    if (path[ii].second > path[index_of_largest_edge_weight].second) {
      index_of_largest_edge_weight = ii;
    }
  }
  highest_edge_weight = path[index_of_largest_edge_weight].second;
  REQUIRE(sum_of_swap_edge_weights >= highest_edge_weight);
  sum_of_swap_edge_weights -= highest_edge_weight;
}

void PlacementCostModelInterface::enact_swaps(
    const Path& path, const WeightsData& weights_data) {
  require_valid();

  // On both pieces of the path, we do a cyclic shift.
  //
  // E.g., consider   index_of_largest_edge_weight=3,
  //      path vertices [v0, v1, v2, v3, v4, v5, v6],
  // with current tokens [t0, t1, t2, t3, t4, t5, t6],
  // so that Weight(v2, v3) is the largest.
  // (Since W[i] = Weight(v[i-1], v[i]) by definition).
  //

  // We will cyclically shift the head and tail, so that the new tokens
  // at those same vertices afterwards are
  //    [t1, t2, t0, t6, t3, t4, t5].
  // (so that the end tokens t0, t6 now sit at the vertices v2, v3
  // forming the largest edge weight).
  const VertexWSM start_v = path[0].first;
  const VertexWSM destination_v =
      path[weights_data.index_of_largest_edge_weight - 1].first;

  // First, ensure that all token data exists.
  // We cannot be swapping dummy tokens.
  REQUIRE(m_current_tokens.at(path[0].first) < m_invalid_token);
  REQUIRE(m_current_tokens.at(path.back().first) < m_invalid_token);

  for (const auto& entry : path) {
    const auto token_opt = get_optional_value(m_current_tokens, entry.first);
    if (!token_opt) {
      // For debugging, it's handy to have distinct token values.
      VertexWSM dummy_token = m_current_tokens.crbegin()->second;
      if (dummy_token < m_invalid_token) {
        dummy_token += m_invalid_token;
      }
      m_current_tokens[entry.first] = dummy_token;
    }
  }

  if (start_v != destination_v) {
    const VertexWSM start_token = m_current_tokens.at(start_v);
    for (unsigned ii = 0; ii + 1 < weights_data.index_of_largest_edge_weight;
         ++ii) {
      const VertexWSM new_token_at_v = m_current_tokens.at(path[ii + 1].first);
      m_current_tokens.at(path[ii].first) = new_token_at_v;
      if (new_token_at_v < m_invalid_token) {
        m_current_placement.at(new_token_at_v) = path[ii].first;
      }
    }
    m_current_tokens.at(destination_v) = start_token;
    m_current_placement.at(start_token) = destination_v;
  }

  const VertexWSM tail_start_v = path.back().first;
  const VertexWSM tail_destination_v =
      path[weights_data.index_of_largest_edge_weight].first;
  if (tail_start_v != tail_destination_v) {
    const VertexWSM tail_start_token = m_current_tokens.at(tail_start_v);
    for (unsigned ii = path.size() - 1;
         ii >= weights_data.index_of_largest_edge_weight + 1; --ii) {
      const VertexWSM new_token_at_v = m_current_tokens.at(path[ii - 1].first);
      m_current_tokens.at(path[ii].first) = new_token_at_v;
      if (new_token_at_v < m_invalid_token) {
        m_current_placement.at(new_token_at_v) = path[ii].first;
      }
    }
    m_current_tokens.at(tail_destination_v) = tail_start_token;
    m_current_placement.at(tail_start_token) = tail_destination_v;
  }
  require_valid();
}

const std::map<VertexWSM, VertexWSM>&
PlacementCostModelInterface::get_current_placement() const {
  return m_current_placement;
}

const std::map<VertexWSM, VertexWSM>&
PlacementCostModelInterface::get_current_tokens() const {
  return m_current_tokens;
}

void PlacementCostModelInterface::require_valid() const {
  for (const auto& entry : m_current_placement) {
    const auto& token = entry.first;
    const auto& vertex = entry.second;
    REQUIRE(token < m_invalid_token);
    REQUIRE(m_current_tokens.at(vertex) == token);
  }
  for (const auto& entry : m_current_tokens) {
    const auto& vertex = entry.first;
    const auto& token = entry.second;
    if (token < m_invalid_token) {
      REQUIRE(m_current_placement.at(token) == vertex);
    }
  }
}

void PlacementCostModelInterface::try_random_placements(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gates) {
  auto vertices_vector = get_vertices(get_graph_data());
  std::set<VertexWSM> gate_pv_used;
  for (const auto& entry : gates) {
    gate_pv_used.insert(entry.first);
    gate_pv_used.insert(entry.second);
  }
  unsigned iterations_report = 100;
  WeightWSM best_cost;
  set_maximum(best_cost);
  RNG rng;
  std::vector<EdgeWSM> placement_vector;
  WeightWSM current_cost;
  set_maximum(current_cost);

  const auto start = Clock::now();

  // This is, of course, a simplified version of the Monte Carlo algorithm
  // for complete target graphs - easy because we don't need to terminate!!

  const unsigned max_iters_without_progress = 10 * vertices_vector.size();
  const unsigned max_iters_without_record_breaker =
      100 * vertices_vector.size();
  unsigned next_reset_iter_if_no_progress = 0;
  unsigned next_reset_iter_if_no_record_breaker = 0;

  for (unsigned iterations = 0;; ++iterations) {
    if (iterations >= next_reset_iter_if_no_progress ||
        iterations >= next_reset_iter_if_no_record_breaker) {
      // We'll reset, try a new solution.
      set_maximum(current_cost);
      rng.do_shuffle(vertices_vector);
      next_reset_iter_if_no_progress = iterations + max_iters_without_progress;
      next_reset_iter_if_no_record_breaker =
          iterations + max_iters_without_record_breaker;
      continue;
    }
    const auto v_index1 = rng.get_size_t(vertices_vector.size() - 1);
    const auto v_index2 = rng.get_size_t(vertices_vector.size() - 1);
    if (v_index1 == v_index2) {
      continue;
    }
    std::swap(vertices_vector[v_index1], vertices_vector[v_index2]);
    placement_vector.clear();
    {
      unsigned next_ii = 0;
      for (auto pv : gate_pv_used) {
        placement_vector.emplace_back(pv, vertices_vector[next_ii]);
        ++next_ii;
      }
    }
    initialise_with_qubit_placement(placement_vector);
    WeightWSM cost = get_cost(gates);
    if (cost < best_cost) {
      std::cerr << "\nIter=" << iterations << "; after "
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       Clock::now() - start)
                       .count()
                << " ms, found new best cost " << cost << " for placement { ";
      for (const auto& entry : placement_vector) {
        std::cerr << "{" << entry.first << "," << entry.second << "}, ";
      }
      std::cerr << "}\n";
      best_cost = cost;
      current_cost = cost;
      next_reset_iter_if_no_progress = iterations + max_iters_without_progress;
      next_reset_iter_if_no_record_breaker =
          iterations + max_iters_without_record_breaker;
      continue;
    }
    if (cost <= current_cost) {
      // Accept the step.
      next_reset_iter_if_no_progress = iterations + max_iters_without_progress;
    } else {
      // Reject!
      std::swap(vertices_vector[v_index1], vertices_vector[v_index2]);
    }
    // Print periodic updates, so the caller can see something is happening...
    if (iterations >= iterations_report) {
      std::cerr << "\ni=" << iterations;
      iterations_report = iterations * 2;
    }
  }
}

void PlacementCostModelInterface::reverse_path(Path& path) {
  std::reverse(path.begin(), path.end());

  // But remember, the path[i] weight is the weight FROM v[i-1] TO v[i],
  // so we must now do a cyclic shift on the weights.
  //
  // E.g., if we std::reverse  [ (v0, 0)  (v1, w01)  (v2, w12)  (v3, w23)]
  // we get   [ (v3, w23)  (v2, w12)  (v1, w01)  (v0, 0) ].
  // But we actually want
  //
  //  [ (v3, 0)  (v2, w23)  (v1, w12)  (v0, w01) ].
  //
  REQUIRE(path.back().second == 0);
  for (unsigned ii = path.size() - 1; ii > 0; --ii) {
    path[ii].second = path[ii - 1].second;
  }
  path[0].second = 0;
}

WeightWSM PlacementCostModelInterface::get_total_weight(const Path& path) {
  WeightWSM weight = 0;
  for (const auto& entry : path) {
    weight += entry.second;
  }
  return weight;
}

std::string PlacementCostModelInterface::get_path_str(
    VertexWSM v1, VertexWSM v2) const {
  const auto& path = get_path_to_use(v1, v2);
  REQUIRE(path.size() >= 2);
  REQUIRE(path[0].first == v1);
  REQUIRE(path[0].second == 0);
  REQUIRE(path.back().first == v2);

  std::stringstream ss;
  ss << "\nPath vertices: [ ";
  for (const auto& entry : path) {
    ss << entry.first << " ";
  }
  ss << "]\nEdge weights: [ ";
  WeightWSM weight = 0;
  for (unsigned ii = 1; ii < path.size(); ++ii) {
    ss << path[ii].second << " ";
    weight += path[ii].second;
  }
  ss << "] (total weight " << weight << ")";
  REQUIRE(weight == get_total_weight(path));
  return ss.str();
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
