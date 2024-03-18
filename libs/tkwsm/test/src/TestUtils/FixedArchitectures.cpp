// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "FixedArchitectures.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

GraphEdgeWeights FixedArchitectures::get_ibm_brooklyn_65_qubits() {
  auto result = get_path_qubits(0, 9);
  merge_first_with_second(result, get_path_qubits(13, 23));
  merge_first_with_second(result, get_path_qubits(27, 37));
  merge_first_with_second(result, get_path_qubits(41, 51));
  merge_first_with_second(result, get_path_qubits(55, 64));

  add_paths(
      result, {{0, 10, 13},
               {4, 11, 17},
               {8, 12, 21},
               {15, 24, 29},
               {19, 25, 33},
               {23, 26, 37},
               {27, 38, 41},
               {31, 39, 45},
               {35, 40, 49},
               {43, 52, 56},
               {47, 53, 60},
               {51, 54, 64}});
  return result;
}

GraphEdgeWeights FixedArchitectures::get_ibm_montreal_27_qubits() {
  auto result = get_path_qubits({0,  1,  4,  7,  10, 12, 15, 18, 21, 23, 24,
                                 25, 22, 19, 16, 14, 11, 8,  5,  3,  2,  1});
  const std::vector<EdgeWSM> extra_edges{{6, 7}, {17, 18}, {12, 13}, {13, 14},
                                         {8, 9}, {19, 20}, {25, 26}};
  add_edges(result, extra_edges);
  return result;
}

GraphEdgeWeights FixedArchitectures::get_ibm_guadalupe_16_qubits() {
  auto result =
      get_path_qubits({0, 1, 4, 7, 10, 12, 13, 14, 11, 8, 5, 3, 2, 1});
  const std::vector<EdgeWSM> extra_edges{{6, 7}, {12, 15}, {8, 9}};
  add_edges(result, extra_edges);
  return result;
}

GraphEdgeWeights FixedArchitectures::get_ibm_perth_7_qubits() {
  auto result = get_path_qubits({0, 1, 3, 5, 6});
  const std::vector<EdgeWSM> extra_edges{{1, 2}, {4, 5}};
  add_edges(result, extra_edges);
  return result;
}

GraphEdgeWeights FixedArchitectures::get_path_qubits(
    const std::vector<VertexWSM>& vertices, bool allow_cycles) {
  if (!allow_cycles) {
    auto vertices_copy = vertices;
    std::sort(vertices_copy.begin(), vertices_copy.end());
    REQUIRE(is_sorted_and_unique(vertices_copy));
  }
  GraphEdgeWeights result;
  for (unsigned ii = 1; ii < vertices.size(); ++ii) {
    const auto& v1 = vertices[ii - 1];
    const auto& v2 = vertices[ii];
    REQUIRE(v1 != v2);
    result[get_edge(v1, v2)] = 1;
  }
  return result;
}

GraphEdgeWeights FixedArchitectures::get_path_qubits(
    VertexWSM first, VertexWSM last) {
  REQUIRE(first != last);
  int step = 1;
  VertexWSM prev_v = first;
  VertexWSM final_v = last;
  if (first > last) {
    prev_v = last;
    final_v = first;
    step = -1;
  }
  GraphEdgeWeights result;
  for (;;) {
    VertexWSM current_v = prev_v + step;
    result[get_edge(prev_v, current_v)] = 1;
    if (current_v == final_v) {
      break;
    }
    prev_v = current_v;
  }
  return result;
}

void FixedArchitectures::add_edges(
    GraphEdgeWeights& data, const std::vector<EdgeWSM>& edges,
    bool require_edges_to_be_new) {
  for (const auto& entry : edges) {
    const auto edge = get_edge(entry.first, entry.second);
    if (require_edges_to_be_new) {
      REQUIRE(data.count(edge) == 0);
    }
    data[edge] = 1;
  }
}

void FixedArchitectures::add_paths(
    GraphEdgeWeights& data, const std::vector<std::vector<VertexWSM>>& paths,
    bool allow_cycles, bool require_edges_to_be_new) {
  for (const auto& path : paths) {
    const auto path_data = get_path_qubits(path, allow_cycles);
    merge_first_with_second(data, path_data, require_edges_to_be_new);
  }
}

void FixedArchitectures::merge_first_with_second(
    GraphEdgeWeights& first, const GraphEdgeWeights& second,
    bool require_edges_to_be_new) {
  for (const auto& entry : second) {
    const auto& edge = entry.first;
    const auto& weight = entry.second;
    if (require_edges_to_be_new) {
      REQUIRE(first.count(edge) == 0);
    }
    first[edge] = weight;
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
