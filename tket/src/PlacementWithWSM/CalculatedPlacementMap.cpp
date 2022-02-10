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

#include "PlacementWithWSM/CalculatedPlacementMap.hpp"

#include <algorithm>

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "PlacementWithWSM/FullPlacementResult.hpp"
#include "PlacementWithWSM/PatternGraphTimeSlices.hpp"
#include "PlacementWithWSM/TargetGraphData.hpp"
#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

// Call this repeatedly to assign the vertex numbers 0,1,2,...
// to the qubits.
// Returns the vertex number of the qubit.
template <typename T>
static VertexWSM add_qubit_to_map(
    std::map<T, VertexWSM>& qubit_to_vertex_map, const T& qubit) {
  const auto vertex_opt = get_optional_value(qubit_to_vertex_map, qubit);
  if (vertex_opt) {
    return vertex_opt.value();
  }
  const VertexWSM new_v = qubit_to_vertex_map.size();
  qubit_to_vertex_map[qubit] = new_v;
  return new_v;
}

namespace {
struct LogicalQubitInformation {
  std::vector<std::set<VertexWSM>> gates_in_order;
  std::map<Qubit, VertexWSM> logical_qubit_to_p_vertex_map;

  explicit LogicalQubitInformation(const Circuit& circ) {
    std::set<VertexWSM> vertex_set;
    const auto commands = circ.get_commands();
    for (const auto& com : commands) {
      const auto logical_qubits = com.get_qubits();
      vertex_set.clear();
      for (const auto& qubit : logical_qubits) {
        vertex_set.insert(
            add_qubit_to_map(logical_qubit_to_p_vertex_map, qubit));
      }
      TKET_ASSERT(logical_qubits.size() == vertex_set.size());
      gates_in_order.emplace_back(std::move(vertex_set));
    }
  }
};

struct PhysicalQubitInformation {
  std::map<Node, VertexWSM> physical_qubit_to_t_vertex_map;
  // Initially, let all edges have weight 1.
  GraphEdgeWeights edges_and_weights;

  explicit PhysicalQubitInformation(const Architecture& arch) {
    unsigned edges_count = 0;
    for (auto [n1, n2] : arch.get_all_edges_vec()) {
      ++edges_count;
      const auto tv1 = add_qubit_to_map(physical_qubit_to_t_vertex_map, n1);
      const auto tv2 = add_qubit_to_map(physical_qubit_to_t_vertex_map, n2);
      edges_and_weights[get_edge(tv1, tv2)] = 1;
    }
    TKET_ASSERT(edges_and_weights.size() == edges_count);
  }
};
}  // namespace

CalculatedPlacementMap::CalculatedPlacementMap(
    const Circuit& circ, const Architecture& arch,
    const Parameters& parameters) {
  const LogicalQubitInformation p_vertex_information(circ);
  const PhysicalQubitInformation t_vertex_information(arch);
  const PatternGraphTimeSlices slices(p_vertex_information.gates_in_order);
  const PatternGraphTimeSlices::WeightParameters p_parameters;
  const auto pattern_graph = slices.get_weights(p_parameters);
  const TargetGraphData::Parameters t_parameters;
  const TargetGraphData target_full_graph(
      t_vertex_information.edges_and_weights, t_parameters);
  FullPlacementResult::Parameters full_result_parameters;
  full_result_parameters.timeout_ms = parameters.timeout_ms;

  full_placement_result = FullPlacementResult(
      pattern_graph, t_vertex_information.edges_and_weights,
      target_full_graph.final_data, p_vertex_information.gates_in_order,
      full_result_parameters);

  const auto pv_to_qubit_map =
      get_reversed_map(p_vertex_information.logical_qubit_to_p_vertex_map);
  const auto tv_to_qubit_map =
      get_reversed_map(t_vertex_information.physical_qubit_to_t_vertex_map);

  for (const auto& vertex_pair :
       full_placement_result.result.valid_assignments) {
    const auto& pv = vertex_pair.first;
    const auto& tv = vertex_pair.second;
    placement_map[pv_to_qubit_map.at(pv)] = tv_to_qubit_map.at(tv);
  }
}

}  // namespace tket
