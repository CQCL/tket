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

#include "PlacementWithWSM/PlacementAndStatistics.hpp"

#include <algorithm>
#include <sstream>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SpecialExceptions.hpp"
#include "WeightSubgrMono/EndToEndWrappers/SolutionStatistics.hpp"
#include "WeightSubgrMono/Searching/SolutionWSM.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

PlacementAndStatistics::PlacementAndStatistics() {}

PlacementAndStatistics::PlacementAndStatistics(
    const GraphEdgeWeights& pattern_graph,
    const GraphEdgeWeights& original_target_graph,
    const GraphEdgeWeights& enlarged_target_graph,
    const std::vector<std::set<VertexWSM>>& gates,
    const SolutionWSM& solution) {
  // Check whatever solution was returned.
  // TODO: use bipartite matching to match
  // as many vertices as possible!
  std::set<VertexWSM> t_vertices_used;
  for (const auto& entry : solution.assignments) {
    const auto& pv = entry.first;
    const auto& tv = entry.second;
    if (valid_assignments.count(pv) == 0 && t_vertices_used.count(tv) == 0) {
      valid_assignments[pv] = tv;
      t_vertices_used.insert(tv);
    }
  }

  // Now, with these assignments, check how many gates etc. worked out.
  for (const auto& gate : gates) {
    if (gate.size() < 2) {
      ++single_qubit_gates;
      continue;
    }
    if (gate.size() > 2) {
      ++n_many_qubit_gates;
      // Are they all assigned?
      for (auto pv : gate) {
        if (valid_assignments.count(pv) == 0) {
          ++n_many_qubit_gates_unassigned;
          break;
        }
      }
      continue;
    }
    // Now, we have a 2-qubit gate.
    const auto pv1 = *gate.cbegin();
    const auto pv2 = *gate.crbegin();
    const auto p_edge = get_edge(pv1, pv2);
    TKET_ASSERT(pattern_graph.count(p_edge) != 0);
    const auto tv1_opt = get_optional_value(valid_assignments, pv1);
    if (tv1_opt) {
      const auto tv2_opt = get_optional_value(valid_assignments, pv2);
      if (tv2_opt) {
        const auto t_edge = get_edge(tv1_opt.value(), tv2_opt.value());
        if (original_target_graph.count(t_edge) != 0) {
          ++n_gates_with_original_edges;
        } else {
          const auto t_weight_opt =
              get_optional_value(enlarged_target_graph, t_edge);
          if (t_weight_opt) {
            ++n_gates_with_some_token_swapping;
            total_weights_with_token_swapping += t_weight_opt.value();
          } else {
            ++n_poor_gates;
          }
        }
        continue;
      }
    }
    ++n_unassigned_gates;
  }
}

bool PlacementAndStatistics::prefer_other_solution(
    const PlacementAndStatistics& other) const {
  if (valid_assignments.size() < other.valid_assignments.size()) {
    return true;
  }
  if (valid_assignments.size() > other.valid_assignments.size()) {
    return false;
  }
  // Both have the same number of assignments.
  if (n_gates_with_original_edges < other.n_gates_with_original_edges) {
    return true;
  }
  if (n_gates_with_original_edges > other.n_gates_with_original_edges) {
    return false;
  }
  // ...and also the same number of 2-qubit gates nicely assigned.
  if (n_gates_with_some_token_swapping <
      other.n_gates_with_some_token_swapping) {
    return true;
  }
  if (n_gates_with_some_token_swapping >
      other.n_gates_with_some_token_swapping) {
    return false;
  }
  // Higher target weights are worse.
  return total_weights_with_token_swapping >
         other.total_weights_with_token_swapping;
}

std::string PlacementAndStatistics::str() const {
  std::stringstream ss;
  ss << "assigned " << valid_assignments.size() << " qubits; "
     << n_gates_with_original_edges << " twoQ gates in place; "
     << n_gates_with_some_token_swapping << " twoQ gates nearby; "
     << total_weights_with_token_swapping << " total swap weights; "
     << n_poor_gates << " twoQ bad gates; " << n_unassigned_gates
     << " twoQ gates unassigned; " << single_qubit_gates << " oneQ gates; "
     << n_many_qubit_gates << " nQ gates; " << n_many_qubit_gates_unassigned
     << " nQ gates unassigned.";
  return ss.str();
}

}  // namespace tket
