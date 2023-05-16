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

#include "tkwsm/EndToEndWrappers/SolutionWSM.hpp"

#include <sstream>

#include "tkwsm/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// Fills "assignments_map" with the actual assignments made;
// fills "values" with the TV values seen.
static void check_assignments_for_value_clashes(
    Assignments& assignments_map, const SolutionWSM& solution,
    std::set<VertexWSM>& values, std::stringstream& ss) {
  values.clear();
  for (const std::pair<VertexWSM, VertexWSM>& entry : solution.assignments) {
    const auto value_opt = get_optional_value(assignments_map, entry.first);
    if (value_opt) {
      ss << "\nRepeated PV";
    }
    assignments_map[entry.first] = entry.second;

    if (values.count(entry.second) != 0) {
      ss << "\nDuplicate value " << entry.second << " seen, when adding "
         << entry.first << "->" << entry.second;
    }
    values.insert(entry.second);
  }
  if (assignments_map.size() != solution.assignments.size() ||
      (values.size() != solution.assignments.size())) {
    ss << "\nSizes mismatch: " << assignments_map.size() << ","
       << solution.assignments.size() << "," << values.size();
  }
}

std::string SolutionWSM::get_errors(
    const GraphEdgeWeights& pattern_edges_and_weights,
    const GraphEdgeWeights& target_edges_and_weights) const {
  std::stringstream ss;
  if (assignments.empty()) {
    if (scalar_product != 0 || total_p_edges_weight != 0) {
      ss << "empty assignments, but sc.prod=" << scalar_product
         << ", total p.edge.weights=" << total_p_edges_weight;
    }
    return ss.str();
  }
  Assignments assignments_map;
  std::set<VertexWSM> work_set;

  check_assignments_for_value_clashes(assignments_map, *this, work_set, ss);

  // Now, recalculate the weights, checking carefully for overflow.
  WeightWSM total_expected_p_edge_weight = 0;
  WeightWSM expected_scalar_product = 0;
  std::set<VertexWSM>& p_vertices_used = work_set;
  p_vertices_used.clear();

  for (const auto& p_edge_data : pattern_edges_and_weights) {
    const WeightWSM& p_edge_weight = p_edge_data.second;
    total_expected_p_edge_weight =
        get_sum_or_throw(total_expected_p_edge_weight, p_edge_weight);

    const EdgeWSM& p_edge = p_edge_data.first;
    if (p_edge.first == p_edge.second) {
      ss << "\nInvalid loop at PV=" << p_edge.first;
    }
    if (pattern_edges_and_weights.count(
            std::make_pair(p_edge.second, p_edge.first)) != 0) {
      ss << "\nRepeated pattern edge";
    }
    p_vertices_used.insert(p_edge.first);
    p_vertices_used.insert(p_edge.second);
    const auto tv1_opt = get_optional_value(assignments_map, p_edge.first);
    const auto tv2_opt = get_optional_value(assignments_map, p_edge.second);
    if (!tv1_opt || !tv2_opt) {
      ss << "\nP-edge (" << p_edge.first << "," << p_edge.second
         << ") has unassigned vertices";
      break;
    }
    const auto tv1 = tv1_opt.value();
    const auto tv2 = tv2_opt.value();
    if (tv1 == tv2) {
      ss << "\nP vertices " << p_edge.first << "," << p_edge.second
         << " both map to " << tv1;
      break;
    }
    const auto t_edge = get_edge(tv1, tv2);
    const auto t_weight_opt =
        get_optional_value(target_edges_and_weights, t_edge);
    if (!t_weight_opt) {
      ss << "\nP-edge [" << p_edge.first << "," << p_edge.second
         << "] maps to nonexistent target edge [" << tv1 << "," << tv2 << "]";
      break;
    }
    const auto extra_sc_product =
        get_product_or_throw(p_edge_weight, t_weight_opt.value());
    expected_scalar_product =
        get_sum_or_throw(expected_scalar_product, extra_sc_product);
  }
  if (expected_scalar_product != scalar_product ||
      total_expected_p_edge_weight != total_p_edges_weight) {
    ss << "\nWeights mismatch: scalar products " << expected_scalar_product
       << "," << scalar_product << "; total p-edge weights "
       << total_expected_p_edge_weight << "," << total_p_edges_weight;
  }
  if (p_vertices_used.size() != assignments.size()) {
    ss << "\nNumber of used p vertices mismatch: " << p_vertices_used.size()
       << "," << assignments.size();
  }
  return ss.str();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
