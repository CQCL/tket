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

#include "WeightSubgrMono/Searching/SolutionWSM.hpp"

#include <sstream>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

static void check_assignments_for_value_clashes(
    Assignments& assignments_map, const SolutionWSM& solution,
    std::set<VertexWSM>& values, std::stringstream& ss) {
  for (const auto& entry : solution.assignments) {
    const auto value_opt = get_optional_value(assignments_map, entry.first);
    if (value_opt) {
      ss << "\nRepeated assignments " << entry.first << "-> " << entry.second
         << " and " << value_opt.value();
    }
    assignments_map[entry.first] = entry.second;

    // Incomplete solutions are allowed to have t-vertex clashes.
    if (values.count(entry.second) != 0 && solution.complete) {
      ss << "\nDuplicate value " << entry.second << " seen, when trying "
         << entry.first << "->" << entry.second;
    }
    values.insert(entry.second);
  }
  if (assignments_map.size() != solution.assignments.size() ||
      (solution.complete && (values.size() != solution.assignments.size()))) {
    ss << "\nSizes mismatch: " << assignments_map.size() << ","
       << solution.assignments.size() << "," << values.size();
  }
}

namespace {
struct WeightChecks {
  std::stringstream& ss;
  WeightWSM recalc_p_edges_weight = 0;
  WeightWSM recalc_total_weight = 0;
  bool total_p_edges_overflow = false;
  bool total_weight_overflow = false;

  explicit WeightChecks(std::stringstream& sstream) : ss(sstream) {}

  void add_p_edge_weights(WeightWSM p_weight) {
    if (total_p_edges_overflow) {
      return;
    }
    const auto sum = get_checked_sum(recalc_p_edges_weight, p_weight);
    if (sum) {
      recalc_p_edges_weight = sum.value();
    } else {
      ss << "\nOverflow calculating total p-weight: " << recalc_p_edges_weight
         << "+" << p_weight;
      total_p_edges_overflow = true;
    }
  }

  void add_scalar_product(WeightWSM p_edge_weight, WeightWSM t_edge_weight) {
    if (total_weight_overflow) {
      return;
    }
    const auto product = get_checked_product(p_edge_weight, t_edge_weight);
    if (!product) {
      ss << "\nOverflow: w(p-edge) * w(t-edge): " << p_edge_weight << "*"
         << t_edge_weight;
      total_weight_overflow = true;
      return;
    }
    const auto prod_value = product.value();
    const auto sum = get_checked_sum(recalc_total_weight, prod_value);
    if (!sum) {
      ss << "\nOverflow calculating total weight: " << recalc_total_weight
         << "+" << prod_value;
      total_weight_overflow = true;
      return;
    }
    recalc_total_weight = sum.value();
  }

  void final_check(const SolutionWSM& solution) {
    if (total_p_edges_overflow || total_weight_overflow ||
        (recalc_p_edges_weight == solution.total_p_edges_weight &&
         recalc_total_weight == solution.total_scalar_product_weight)) {
      return;
    }
    ss << "\nRecalc/orig weights mismatch: p-edges: " << recalc_p_edges_weight
       << "," << solution.total_p_edges_weight << "; scalar product "
       << recalc_total_weight << "," << solution.total_scalar_product_weight;
  }
};

}  // namespace

std::string SolutionWSM::get_errors(
    const GraphEdgeWeights& pattern_edges_and_weights,
    const GraphEdgeWeights& target_edges_and_weights) const {
  std::stringstream ss;
  Assignments assignments_map;
  std::set<VertexWSM> work_set;

  check_assignments_for_value_clashes(assignments_map, *this, work_set, ss);
  if (!complete) {
    return ss.str();
  }

  // Now, recalculate the weights, checking carefully for overflow.

  auto& p_vertices_used = work_set;
  p_vertices_used.clear();

  WeightChecks weight_checks(ss);

  for (const auto& p_edge_data : pattern_edges_and_weights) {
    const auto& p_weight = p_edge_data.second;
    weight_checks.add_p_edge_weights(p_weight);

    const auto& p_edge = p_edge_data.first;
    p_vertices_used.insert(p_edge.first);
    p_vertices_used.insert(p_edge.second);
    const auto tv1_opt = get_optional_value(assignments_map, p_edge.first);
    const auto tv2_opt = get_optional_value(assignments_map, p_edge.second);
    if (!tv1_opt || !tv2_opt) {
      ss << "\nP-edge (" << p_edge.first << "," << p_edge.second
         << ") has unknown vertices";
      continue;
    }
    const auto tv1 = tv1_opt.value();
    const auto tv2 = tv2_opt.value();
    if (tv1 == tv2) {
      ss << "\nP vertices " << p_edge.first << "," << p_edge.second
         << " both map to " << tv1;
      continue;
    }
    const auto t_edge = get_edge(tv1, tv2);
    const auto t_weight_opt =
        get_optional_value(target_edges_and_weights, t_edge);
    if (!t_weight_opt) {
      ss << "\nP-edge [" << p_edge.first << "," << p_edge.second
         << "] maps to nonexistent target edge [" << tv1 << "," << tv2 << "]";
      continue;
    }
    weight_checks.add_scalar_product(p_weight, t_weight_opt.value());
  }
  weight_checks.final_check(*this);

  if (p_vertices_used.size() != assignments.size()) {
    ss << "\nnumber of used p vertices mismatch: " << p_vertices_used.size()
       << "," << assignments.size();
  }
  return ss.str();
}

SolutionWSM::SolutionWSM()
    : complete(false),
      total_scalar_product_weight(0),
      total_p_edges_weight(0) {}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
