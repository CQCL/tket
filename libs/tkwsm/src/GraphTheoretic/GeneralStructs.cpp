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

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

#include <set>
#include <sstream>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

unsigned get_number_of_vertices(const GraphEdgeWeights& edges_and_weights) {
  std::set<VertexWSM> vertices;
  for (const auto& entry : edges_and_weights) {
    const EdgeWSM& edge = entry.first;
    const VertexWSM& v1 = edge.first;
    const VertexWSM& v2 = edge.second;
    vertices.insert(v1);
    vertices.insert(v2);
  }
  return vertices.size();
}

std::vector<VertexWSM> get_vertices(
    const GraphEdgeWeights& edges_and_weights,
    const GetVerticesOptions& options) {
  std::set<VertexWSM> vertices;
  for (const auto& entry : edges_and_weights) {
    const EdgeWSM& edge = entry.first;
    const WeightWSM& weight = entry.second;
    const VertexWSM& v1 = edge.first;
    const VertexWSM& v2 = edge.second;
    vertices.insert(v1);
    vertices.insert(v2);
    try {
      if (weight == 0 && !options.allow_zero_weights) {
        throw std::runtime_error("Zero weight not allowed");
      }
      if (v1 == v2 && !options.allow_loops) {
        throw std::runtime_error("Loop not allowed");
      }
      const auto reversed_edge = std::make_pair(v2, v1);
      const auto reversed_weight_opt =
          get_optional_value(edges_and_weights, reversed_edge);
      if (reversed_weight_opt) {
        if (reversed_weight_opt.value() != weight) {
          throw std::runtime_error("reversed edge has different weight");
        }
        if (!options.allow_duplicate_edges) {
          throw std::runtime_error("duplicate edges not allowed");
        }
      }
      if (!options.allow_edge_vertices_not_in_order && v2 < v1) {
        throw std::runtime_error("we do not allow v2<v1 in edge (v1,v2)");
      }
    } catch (const std::exception& e) {
      std::stringstream ss;
      ss << "get_vertices called for edge->weight map, size "
         << edges_and_weights.size() << "; for edge (" << v1 << "," << v2
         << "), weight " << weight << ": " << e.what();
      throw std::runtime_error(ss.str());
    }
  }
  return {vertices.cbegin(), vertices.cend()};
}

WeightWSM get_max_weight(const GraphEdgeWeights& graph_data) {
  WeightWSM weight = 0;
  for (const auto& entry : graph_data) {
    weight = std::max(weight, entry.second);
  }
  return weight;
}

WeightWSM get_checked_scalar_product(
    const GraphEdgeWeights& pdata, const GraphEdgeWeights& tdata,
    const std::vector<std::pair<VertexWSM, VertexWSM>>& solution) {
  // Is target data valid?
  for (const auto& t_edge_and_weight : tdata) {
    const EdgeWSM& t_edge = t_edge_and_weight.first;
    TKET_ASSERT(t_edge.first != t_edge.second);
    const EdgeWSM reversed_edge{t_edge.second, t_edge.first};
    const auto weight_opt = get_optional_value(tdata, reversed_edge);
    if (weight_opt) {
      TKET_ASSERT(weight_opt.value() == t_edge_and_weight.second);
    }
  }

  const auto sorted_pv = get_vertices(pdata);
  TKET_ASSERT(solution.size() == sorted_pv.size());

  const auto sorted_tv = get_vertices(tdata);

  // Build up an assignment map.
  std::map<VertexWSM, VertexWSM> assignments_map;
  std::set<VertexWSM> used_tv;

  for (const auto& assignment : solution) {
    TKET_ASSERT(std::binary_search(
        sorted_pv.cbegin(), sorted_pv.cend(), assignment.first));
    TKET_ASSERT(std::binary_search(
        sorted_tv.cbegin(), sorted_tv.cend(), assignment.second));
    TKET_ASSERT(used_tv.count(assignment.second) == 0);
    used_tv.insert(assignment.second);
    assignments_map.emplace(assignment);
  }
  TKET_ASSERT(assignments_map.size() == solution.size());
  TKET_ASSERT(assignments_map.size() == used_tv.size());

  // Now build up the scalar product.
  WeightWSM scalar_product = 0;
  for (const auto& entry : pdata) {
    const EdgeWSM& p_edge = entry.first;
    TKET_ASSERT(p_edge.first != p_edge.second);

    const WeightWSM& p_weight = entry.second;
    {
      const auto reversed_p_edge = std::make_pair(p_edge.second, p_edge.first);
      const auto reversed_p_weight_opt =
          get_optional_value(pdata, reversed_p_edge);
      if (reversed_p_weight_opt) {
        TKET_ASSERT(reversed_p_weight_opt.value() == p_weight);
        if (p_edge.first > p_edge.second) {
          // If (pv1, pv2) and (pv2, pv1) BOTH exist, with pv1<pv2,
          // add only the contribution for (pv1, pv2).
          // This avoids counting twice.
          continue;
        }
      }
    }

    const auto t_edge = std::make_pair(
        assignments_map.at(p_edge.first), assignments_map.at(p_edge.second));

    // We already checked that the target data is at least valid.
    auto t_weight_opt = get_optional_value(tdata, t_edge);
    if (!t_weight_opt) {
      t_weight_opt = get_optional_value(
          tdata, std::make_pair(t_edge.second, t_edge.first));
    }
    TKET_ASSERT(t_weight_opt);
    const WeightWSM sc_prod_contrib =
        get_product_or_throw(p_weight, t_weight_opt.value());
    scalar_product = get_sum_or_throw(scalar_product, sc_prod_contrib);
  }
  return scalar_product;
}

EdgeWSM get_edge(VertexWSM v1, VertexWSM v2) {
  if (v1 > v2) {
    std::swap(v1, v2);
  }
  if (!(v1 < v2)) {
    throw std::runtime_error(
        std::string("get_edge on equal vertex v1=") + std::to_string(v1));
  }
  return std::make_pair(v1, v2);
}

std::string str(const GraphEdgeWeights& gdata) {
  std::stringstream ss;
  ss << gdata.size() << " edges with weights: [ ";
  for (const auto& entry : gdata) {
    ss << " (" << entry.first.first << "," << entry.first.second << ": "
       << entry.second << "), ";
  }
  ss << "]\n";
  GetVerticesOptions options;
  // Allow everything!
  options.allow_duplicate_edges = true;
  options.allow_edge_vertices_not_in_order = true;
  options.allow_loops = true;
  options.allow_zero_weights = true;
  const auto vertices = get_vertices(gdata, options);

  ss << vertices.size() << " vertices: {";
  for (auto vv : vertices) {
    ss << vv << " ";
  }
  ss << "}\n";
  return ss.str();
}

std::string str(const std::vector<EdgeWSM>& assignments) {
  std::stringstream ss;
  ss << "[";
  for (const EdgeWSM& entry : assignments) {
    ss << " " << entry.first << ":" << entry.second << " ";
  }
  ss << "]";
  return ss.str();
}

std::string str(const Assignments& assignments) {
  std::stringstream ss;
  ss << "[";
  for (const std::pair<const VertexWSM, VertexWSM>& entry : assignments) {
    ss << " " << entry.first << ":" << entry.second << " ";
  }
  ss << "]";
  return ss.str();
}

BitsetInformation::BitsetInformation(const boost::dynamic_bitset<>& domain) {
  // NOTE: we really do want auto here.
  // If no bit is set, this value will be the max value
  // of whatever type find_first() returns...
  // which is very incovenient to find without auto...
  auto first_val = domain.find_first();
  if (first_val >= domain.size()) {
    empty = true;
    return;
  }
  empty = false;
  if (domain.find_next(first_val) >= domain.size()) {
    // It has exactly one bit.
    single_element = static_cast<VertexWSM>(first_val);
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
