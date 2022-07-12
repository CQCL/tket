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

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

#include <set>
#include <sstream>

namespace tket {
namespace WeightedSubgraphMonomorphism {

std::vector<VertexWSM> get_vertices(const GraphEdgeWeights& edges_and_weights) {
  std::set<VertexWSM> vertices;
  for (const auto& entry : edges_and_weights) {
    if (!(entry.first.first < entry.first.second)) {
      throw std::runtime_error(
          "get_vertices called on invalid raw edges_and_weights");
    }
    vertices.insert(entry.first.first);
    vertices.insert(entry.first.second);
  }
  return {vertices.cbegin(), vertices.cend()};
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
  const auto vertices = get_vertices(gdata);
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

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
