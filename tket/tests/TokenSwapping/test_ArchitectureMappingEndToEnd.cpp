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

#include <catch2/catch_test_macros.hpp>
#include <sstream>

#include "Architecture/ArchitectureMapping.hpp"
#include "Architecture/DistancesFromArchitecture.hpp"
#include "Architecture/NeighboursFromArchitecture.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

SCENARIO("Simple path") {
  const vector<std::pair<unsigned, unsigned>> edges{
      {111, 222}, {555, 444}, {333, 222}, {777, 666}, {333, 444}, {666, 555}};
  const unsigned n_verts = edges.size() + 1;
  std::stringstream ss;
  ss << "Original input edges:\n";
  for (auto edge : edges) {
    ss << "(" << edge.first << "," << edge.second << ") ";
  }
  const Architecture arch(edges);
  const ArchitectureMapping arch_mapping(arch, edges);

  ss << "...\nEdges from arch.mapping:\n";
  for (auto edge : arch_mapping.get_edges()) {
    ss << "(" << edge.first << "," << edge.second << ") ";
  }
  ss << "...\nVertex-to-node:";

  for (unsigned vv = 0; vv < n_verts; ++vv) {
    const auto node = arch_mapping.get_node(vv);
    REQUIRE(vv == arch_mapping.get_vertex(node));
    ss << "\n" << vv << " == " << node.repr();
  }
  ss << "...\nDistances:";

  DistancesFromArchitecture distances(arch_mapping);
  NeighboursFromArchitecture neighbours(arch_mapping);

  for (unsigned ii = 0; ii < n_verts; ++ii) {
    ss << "\n" << ii << ": [";
    for (unsigned jj = ii + 1; jj < n_verts; ++jj) {
      REQUIRE(0 == distances(ii, ii));
      const auto dist = distances(ii, jj);
      ss << " " << dist;
      REQUIRE(dist == distances(jj, ii));
    }
    ss << "]";
  }
  ss << "\nNeighbours:";
  for (unsigned ii = 0; ii < n_verts; ++ii) {
    ss << "\n" << ii << ": [";
    const auto& neighb = neighbours(ii);
    for (auto nn : neighb) {
      ss << " " << nn;
    }
    ss << " ]";
  }
  CHECK(
      ss.str() ==
      "Original input edges:\n"
      "(111,222) (555,444) (333,222) (777,666) (333,444) (666,555) ...\n"
      "Edges from arch.mapping:\n"
      "(0,1) (2,3) (1,4) (5,6) (3,4) (2,6) ...\n"
      "Vertex-to-node:\n"
      "0 == node[111]\n"
      "1 == node[222]\n"
      "2 == node[555]\n"
      "3 == node[444]\n"
      "4 == node[333]\n"
      "5 == node[777]\n"
      "6 == node[666]...\n"
      "Distances:\n"
      "0: [ 1 4 3 2 6 5]\n"
      "1: [ 3 2 1 5 4]\n"
      "2: [ 1 2 2 1]\n"
      "3: [ 1 3 2]\n"
      "4: [ 4 3]\n"
      "5: [ 1]\n"
      "6: []\n"
      "Neighbours:\n"
      "0: [ 1 ]\n"
      "1: [ 0 4 ]\n"
      "2: [ 3 6 ]\n"
      "3: [ 2 4 ]\n"
      "4: [ 1 3 ]\n"
      "5: [ 6 ]\n"
      "6: [ 2 5 ]");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
