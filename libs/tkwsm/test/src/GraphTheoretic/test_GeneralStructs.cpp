// Copyright Quantinuum
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
#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("get_edge, get_vertices on invalid input") {
  for (VertexWSM ii = 0; ii <= 1; ++ii) {
    for (VertexWSM jj = 0; jj <= 1; ++jj) {
      if (ii == jj) {
        REQUIRE_THROWS_AS(get_edge(ii, jj), std::runtime_error);
        continue;
      }
      const auto edge = get_edge(ii, jj);
      REQUIRE(edge.first == std::min(ii, jj));
      REQUIRE(edge.second == std::max(ii, jj));
    }
  }
  GraphEdgeWeights gdata;
  gdata[std::make_pair<VertexWSM, VertexWSM>(1, 0)] = 0;
  REQUIRE_THROWS_AS(get_vertices(gdata), std::runtime_error);
}

SCENARIO("general structs string functions") {
  GraphEdgeWeights gdata;
  gdata[std::make_pair<VertexWSM, VertexWSM>(0, 1)] = 2;
  CHECK(
      str(gdata) ==
      "1 edges with weights: [  (0,1: 2), ]\n2 vertices: {0 1 }\n");

  std::vector<EdgeWSM> assignments;
  assignments.emplace_back(11, 22);
  assignments.emplace_back(33, 44);
  CHECK(str(assignments) == "[ 11:22  33:44 ]");

  Assignments assignments_map;
  assignments_map[7] = 8;
  CHECK(str(assignments_map) == "[ 7:8 ]");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
