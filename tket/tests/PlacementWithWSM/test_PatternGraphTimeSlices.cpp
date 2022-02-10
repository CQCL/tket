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

#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <sstream>

#include "PlacementWithWSM/PatternGraphTimeSlices.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

SCENARIO("Fixed examples") {
  const std::vector<std::set<VertexWSM>> gates{
      {0, 1}, {0, 2}, {1, 3}, {0, 4}, {1, 3, 4}, {2, 3}, {3, 4}};
  const PatternGraphTimeSlices slices(gates);
  std::stringstream ss;
  for (const auto& pair_list : slices.time_sliced_data) {
    ss << "[";
    for (const auto& pair : pair_list) {
      ss << pair.first << pair.second << " ";
    }
    ss << "]";
  }
  CHECK(ss.str() == "[01 ][02 13 ][04 ][13 34 ][23 ][34 ]");
  PatternGraphTimeSlices::WeightParameters parameters;
  CHECK(
      str(slices.get_weights(parameters)) ==
      "6 edges with weights: [  (0,1: 1000),  (0,2: 840),  (0,4: 680),  "
      "(1,3: 1360),  (2,3: 360),  (3,4: 720), ]\n"
      "5 vertices: {0 1 2 3 4 }\n");
}

}  // namespace tket
