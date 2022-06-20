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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <sstream>
#include <stdexcept>

#include "Architecture/DistancesFromArchitecture.hpp"

using Catch::Matchers::ContainsSubstring;
using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

SCENARIO("Architecture with disconnected graph") {
  // Check that distance(v1, v2) does indeed give an error if v1, v2 are in
  // different connected components.
  const std::vector<std::pair<unsigned, unsigned>> edges{
      {0, 1}, {0, 2}, {1, 3}, {4, 5}};
  const size_t number_of_vertices = 6;
  const Architecture arch(edges);
  // Note: it's a "coincidence" that the vertex numbers are unchanged,
  // because 0,1,2,3,4,5 are first seen in this order.
  const ArchitectureMapping mapping(arch, edges);
  REQUIRE(mapping.number_of_vertices() == number_of_vertices);
  DistancesFromArchitecture dist_calculator(mapping);
  std::stringstream summary;
  for (size_t v1 = 0; v1 < number_of_vertices; ++v1) {
    for (size_t v2 = 0; v2 < number_of_vertices; ++v2) {
      summary << "d(" << v1 << "," << v2 << ")=";
      try {
        const auto distance = dist_calculator(v1, v2);
        summary << distance << ";";
        if (distance == 0) {
          CHECK(v1 == v2);
        } else {
          CHECK(v1 != v2);
        }
      } catch (const std::exception& e) {
        // 4 or 5 is involved, but not (4,5).
        const bool four_or_five_occurs =
            (v1 == 4 || v2 == 4 || v1 == 5 || v2 == 5);
        CHECK(four_or_five_occurs);
        // ...but not (4,5).
        CHECK(v1 + v2 != 9);
        summary << "INF;";
        const std::string message = e.what();
        CHECK_THAT(message, ContainsSubstring("are not connected"));
      }
    }
  }
  CHECK(
      summary.str() ==
      "d(0,0)=0;d(0,1)=1;d(0,2)=1;d(0,3)=2;d(0,4)=INF;d(0,5)=INF;d(1,0)=1;"
      "d(1,1)=0;"
      "d(1,2)=2;d(1,3)=1;d(1,4)=INF;d(1,5)=INF;d(2,0)=1;d(2,1)=2;d(2,2)=0;"
      "d(2,3)=3;d"
      "(2,4)=INF;d(2,5)=INF;d(3,0)=2;d(3,1)=1;d(3,2)=3;d(3,3)=0;d(3,4)=INF;"
      "d(3,5)="
      "INF;d(4,0)=INF;d(4,1)=INF;d(4,2)=INF;d(4,3)=INF;d(4,4)=0;d(4,5)=1;d("
      "5,0)=INF;"
      "d(5,1)=INF;d(5,2)=INF;d(5,3)=INF;d(5,4)=1;d(5,5)=0;");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
