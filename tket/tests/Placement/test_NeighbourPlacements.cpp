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

#include <Placement/NeighbourPlacements.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>
#include <utility>

#include "../testutil.hpp"

namespace tket {
namespace test_NeighbourPlacements {

using Connection = Architecture::Connection;

SCENARIO("class NeighbourPlacments") {
  GIVEN("a realistic-ish instance") {
    Architecture arc(
        {{Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)},
         {Node(5), Node(7)}});
    qubit_mapping_t map(
        {{Qubit(0), Node(4)},
         {Qubit(1), Node(5)},
         {Qubit(2), Node(6)},
         {Qubit(3), Node(7)}});
    NeighbourPlacements np(arc, map);

    WHEN("Getting a placement dist=0") {
      auto res = np.get(0, 1);
      THEN("There is a single result") { REQUIRE(res.size() == 1); }
      auto [new_map, swaps] = res.front();
      THEN("The resulting map is identical") {
        for (auto [k, v] : map) {
          REQUIRE(new_map.contains(k));
          REQUIRE(new_map[k] == v);
        }
      }
    }

    WHEN("Getting a placement dist=2, optimise=true") {
      auto res = np.get(2, 1);
      THEN("There is a single result") { REQUIRE(res.size() == 1); }
      auto [new_map, swaps] = res.front();
      THEN("The results are valid") {
        REQUIRE(new_map.size() == 4);
        REQUIRE(swaps.size() == 2);
        for (unsigned i = 0; i < 4; ++i) {
          REQUIRE(new_map.contains(Qubit(i)));
        }
      }
      THEN("The resulting map is correct") {
        REQUIRE(new_map[Qubit(0)] == Node(4));
        REQUIRE(new_map[Qubit(1)] == Node(7));
        REQUIRE(new_map[Qubit(2)] == Node(5));
        REQUIRE(new_map[Qubit(3)] == Node(6));
      }
      THEN("The swaps are correct") {
        REQUIRE(swaps[0] == std::pair<Node, Node>{Node(5), Node(7)});
        REQUIRE(swaps[1] == std::pair<Node, Node>{Node(5), Node(6)});
      }
    }
    WHEN("Getting 10 placement dist=3, optimise=true") {
      auto res = np.get(3, 10);
      THEN("There are 10 resulting placements") { REQUIRE(res.size() == 10); }
    }
  }
  GIVEN("the simplest possible instance") {
    Architecture arc(std::vector<std::pair<Node, Node>>{{Node(0), Node(1)}});
    qubit_mapping_t map({{Qubit(0), Node(0)}, {Qubit(1), Node(1)}});
    NeighbourPlacements np(arc, map);
    WHEN("Getting a placement dist=2, optimise=false") {
      auto res = np.get(2, 1, false);
      THEN("There is a single result") { REQUIRE(res.size() == 1); }
      auto [new_map, swaps] = res.front();
      THEN("Both swaps are identical") {
        REQUIRE(swaps.size() == 2);
        REQUIRE(swaps[0] == swaps[1]);
      }
    }
    WHEN("Getting a placement dist=2, optimise=true") {
      THEN("Con only find a solution with dist=1") {
        auto res = np.get(2, 1, true);
        REQUIRE(res.size() == 1);
        REQUIRE(res.front().swaps.size() == 1);
      }
    }
    WHEN("Getting two placements of dist=1") {
      THEN("Can only find one result") {
        REQUIRE(np.get(1, 2, false, 100).size() == 1);
      }
    }
  }
  GIVEN("an instance with unlucky seed") {
    Architecture arc({{Node(0), Node(1)}, {Node(1), Node(2)}});
    qubit_mapping_t map(
        {{Qubit(0), Node(0)}, {Qubit(1), Node(1)}, {Qubit(2), Node(2)}});
    NeighbourPlacements np(arc, map);

    // find unlucky seed
    unsigned seed;
    for (seed = 0; seed < 10; ++seed) {
      auto res = np.get(2, 1, false, seed);
      THEN("There is a single result") { REQUIRE(res.size() == 1); }
      auto [new_map, swaps] = res.front();
      REQUIRE(swaps.size() == 2);
      if (swaps[0] == swaps[1]) {
        break;
      }
    }
    THEN("There is an unlucky seed") { REQUIRE(seed < 10u); }

    WHEN("Getting a placement dist=2, optimise=false and fixed seed") {
      auto res = np.get(2, 1, false, seed);
      THEN("There is a single result") { REQUIRE(res.size() == 1); }
      auto [new_map, swaps] = res.front();
      THEN("Both swaps are identical") {
        REQUIRE(swaps.size() == 2);
        REQUIRE(swaps[0] == swaps[1]);
      }
    }
    WHEN("Getting a placement dist=2, optimise=true and fixed seed") {
      auto res = np.get(2, 1, true, seed);
      THEN("There is a single result") { REQUIRE(res.size() == 1); }
      auto [new_map, swaps] = res.front();
      THEN("Both swaps are now different") {
        REQUIRE(swaps.size() == 2);
        REQUIRE(swaps[0] != swaps[1]);
      }
    }
  }
}

}  // namespace test_NeighbourPlacements
}  // namespace tket
