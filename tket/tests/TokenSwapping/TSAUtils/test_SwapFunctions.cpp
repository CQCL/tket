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

#include "TokenSwapping/SwapFunctions.hpp"

using Catch::Matchers::ContainsSubstring;

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

SCENARIO("Get swaps, with exceptions") {
  for (size_t ii = 0; ii < 5; ++ii) {
    for (size_t jj = 0; jj < 5; ++jj) {
      try {
        const auto swap = get_swap(ii, jj);
        CHECK(ii != jj);
        CHECK(swap.first == std::min(ii, jj));
        CHECK(swap.second == std::max(ii, jj));
      } catch (const std::exception& e) {
        CHECK(ii == jj);
        CHECK_THAT(std::string(e.what()), ContainsSubstring("equal vertices"));
      }
    }
  }
}

SCENARIO("Disjoint swaps") {
  std::vector<Swap> swaps;
  for (size_t ii = 0; ii < 5; ++ii) {
    for (size_t jj = ii + 1; jj < 5; ++jj) {
      swaps.push_back(get_swap(ii, jj));
    }
  }
  std::stringstream disjoint_pairs;
  std::stringstream non_disjoint_pairs;
  for (const auto& swap1 : swaps) {
    for (const auto& swap2 : swaps) {
      auto& ss = disjoint(swap1, swap2) ? disjoint_pairs : non_disjoint_pairs;
      ss << "[" << swap1.first << swap1.second << " " << swap2.first
         << swap2.second << "] ";
    }
  }
  CHECK(
      disjoint_pairs.str() ==
      "[01 23] [01 24] [01 34] [02 13] [02 14] [02 34] [03 12] [03 14] [03 24] "
      "[04 "
      "12] [04 13] [04 23] [12 03] [12 04] [12 34] [13 02] [13 04] [13 24] [14 "
      "02] "
      "[14 03] [14 23] [23 01] [23 04] [23 14] [24 01] [24 03] [24 13] [34 01] "
      "[34 "
      "02] [34 12] ");
  CHECK(
      non_disjoint_pairs.str() ==
      "[01 01] [01 02] [01 03] [01 04] [01 12] [01 13] [01 14] [02 01] [02 02] "
      "[02 "
      "03] [02 04] [02 12] [02 23] [02 24] [03 01] [03 02] [03 03] [03 04] [03 "
      "13] "
      "[03 23] [03 34] [04 01] [04 02] [04 03] [04 04] [04 14] [04 24] [04 34] "
      "[12 "
      "01] [12 02] [12 12] [12 13] [12 14] [12 23] [12 24] [13 01] [13 03] [13 "
      "12] "
      "[13 13] [13 14] [13 23] [13 34] [14 01] [14 04] [14 12] [14 13] [14 14] "
      "[14 "
      "24] [14 34] [23 02] [23 03] [23 12] [23 13] [23 23] [23 24] [23 34] [24 "
      "02] "
      "[24 04] [24 12] [24 14] [24 23] [24 24] [24 34] [34 03] [34 04] [34 13] "
      "[34 "
      "14] [34 23] [34 24] [34 34] ");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
