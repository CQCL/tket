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
#include <string>

#include "TokenSwapping/SwapFunctions.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

std::string get_swaps_str(const SwapList& swap_list) {
  std::stringstream ss;
  const auto svect = swap_list.to_vector();
  ss << "[" << svect.size() << " swaps:";
  for (auto swap : svect) {
    ss << " (" << swap.first << " " << swap.second << ") ";
  }
  ss << "]";
  return ss.str();
}

SCENARIO("simple swap list") {
  SwapList swap_list;
  CHECK(get_swaps_str(swap_list) == "[0 swaps:]");
  swap_list.clear();
  CHECK(get_swaps_str(swap_list) == "[0 swaps:]");
  swap_list.push_front(get_swap(0, 1));
  CHECK(get_swaps_str(swap_list) == "[1 swaps: (0 1) ]");
  const auto current_front = swap_list.front_id().value();
  const auto new_front = swap_list.emplace_front();
  CHECK(current_front != new_front);
  CHECK(new_front == swap_list.front_id().value());
  swap_list.front() = get_swap(998, 999);
  CHECK(get_swaps_str(swap_list) == "[2 swaps: (998 999)  (0 1) ]");
  swap_list.pop_front();
  CHECK(get_swaps_str(swap_list) == "[1 swaps: (0 1) ]");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
