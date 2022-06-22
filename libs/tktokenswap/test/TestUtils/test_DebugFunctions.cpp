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

#include "DebugFunctions.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

SCENARIO("debug functions - string functions") {
  const VertexMapping vm{{0, 1}, {1, 2}, {3, 5}};
  CHECK(str(vm) == "VM: 0->1  1->2  3->5 ");

  vector<Swap> swaps_vect;
  swaps_vect.push_back(get_swap(111, 222));
  swaps_vect.push_back(get_swap(5555, 4444));
  const auto swaps_vect_str = str(swaps_vect);
  CHECK(swaps_vect_str == " (111,222)  (4444,5555) ");

  SwapList swaps;
  for (const auto& swap : swaps_vect) {
    swaps.push_back(swap);
  }
  CHECK(swaps_vect_str == str(swaps));
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
