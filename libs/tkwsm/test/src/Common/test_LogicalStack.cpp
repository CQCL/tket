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

#include <catch2/catch_test_macros.hpp>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/LogicalStack.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("Test random logical stack ops") {
  RNG rng;
  LogicalStack<int> stack;
  std::vector<int> shadowing_stack;

  for (int value = 0; value < 1000; ++value) {
    if (!shadowing_stack.empty() && rng.check_percentage(50)) {
      stack.pop();
      shadowing_stack.pop_back();
    } else {
      stack.push();
      stack.top() = value;
      shadowing_stack.push_back(value);
    }
    REQUIRE(stack.empty() == shadowing_stack.empty());
    REQUIRE(stack.size() == shadowing_stack.size());
    for (unsigned ii = 0; ii < shadowing_stack.size(); ++ii) {
      REQUIRE(stack[ii] == shadowing_stack[ii]);
    }
    if (shadowing_stack.size() >= 1) {
      REQUIRE(stack.top() == shadowing_stack.back());
      if (shadowing_stack.size() >= 2) {
        REQUIRE(
            stack.one_below_top() ==
            shadowing_stack[shadowing_stack.size() - 2]);
      }
    }
    if ((value % 20) == 0) {
      stack.clear();
      shadowing_stack.clear();
      REQUIRE(stack.size() == 0);
    }
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
