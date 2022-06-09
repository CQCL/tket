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

#include "Utils/HelperFunctions.hpp"

namespace tket {
namespace test_Utils {

SCENARIO("Test GrayCode generation") {
  WHEN("Try 2 control graycode") {
    GrayCode gc = gen_graycode(2);
    GrayCode correct_gc = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    REQUIRE(gc == correct_gc);
  }
  WHEN("Try 3 control graycode") {
    GrayCode gc = gen_graycode(3);
    GrayCode correct_gc = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                           {0, 1, 1}, {1, 1, 1}, {1, 0, 1}, {0, 0, 1}};
    REQUIRE(gc == correct_gc);
  }
  WHEN("Try 4 control graycode") {
    GrayCode gc = gen_graycode(4);
    GrayCode correct_gc = {
        {0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, 0},
        {0, 1, 1, 0}, {1, 1, 1, 0}, {1, 0, 1, 0}, {0, 0, 1, 0},
        {0, 0, 1, 1}, {1, 0, 1, 1}, {1, 1, 1, 1}, {0, 1, 1, 1},
        {0, 1, 0, 1}, {1, 1, 0, 1}, {1, 0, 0, 1}, {0, 0, 0, 1}};
    REQUIRE(gc == correct_gc);
  }
}

}  // namespace test_Utils
}  // namespace tket
