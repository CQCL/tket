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
#include <tkassert/Assert.hpp>
#include <tkassert/AssertMessage.hpp>

namespace tket {

SCENARIO("Assertion") {
  for (double x = -10.1; x < 10.5; x += 1.) {
    for (double y = -10.1; y < 10.5; y += 1.) {
      double z = (x + y) * (x - y) - x * x + y * y;
      double w = z * z;
      TKET_ASSERT(w >= 0);
      TKET_ASSERT(
          w < 1e-12 || AssertMessage()
                           << "Maths failure for x = " << x << ", y = " << y);
    }
  }
  REQUIRE(true);
}

}  // namespace tket
