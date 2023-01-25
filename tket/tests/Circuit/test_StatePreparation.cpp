// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <boost/dynamic_bitset.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../testutil.hpp"
#include "Circuit/Boxes.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Gate/Rotation.hpp"
#include "Simulation/CircuitSimulator.hpp"

namespace tket {
namespace test_StatePreparation {

SCENARIO("Test bloch sphere coordinates") {
  GIVEN("Test decomposition is correct") {
    std::vector<Eigen::Vector2cd> test_states;
    test_states.push_back({1, 0});
    test_states.push_back({0, 1});
    test_states.push_back({std::sqrt(0.5), std::sqrt(0.5)});
    test_states.push_back({std::sqrt(0.5), -std::sqrt(0.5)});
    for (unsigned i = 0; i < 10; i++) {
      test_states.push_back(random_state(2, i));
    }
    for (auto state : test_states) {
      Complex a = state[0];
      Complex b = state[1];
      double theta, phi, t;
      std::tie(theta, phi, t) = get_bloch_coordinate_from_state(a, b);
      REQUIRE(
          std::abs(std::cos(theta * PI * 0.5) * std::exp(i_ * t * PI) - a) <
          EPS);
      REQUIRE(
          std::abs(
              std::sin(theta * PI * 0.5) * std::exp(i_ * (t + phi) * PI) - b) <
          EPS);
    }
  }
  GIVEN("Unnormalised vector") {
    Complex a = 0.6;
    Complex b = 2;
    REQUIRE_THROWS_MATCHES(
        get_bloch_coordinate_from_state(a, b), std::invalid_argument,
        MessageContains("unnormalised"));
  }
}

}  // namespace test_StatePreparation
}  // namespace tket
