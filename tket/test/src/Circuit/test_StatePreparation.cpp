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

#include <Eigen/Core>
#include <boost/dynamic_bitset.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

#include "../testutil.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Circuit/StatePreparation.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/OpType/OpType.hpp"

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

SCENARIO("Test state preparation") {
  GIVEN("n qubit states") {
    std::vector<Eigen::VectorXcd> test_states;
    test_states.push_back(Eigen::Vector2cd(0, 1));
    test_states.push_back(Eigen::Vector2cd(1, 0));
    test_states.push_back(Eigen::Vector2cd(std::sqrt(0.5), std::sqrt(0.5)));
    test_states.push_back(Eigen::Vector2cd(std::sqrt(0.5), -std::sqrt(0.5)));
    test_states.push_back(Eigen::Vector4cd(0, 1, 0, 0));
    test_states.push_back(Eigen::Vector4cd(0, 0, 0, 1));
    test_states.push_back(Eigen::Vector4cd(1, 0, 0, 0));
    test_states.push_back(Eigen::Vector4cd(
        -std::sqrt(0.25), std::sqrt(0.25), std::sqrt(0.25), -std::sqrt(0.25)));
    for (unsigned i = 0; i < 5; i++) {
      test_states.push_back(random_state(8, i));
      test_states.push_back(random_state(16, i));
      test_states.push_back(random_state(32, i));
      test_states.push_back(random_state(64, i));
    }
    for (auto psi : test_states) {
      // test state prep is correct
      StatePreparationBox prep(psi);
      std::shared_ptr<Circuit> c = prep.to_circuit();
      const Eigen::VectorXcd sv = tket_sim::get_statevector(*c);
      REQUIRE((psi - sv).cwiseAbs().sum() < ERR_EPS);
      // test the inverse of state prep is correct
      StatePreparationBox inverse_prep(psi, true);
      std::shared_ptr<Circuit> d = inverse_prep.to_circuit();
      const Eigen::MatrixXcd U = tket_sim::get_unitary(*d);
      const Eigen::VectorXcd final_state = U * psi;
      // the final state should be the computational basis 0 state
      REQUIRE(std::abs(Complex(1, 0) - final_state[0]) < ERR_EPS);
      for (unsigned i = 1; i < psi.size(); i++) {
        REQUIRE(std::abs(final_state[i]) < ERR_EPS);
      }
    }
  }
  GIVEN("unnormalised vector") {
    Eigen::Vector2cd state(1, 1);
    REQUIRE_THROWS_MATCHES(
        StatePreparationBox(state), std::invalid_argument,
        MessageContains("not normalised"));
  }
  GIVEN("vectors with wrong size") {
    Eigen::VectorXcd state(1);
    state << 1;
    REQUIRE_THROWS_MATCHES(
        StatePreparationBox(state), std::invalid_argument,
        MessageContains("not a power of 2"));
    Eigen::Vector3cd state2(1, 0, 0);
    REQUIRE_THROWS_MATCHES(
        StatePreparationBox(state2), std::invalid_argument,
        MessageContains("not a power of 2"));
  }
  GIVEN("test dagger") {
    Eigen::Vector2cd state(std::sqrt(0.5), -std::sqrt(0.5));
    StatePreparationBox prep(state);
    const StatePreparationBox dag_box =
        static_cast<const StatePreparationBox &>(*prep.dagger());
    REQUIRE((state - dag_box.get_statevector()).cwiseAbs().sum() < ERR_EPS);
    REQUIRE(dag_box.is_inverse() == true);
  }
  GIVEN("copy constructor") {
    Eigen::Vector2cd state(std::sqrt(0.5), -std::sqrt(0.5));
    StatePreparationBox prep(state, true);
    StatePreparationBox prep2(prep);
    REQUIRE((state - prep2.get_statevector()).cwiseAbs().sum() < ERR_EPS);
    REQUIRE(prep2.is_inverse() == true);
  }
  GIVEN("box with initial reset") {
    Eigen::Vector4cd state(0, 1, 0, 0);
    StatePreparationBox prep(state, false, true);
    REQUIRE(prep.with_initial_reset());
    REQUIRE_THROWS_MATCHES(
        prep.dagger(), std::logic_error, MessageContains("with initial reset"));
    Circuit c(3);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_box(prep, {0, 1});
    REQUIRE(c.count_gates(OpType::Reset) == 0);
    c.decompose_boxes_recursively();
    REQUIRE(c.count_gates(OpType::Reset) == 2);
  }
}

}  // namespace test_StatePreparation
}  // namespace tket
