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

#include <Eigen/Core>
#include <boost/dynamic_bitset.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>

#include "../testutil.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/ConjugationBox.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Gate/Rotation.hpp"

namespace tket {
namespace test_ConjugationBox {

SCENARIO("Test ConjugationBox") {
  GIVEN("Constructor with default uncompute") {
    Circuit compute(2);
    compute.add_op<unsigned>(OpType::CRx, 0.5, {1, 0});
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    action.add_op<unsigned>(OpType::H, {0});
    Op_ptr action_op = std::make_shared<CircBox>(action);
    ConjugationBox box(compute_op, action_op);
    std::shared_ptr<Circuit> c = box.to_circuit();
    Circuit d(2);
    d.add_op<unsigned>(compute_op, {0, 1});
    d.add_op<unsigned>(action_op, {0, 1});
    d.add_op<unsigned>(compute_op->dagger(), {0, 1});
    REQUIRE(*c == d);
  }
  GIVEN("Constructor with explicit uncompute op") {
    Circuit compute(2);
    compute.add_op<unsigned>(OpType::CX, {0, 1});
    compute.add_op<unsigned>(OpType::CX, {1, 0});
    compute.add_op<unsigned>(OpType::CX, {0, 1});
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    action.add_op<unsigned>(OpType::H, {0});
    Op_ptr action_op = std::make_shared<CircBox>(action);
    Circuit uncompute(2);
    uncompute.add_op<unsigned>(OpType::CX, {1, 0});
    uncompute.add_op<unsigned>(OpType::CX, {0, 1});
    uncompute.add_op<unsigned>(OpType::CX, {1, 0});
    Op_ptr uncompute_op = std::make_shared<CircBox>(uncompute);
    ConjugationBox box(compute_op, action_op, uncompute_op);
    std::shared_ptr<Circuit> c = box.to_circuit();
    Circuit d(2);
    d.add_op<unsigned>(compute_op, {0, 1});
    d.add_op<unsigned>(action_op, {0, 1});
    d.add_op<unsigned>(uncompute_op, {0, 1});
    REQUIRE(*c == d);
  }
  GIVEN("Test dagger") {
    Circuit compute(2);
    compute.add_op<unsigned>(OpType::CRx, 0.5, {1, 0});
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    action.add_op<unsigned>(OpType::Rz, 0.5, {0});
    Op_ptr action_op = std::make_shared<CircBox>(action);
    WHEN("with default uncompute") {
      ConjugationBox box(compute_op, action_op);
      ConjugationBox correct_box_dagger(compute_op, action_op->dagger());
      const ConjugationBox box_dagger =
          static_cast<const ConjugationBox &>(*box.dagger());
      REQUIRE(box_dagger == correct_box_dagger);
    }
    WHEN("with explicit uncompute") {
      ConjugationBox box(compute_op, action_op, compute_op->dagger());
      ConjugationBox correct_box_dagger(
          compute_op, action_op->dagger(), compute_op->dagger());
      const ConjugationBox box_dagger =
          static_cast<const ConjugationBox &>(*box.dagger());
      REQUIRE(box_dagger == correct_box_dagger);
    }
  }
  GIVEN("Test transpose") {
    Circuit compute(1);
    compute.add_op<unsigned>(OpType::TK1, {0.1, 0.2, 0.3}, {0});
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(1);
    action.add_op<unsigned>(OpType::TK1, {1.1, 1.2, 1.3}, {0});
    Op_ptr action_op = std::make_shared<CircBox>(action);
    WHEN("with default uncompute") {
      ConjugationBox box(compute_op, action_op);
      ConjugationBox correct_box_transpose(
          compute_op->dagger()->transpose(), action_op->transpose(),
          compute_op->transpose());
      const ConjugationBox box_transpose =
          static_cast<const ConjugationBox &>(*box.transpose());
      REQUIRE(box_transpose == correct_box_transpose);
    }
    WHEN("with explicit uncompute") {
      ConjugationBox box(compute_op, action_op, compute_op->dagger());
      ConjugationBox correct_box_transpose(
          compute_op->dagger()->transpose(), action_op->transpose(),
          compute_op->transpose());
      const ConjugationBox box_transpose =
          static_cast<const ConjugationBox &>(*box.transpose());
      REQUIRE(box_transpose == correct_box_transpose);
    }
  }
}
SCENARIO("Test ConjugationBox Exceptions") {
  GIVEN("Ops with classical wires") {
    Circuit compute(2, 1);
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    Op_ptr action_op = std::make_shared<CircBox>(action);
    REQUIRE_THROWS_MATCHES(
        ConjugationBox(compute_op, action_op), std::invalid_argument,
        MessageContains("only supports quantum operations"));
  }
  GIVEN("Uncompute with classical wires") {
    Circuit compute(2);
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    Op_ptr action_op = std::make_shared<CircBox>(action);
    Circuit uncompute(2, 1);
    Op_ptr uncompute_op = std::make_shared<CircBox>(uncompute);
    REQUIRE_THROWS_MATCHES(
        ConjugationBox(compute_op, action_op, uncompute_op),
        std::invalid_argument,
        MessageContains("only supports quantum operations"));
  }
  GIVEN("Unmatched size") {
    Circuit compute(3);
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    Op_ptr action_op = std::make_shared<CircBox>(action);
    REQUIRE_THROWS_MATCHES(
        ConjugationBox(compute_op, action_op), std::invalid_argument,
        MessageContains("have the same number of qubits"));
  }
  GIVEN("Unmatched size caused by uncompute") {
    Circuit compute(2);
    Op_ptr compute_op = std::make_shared<CircBox>(compute);
    Circuit action(2);
    Op_ptr action_op = std::make_shared<CircBox>(action);
    Circuit uncompute(3);
    Op_ptr uncompute_op = std::make_shared<CircBox>(uncompute);
    REQUIRE_THROWS_MATCHES(
        ConjugationBox(compute_op, action_op, uncompute_op),
        std::invalid_argument,
        MessageContains("have the same number of qubits"));
  }
}
}  // namespace test_ConjugationBox
}  // namespace tket
