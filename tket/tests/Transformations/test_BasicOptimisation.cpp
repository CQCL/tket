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
#include <iostream>

#include "Circuit/Circuit.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/PQPSquash.hpp"

namespace tket {
namespace test_BasicOptimisation {

SCENARIO(
    "Transforms::remove_redundancies removes 1 and 2 qubit identities from a "
    "simple two qubit circuit") {
  Circuit original_circuit(2);
  Circuit test_circuit(original_circuit);

  GIVEN("noop") {
    test_circuit.add_op<unsigned>(OpType::noop, {0});
    WHEN("calling Transforms::remove_redundancies on rcuit") {
      Transforms::remove_redundancies().apply(test_circuit);
      THEN("added gates should be removed") {
        REQUIRE(test_circuit.circuit_equality(original_circuit));
      }
    }
  }
}

SCENARIO("Transforms::remove_redundancies removes swaps") {
  Circuit original_circuit(2);
  Circuit test_circuit(original_circuit);

  GIVEN("two consecutive identical swaps are added") {
    test_circuit.add_op<unsigned>(OpType::SWAP, {1, 0});
    test_circuit.add_op<unsigned>(OpType::SWAP, {1, 0});
    WHEN("calling Transforms::remove_redundancies on rcuit") {
      Transforms::remove_redundancies().apply(test_circuit);
      THEN("added gates should be removed") {
        REQUIRE(test_circuit.circuit_equality(original_circuit));
      }
    }
  }
  //  GIVEN("two consecutive mirrored swaps are added") {
  //    test_circuit.add_op<unsigned>(OpType::SWAP, {0, 1});
  //    test_circuit.add_op<unsigned>(OpType::SWAP, {1, 0});
  //    WHEN("calling Transforms::remove_redundancies on circuit") {
  //      Transforms::remove_redundancies().apply(test_circuit);
  //      THEN("added gates should be removed"){
  //        REQUIRE(test_circuit.circuit_equality(original_circuit));
  //      }
  //    }
  //  }
  //  GIVEN("two consecutive mirrored swaps are added") {
  //    test_circuit.add_op<unsigned>(OpType::SWAP, {0, 1});
  //    test_circuit.add_op<unsigned>(OpType::SWAP, {1, 0});
  //    WHEN("calling Transforms::remove_redundancies on circuit") {
  //      Transforms::remove_redundancies().apply(test_circuit);
  //      THEN("added gates should be removed"){
  //        REQUIRE(test_circuit.circuit_equality(original_circuit));
  //      }
  //    }
  //  }
}

}  // namespace test_BasicOptimisation
}  // namespace tket
