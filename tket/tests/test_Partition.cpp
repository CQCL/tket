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

#include "Diagonalisation/PauliPartition.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_Partition {

SCENARIO("Small sets of Gadgets are partitioned correctly") {
  const std::list<PauliPartitionStrat> strats{
      PauliPartitionStrat::NonConflictingSets,
      PauliPartitionStrat::CommutingSets};

  // NOTE: all methods seem to give the same results,
  // which is not surprising for small sets.
  // It would be good for a subject expert to add
  // more extensive tests with larger sets.
  const std::vector<GraphColourMethod> colouring_methods{
      GraphColourMethod::LargestFirst, GraphColourMethod::Exhaustive,
      GraphColourMethod::Lazy};

  GIVEN("No gadgets") {
    for (auto colouring_method : colouring_methods) {
      for (PauliPartitionStrat strat : strats) {
        std::list<QubitPauliString> tensors;
        std::list<std::list<QubitPauliString>> void_terms =
            term_sequence(tensors, strat, colouring_method);
        REQUIRE(void_terms.empty());
      }
    }
  }
  GIVEN("Two anti-commuting gadgets") {
    /* We know the correct order, as QubitOperator orders
    lexicographically */
    double a = 0.3;
    double b = 53.1;
    Qubit q0(0);
    Qubit q1(1);
    Qubit q2(2);
    QubitPauliString qp_map0({{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}});
    QubitPauliString qp_map1({{q0, Pauli::Z}, {q1, Pauli::Z}, {q2, Pauli::Y}});
    std::list<QubitPauliString> tensors{qp_map0, qp_map1};

    for (auto colouring_method : colouring_methods) {
      for (PauliPartitionStrat strat : strats) {
        std::list<std::list<QubitPauliString>> terms =
            term_sequence(tensors, strat, colouring_method);

        REQUIRE(terms.size() == 2);
        REQUIRE(terms.begin()->size() == 1);
        REQUIRE(*(terms.begin()->begin()) == qp_map0);
        REQUIRE((++terms.begin())->size() == 1);
        REQUIRE(*((++terms.begin())->begin()) == qp_map1);
      }
    }
  }
  GIVEN("Three partitions of four gadgets") {
    double a = 0.3;
    double b = 53.1;
    double c = 1145;
    double d = 0.11111117;
    QubitPauliString qp_map0(Qubit(0), Pauli::I);
    QubitPauliString qp_map1(Qubit(0), Pauli::X);
    QubitPauliString qp_map2(Qubit(0), Pauli::Y);
    QubitPauliString qp_map3(Qubit(0), Pauli::Z);
    std::list<QubitPauliString> tensors{qp_map0, qp_map1, qp_map2, qp_map3};

    for (auto colouring_method : colouring_methods) {
      for (PauliPartitionStrat strat : strats) {
        std::list<std::list<QubitPauliString>> terms =
            term_sequence(tensors, strat, colouring_method);

        REQUIRE(terms.size() == 3);
        unsigned total_terms = 0;
        for (const std::list<QubitPauliString>& g_map : terms) {
          REQUIRE((g_map.size() == 1 || g_map.size() == 2));
          total_terms += g_map.size();
        }
        REQUIRE(total_terms == 4);
      }
    }
  }
}

}  // namespace test_Partition
}  // namespace tket
