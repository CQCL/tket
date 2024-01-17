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
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>

#include "../testutil.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Circuit/ToffoliBox.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/Utils/HelperFunctions.hpp"

namespace tket {
namespace test_ToffoliBox {

state_perm_t random_permutation(unsigned n_qubits, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::vector<unsigned> ints(1 << n_qubits);
  std::iota(std::begin(ints), std::end(ints), 0);
  std::shuffle(ints.begin(), ints.end(), rng);
  state_perm_t perm;
  for (unsigned i = 0; i < ints.size(); i++) {
    perm.insert({dec_to_bin(i, n_qubits), dec_to_bin(ints[i], n_qubits)});
  }
  return perm;
}

Eigen::MatrixXcd permutation_matrix(const state_perm_t &perm) {
  // use xcd for comparison to circuit unitary
  unsigned n_qubits = perm.begin()->first.size();
  Eigen::MatrixXcd u = Eigen::MatrixXcd::Zero(1 << n_qubits, 1 << n_qubits);
  for (unsigned i = 0; i < (unsigned)(1 << n_qubits); i++) {
    auto it = perm.find(dec_to_bin(i, n_qubits));
    if (it == perm.end()) {
      u(i, i) = 1;
    } else {
      u(bin_to_dec(it->second), i) = 1;
    }
  }
  return u;
}

SCENARIO("Test ToffoliBox") {
  state_perm_t perm;
  OpType axis = OpType::Rx;
  ToffoliBoxSynthStrat strat;
  GIVEN("1-q permutation") {
    perm[{0}] = {1};
    perm[{1}] = {0};
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("2-q permutation") {
    perm[{0, 1}] = {1, 1};
    perm[{1, 0}] = {0, 1};
    perm[{1, 1}] = {1, 0};
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("2-q permutation (2)") {
    perm[{0, 0}] = {1, 1};
    perm[{1, 1}] = {0, 0};
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("2-q permutation (3)") {
    perm[{0, 0}] = {1, 1};
    perm[{1, 1}] = {0, 0};
    perm[{0, 1}] = {1, 0};
    perm[{1, 0}] = {0, 1};
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("2-q permutation (4)") {
    perm[{0, 0}] = {1, 1};
    perm[{1, 1}] = {0, 1};
    perm[{0, 1}] = {1, 0};
    perm[{1, 0}] = {0, 0};
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("3-q permutation") {
    perm[{0, 0, 0}] = {1, 0, 0};
    perm[{0, 1, 0}] = {1, 0, 1};
    perm[{0, 1, 1}] = {0, 1, 0};
    perm[{1, 0, 0}] = {0, 0, 0};
    perm[{1, 0, 1}] = {0, 1, 1};
    perm[{1, 1, 0}] = {1, 1, 1};
    perm[{1, 1, 1}] = {1, 1, 0};
    WHEN("Use Matching strat") {
      strat = ToffoliBoxSynthStrat::Matching;
      axis = OpType::Ry;
    }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("3-q permutation (2)") {
    perm[{0, 0, 1}] = {1, 1, 0};
    perm[{1, 1, 0}] = {0, 1, 0};
    perm[{0, 1, 0}] = {1, 0, 1};
    perm[{1, 0, 1}] = {0, 0, 1};
    WHEN("Use Matching strat") {
      strat = ToffoliBoxSynthStrat::Matching;
      axis = OpType::Ry;
    }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("4-q permutation") {
    perm[{0, 0, 0, 0}] = {1, 1, 0, 0};
    perm[{1, 1, 0, 0}] = {1, 1, 0, 1};
    perm[{1, 1, 0, 1}] = {0, 0, 0, 1};
    perm[{0, 0, 0, 1}] = {1, 1, 1, 0};
    perm[{1, 1, 1, 0}] = {0, 0, 1, 1};
    perm[{0, 0, 1, 1}] = {1, 0, 0, 1};
    perm[{1, 0, 0, 1}] = {1, 0, 1, 0};
    perm[{1, 0, 1, 0}] = {0, 0, 0, 0};
    WHEN("Use Matching strat") {
      strat = ToffoliBoxSynthStrat::Matching;
      axis = OpType::Ry;
    }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("Random 4-q permutation") {
    perm = random_permutation(4, 1);
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("Random 5-q permutation") {
    perm = random_permutation(5, 1);
    WHEN("Use Matching strat") { strat = ToffoliBoxSynthStrat::Matching; }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  GIVEN("Random 6-q permutation") {
    perm = random_permutation(6, 1);
    WHEN("Use Matching strat") {
      strat = ToffoliBoxSynthStrat::Matching;
      axis = OpType::Ry;
    }
    WHEN("Use Cycle strat") { strat = ToffoliBoxSynthStrat::Cycle; }
  }
  ToffoliBox box(perm, strat, axis);
  Circuit circ = *box.to_circuit();
  const auto matrix = tket_sim::get_unitary(circ);
  const auto perm_matrix = permutation_matrix(perm);
  REQUIRE((matrix - perm_matrix).cwiseAbs().sum() < ERR_EPS);
}

SCENARIO("Test ToffoliBox Exceptions") {
  GIVEN("Invalid permutation") {
    state_perm_t perm;
    perm[{0, 1}] = {1, 0};
    REQUIRE_THROWS_MATCHES(
        ToffoliBox(perm), std::invalid_argument,
        MessageContains("is not complete"));
  }
  GIVEN("Empty permutation") {
    state_perm_t perm;
    REQUIRE_THROWS_MATCHES(
        ToffoliBox(perm), std::invalid_argument, MessageContains("empty"));
  }
  GIVEN("Wrong axis") {
    state_perm_t perm;
    perm[{0}] = {1};
    perm[{1}] = {0};
    REQUIRE_THROWS_MATCHES(
        ToffoliBox(perm, ToffoliBoxSynthStrat::Matching, OpType::Rz),
        std::invalid_argument, MessageContains("must be Rx or Ry"));
  }
  GIVEN("Invalid entries") {
    state_perm_t perm;
    perm[{0}] = {1, 0};
    REQUIRE_THROWS_MATCHES(
        ToffoliBox(perm), std::invalid_argument,
        MessageContains("with different sizes"));
  }
  GIVEN("Too long") {
    state_perm_t perm;
    std::vector<bool> b(33, 0);
    perm[b] = b;
    REQUIRE_THROWS_MATCHES(
        ToffoliBox(perm), std::invalid_argument,
        MessageContains("up to 32 bits"));
  }
}

SCENARIO("Test constructors & transformations") {
  state_perm_t perm;
  perm[{0, 1}] = {1, 1};
  perm[{1, 0}] = {0, 1};
  perm[{1, 1}] = {1, 0};
  GIVEN("copy constructor") {
    ToffoliBox box(perm, ToffoliBoxSynthStrat::Cycle, OpType::Rx);
    REQUIRE(box.get_rotation_axis() == OpType::Rx);
    REQUIRE(box.get_strat() == ToffoliBoxSynthStrat::Cycle);
    ToffoliBox box_copy(box);
    REQUIRE(box_copy.get_rotation_axis() == OpType::Rx);
    REQUIRE(box_copy.get_strat() == ToffoliBoxSynthStrat::Cycle);
    REQUIRE(box_copy.get_permutation() == perm);
  }
  GIVEN("Dagger") {
    ToffoliBox box(perm);
    Circuit circ1 = *box.to_circuit();
    Circuit circ1_dag = circ1.dagger();
    const ToffoliBox box_dag = static_cast<const ToffoliBox &>(*box.dagger());
    Circuit circ2 = *box_dag.to_circuit();
    REQUIRE(
        (tket_sim::get_unitary(circ1_dag) - tket_sim::get_unitary(circ2))
            .cwiseAbs()
            .sum() < ERR_EPS);
  }
  GIVEN("Transpose") {
    ToffoliBox box(perm);
    Circuit circ1 = *box.to_circuit();
    const ToffoliBox box_transpose =
        static_cast<const ToffoliBox &>(*box.transpose());
    Circuit circ2 = *box_transpose.to_circuit();
    auto matrix1 = tket_sim::get_unitary(circ1);
    auto matrix2 = tket_sim::get_unitary(circ2);
    REQUIRE((matrix1.transpose() - matrix2).cwiseAbs().sum() < ERR_EPS);
  }
}
}  // namespace test_ToffoliBox
}  // namespace tket
