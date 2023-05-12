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
#include <random>

#include "../testutil.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/DiagonalBox.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/Simulation/CircuitSimulator.hpp"

namespace tket {
namespace test_DiagonalBox {

Eigen::VectorXcd random_diagonal(unsigned n, int seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(-10, 10);
  Eigen::VectorXcd diagonal = Eigen::VectorXcd::Zero(n);
  for (unsigned i = 0; i < n; i++) {
    diagonal[i] = std::exp(dist(rng) * i_);
    REQUIRE(std::abs(1 - std::abs(diagonal[i])) < ERR_EPS);
  }
  return diagonal;
}

SCENARIO("Test DiagonalBox") {
  GIVEN("n qubit diagonals") {
    std::vector<Eigen::VectorXcd> test_diagonals;
    test_diagonals.push_back(Eigen::Vector2cd(-1, 1));
    test_diagonals.push_back(Eigen::Vector2cd(1, 1));
    test_diagonals.push_back(Eigen::Vector2cd(std::exp(0.7 * i_), 1));
    test_diagonals.push_back(
        Eigen::Vector2cd(std::exp(3.7 * i_), std::exp(-2. * i_)));
    test_diagonals.push_back(Eigen::Vector4cd(1, 1, i_, i_));
    test_diagonals.push_back(Eigen::Vector4cd(i_, -1, 1, 1));
    Eigen::VectorXcd temp_diag(8);
    temp_diag << i_, i_, i_, i_, i_, i_, i_, i_;
    test_diagonals.push_back(temp_diag);
    for (unsigned i = 0; i < 5; i++) {
      test_diagonals.push_back(random_diagonal(8, i));
      test_diagonals.push_back(random_diagonal(16, i));
      test_diagonals.push_back(random_diagonal(32, i));
      test_diagonals.push_back(random_diagonal(64, i));
    }
    for (auto d : test_diagonals) {
      // test diagonal is correct
      DiagonalBox diag(d);
      std::shared_ptr<Circuit> c = diag.to_circuit();
      const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
      Eigen::DiagonalMatrix<Complex, Eigen::Dynamic> D(d);
      Eigen::MatrixXcd M = D.derived();
      REQUIRE((U - M).cwiseAbs().sum() < ERR_EPS);
    }
  }
  GIVEN("diagonal implemented as a lower triangle") {
    Eigen::VectorXcd d = random_diagonal(32, 0);
    DiagonalBox diag(d, false);
    std::shared_ptr<Circuit> c = diag.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::DiagonalMatrix<Complex, Eigen::Dynamic> D(d);
    Eigen::MatrixXcd M = D.derived();
    REQUIRE((U - M).cwiseAbs().sum() < ERR_EPS);
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 5);
    for (unsigned i = 0; i < 5; i++) {
      unit_vector_t args;
      for (unsigned j = i + 1; j < 5; j++) {
        args.push_back(Qubit(j));
      }
      args.push_back(Qubit(i));
      REQUIRE(cmds[i].get_args() == args);
    }
  }
  GIVEN("Non-unitary diagonal") {
    Eigen::Vector2cd diag(2. * i_, 1);
    REQUIRE_THROWS_MATCHES(
        DiagonalBox(diag), std::invalid_argument,
        MessageContains("not unitary"));
    Eigen::Vector2cd diag2(0, 1);
    REQUIRE_THROWS_MATCHES(
        DiagonalBox(diag2), std::invalid_argument,
        MessageContains("not unitary"));
  }
  GIVEN("vectors with wrong size") {
    Eigen::VectorXcd diag(1);
    diag << 1;
    REQUIRE_THROWS_MATCHES(
        DiagonalBox(diag), std::invalid_argument,
        MessageContains("not a power of 2"));
    Eigen::Vector3cd diag2(1, 0, 0);
    REQUIRE_THROWS_MATCHES(
        DiagonalBox(diag2), std::invalid_argument,
        MessageContains("not a power of 2"));
  }
  GIVEN("test dagger") {
    Eigen::Vector2cd diag(i_, 1);
    Eigen::Vector2cd dag_diag(-1. * i_, 1);
    DiagonalBox diagbox(diag);
    const DiagonalBox dag_box =
        static_cast<const DiagonalBox &>(*diagbox.dagger());
    REQUIRE((dag_diag - dag_box.get_diagonal()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("copy constructor") {
    Eigen::Vector2cd diag(i_, -1);
    DiagonalBox diagbox(diag, false);
    DiagonalBox diagbox2(diagbox);
    REQUIRE((diag - diagbox2.get_diagonal()).cwiseAbs().sum() < ERR_EPS);
    REQUIRE(!diagbox2.is_upper_triangle());
  }
}

}  // namespace test_DiagonalBox
}  // namespace tket
