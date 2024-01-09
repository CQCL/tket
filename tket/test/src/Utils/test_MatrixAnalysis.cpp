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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "../testutil.hpp"
#include "tket/Utils/MatrixAnalysis.hpp"

// Element[n] = 2^n
static std::vector<unsigned> get_powers_of_two() {
  std::vector<unsigned> result;
  result.reserve(31);
  result.push_back(1);
  for (;;) {
    // Unlike left shifting, multiplication by 2 overflow
    // with unsigned int types IS guaranteed by the C++ standard.
    const unsigned new_power = result.back() * 2;
    if (result.back() != new_power / 2) {
      return result;
    }
    result.push_back(new_power);
  }
}

namespace tket {
namespace test_MatrixAnalysis {

SCENARIO("Powers of two") {
  const auto powers_of_two = get_powers_of_two();
  GIVEN("Manually calculated powers of two") {
    for (unsigned nn = 0; nn < powers_of_two.size(); ++nn) {
      CHECK(get_matrix_size(nn) == powers_of_two[nn]);
      CHECK(get_number_of_qubits(powers_of_two[nn]) == nn);
    }
  }
  GIVEN("log2 for small numbers not powers of two") {
    for (unsigned mm = 0; mm < 1000; ++mm) {
      if (!std::binary_search(
              powers_of_two.cbegin(), powers_of_two.cend(), mm)) {
        CHECK_THROWS(get_number_of_qubits(mm));
      }
    }
  }
  GIVEN("log2 for numbers close to, but not equal to, powers of two") {
    for (std::intmax_t power_of_two : powers_of_two) {
      if (power_of_two < 1000) {
        continue;
      }
      for (std::intmax_t dx = -50; dx < 50; ++dx) {
        if (dx == 0) {
          continue;
        }
        CHECK_THROWS(get_number_of_qubits(power_of_two + dx));
      }
    }
  }
  GIVEN("log2 for large numbers very close to the limit") {
    unsigned not_power_of_two = std::numeric_limits<unsigned>::max();
    for (int count = 0; count < 50; ++count) {
      CHECK_THROWS(get_number_of_qubits(not_power_of_two));
      --not_power_of_two;
    }
  }
  GIVEN("2^n for large n, should overflow.") {
    const unsigned min_qubits = std::numeric_limits<unsigned>::digits;
    const unsigned max_qubits = min_qubits + 1000;
    for (unsigned too_many_qubits = min_qubits; too_many_qubits <= max_qubits;
         ++too_many_qubits) {
      CHECK_THROWS(get_matrix_size(too_many_qubits));
    }
  }
}

SCENARIO("Test nth root of a 2x2 unitary") {
  GIVEN("Identity") {
    Eigen::Matrix2cd m = Eigen::Matrix2cd::Identity();
    Eigen::Matrix2cd root1 = nth_root(m, 1);
    REQUIRE(root1.isApprox(Eigen::Matrix2cd::Identity(), ERR_EPS));
    Eigen::Matrix2cd root0 = nth_root(m, 0);
    REQUIRE(root0.isApprox(Eigen::Matrix2cd::Identity(), ERR_EPS));
  }

  GIVEN("Random unitary 0th root") {
    Eigen::Matrix2cd m = random_unitary(2, 1);
    REQUIRE_THROWS(nth_root(m, 0));
  }

  GIVEN("Random unitary") {
    for (unsigned i = 0; i < 100; i++) {
      Eigen::Matrix2cd m = random_unitary(2, i);
      for (unsigned n = 1; n < 10; n++) {
        Eigen::Matrix2cd root = nth_root(m, n);
        Eigen::Matrix2cd u = Eigen::Matrix2cd::Identity();
        for (unsigned k = 0; k < n; k++) {
          u = root * u;
        }
        REQUIRE(u.isApprox(m, ERR_EPS));
      }
    }
  }
}

SCENARIO("Test unitary product") {
  GIVEN("Two perturbed unitaries") {
    Eigen::Matrix4cd U = random_unitary(4, 0);
    U(1, 2) += 0.1;
    Eigen::Matrix4cd V = random_unitary(4, 1);
    V(3, 0) -= 0.1;
    Eigen::Matrix4cd UV = U * V;
    Eigen::Matrix4cd X = unitary_product2(U, V);
    REQUIRE(!is_unitary(UV));
    REQUIRE(is_unitary(X));
    double d = (X - UV).squaredNorm();
    for (int i = 0; i < 10; i++) {
      Eigen::Matrix4cd Y = random_unitary(4, 2 + i);
      CHECK((Y - UV).squaredNorm() >= d);
    }
  }
  GIVEN("Three perturbed unitaries") {
    Eigen::Matrix4cd U = random_unitary(4, 12);
    U(1, 2) += 0.1;
    Eigen::Matrix4cd V = random_unitary(4, 13);
    V(3, 0) -= 0.1;
    Eigen::Matrix4cd W = random_unitary(4, 14);
    W(2, 2) += 0.1;
    Eigen::Matrix4cd UVW = U * V * W;
    Eigen::Matrix4cd X = unitary_product3(U, V, W);
    REQUIRE(!is_unitary(UVW));
    REQUIRE(is_unitary(X));
    double d = (X - UVW).squaredNorm();
    for (int i = 0; i < 10; i++) {
      Eigen::Matrix4cd Y = random_unitary(4, 15 + i);
      CHECK((Y - UVW).squaredNorm() >= d);
    }
  }
}

}  // namespace test_MatrixAnalysis
}  // namespace tket
