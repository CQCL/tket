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

#include <catch2/catch_test_macros.hpp>
#include <random>
#include <tkwsm/Common/BitFunctions.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

static void test_bitlength(std::uint64_t x) {
  const auto bit_length = BitFunctions::get_bit_length(x);
  if (x == 0) {
    REQUIRE(bit_length == 0);
    return;
  }
  REQUIRE(bit_length > 0);
  REQUIRE(bit_length <= 64);
  const auto x_shifted = (x >> (bit_length - 1));
  REQUIRE(x_shifted == 1);
}

static void test_trailing_zeros(std::uint64_t x) {
  const auto zeros = BitFunctions::get_number_of_rightmost_zero_bits(x);
  if (x == 0) {
    REQUIRE(zeros == 64);
    return;
  }
  REQUIRE(zeros < 64);
  const auto x_shifted = x >> zeros;
  const auto x_again = x_shifted << zeros;
  REQUIRE(x_again == x);
}

SCENARIO(
    "Test get_bit_length and get_number_of_rightmost_zero_bits on random "
    "bits") {
  std::mt19937_64 r_engine;
  test_bitlength(0);
  test_trailing_zeros(0);
  for (int nn = 0; nn < 100; ++nn) {
    auto x = r_engine();
    do {
      test_bitlength(x);
      test_trailing_zeros(x);
      x >>= 1;
    } while (x != 0);
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
