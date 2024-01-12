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
#include <random>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

template <typename T>
static void check_product(T x, T y, bool expect_no_overflow = true) {
  const auto product = get_checked_product(x, y);
  INFO(
      "x=" << x << ", y=" << y << ", expect_no_overflow? "
           << expect_no_overflow);
  if (product) {
    CHECK(
        std::uintmax_t(product.value()) ==
        std::uintmax_t(x) * std::uintmax_t(y));
    CHECK(product.value() == x * y);
    CHECK(expect_no_overflow);
  } else {
    CHECK(!expect_no_overflow);
  }
}

template <typename T>
static void test_checked_sum_and_product() {
  const T max = std::numeric_limits<T>::max();
  const std::vector<T> values{0, 1, 2, 3, 4, 5, 10, 20, 50, 100, max / 3};
  for (T x : values) {
    for (T y : values) {
      const auto sum = get_checked_sum(x, y);
      REQUIRE(sum);
      REQUIRE(sum.value() == x + y);
      const auto product = get_checked_product(x, y);
      const T min_v = std::min(x, y);
      const T max_v = std::max(x, y);

      // We happen to know that max = 2^n-1 for n=16,32,64.
      // Mod 3, this is (-1)^2k - 1 = 0 for k, so "max/3"
      // IS actually the value max/3 exactly.
      const bool expect_overflow = (max_v == values.back() && min_v > 3);
      check_product(x, y, !expect_overflow);
    }
  }
  // max/3 is an integer. 3*(max/3) == max is fine, 3*(max/3 + 1) is not.
  check_product<T>(max / 3, 3);
  check_product<T>(max / 4, 4);

  check_product<T>(1 + max / 3, 3, false);
  // max = 4k-1, so   max DIV 4 = k-1, so (max DIV 4 + 1)*4 = k overflows.
  check_product<T>(1 + max / 4, 4, false);

  // max DIV 2 = 2k-1, so if n = max DIV 2 then 2*n = max-1 is fine.
  const T max_div_2 = max / 2;
  check_product<T>(max_div_2, 2);
  CHECK(get_checked_sum(max_div_2, max_div_2));
  for (T ii = 0; ii < 10; ++ii) {
    const T x = max_div_2 - ii;
    const T y = max_div_2 + ii;
    CHECK(get_checked_sum<T>(x, y + 1));
    for (T jj = 2; jj < 10; ++jj) {
      CHECK(!get_checked_sum<T>(x, y + jj));
    }
  }

  const T x(0.5 * std::sqrt(max));

  // x ~ sqrt(M)/2, so x^2 ~ M/4 does not overflow.
  CHECK(x * x > max / 5);
  check_product(x, x);
  const T y = 5 * x;

  // The true value ~5M/4 overflows.
  check_product(x, y, false);

  // ~5M/4 will be reduced to ~ M/4 (mod M).
  // Compilers might give a warning here; but the overflow is intentional!
  // e.g., "warning: unsigned conversion from 'int' to 'short unsigned int'
  // changes value from '80645' to '15109' [-Woverflow]"
  const T xy_overflow = static_cast<T>(x * y);

  // ~((M/4)/3) * 12
  check_product<T>(xy_overflow / 3, 12);
}

SCENARIO("Test sum and product with checked overflows") {
  test_checked_sum_and_product<unsigned>();
  test_checked_sum_and_product<std::size_t>();
  test_checked_sum_and_product<std::uint16_t>();
  test_checked_sum_and_product<std::uint32_t>();
  test_checked_sum_and_product<std::uint64_t>();
  test_checked_sum_and_product<std::uintmax_t>();
}

namespace {
struct UIntMaxResult {
  unsigned sum_overflow_count = 0;
  unsigned sum_normal_count = 0;
  unsigned product_overflow_count = 0;
  unsigned product_normal_count = 0;
};
}  // namespace

// Since uintmax is at least 64 bit,
// we can use it to test 16-bit and 32-bit ints.
template <typename T>
static UIntMaxResult test_checked_sum_and_product_with_uintmax() {
  const std::uintmax_t max_value = std::numeric_limits<T>::max();
  REQUIRE(get_checked_product(std::uintmax_t(max_value + 1), max_value));
  std::vector<std::uintmax_t> numbers{2, 3};
  for (;;) {
    auto next = numbers.back();
    next *= 3;
    next /= 2;
    if (next * next >= 10 * max_value) {
      break;
    }
    numbers.push_back(next);
  }
  const auto old_size = numbers.size();
  for (unsigned ii = 0; ii < old_size; ++ii) {
    numbers.push_back(max_value / numbers[ii]);
  }
  // We have a whole load of numbers ranging in size from 2 to MAX/2.
  {
    const auto [min, max] =
        std::minmax_element(numbers.cbegin(), numbers.cend());
    CHECK(*min >= 2);
    CHECK(*max <= max_value / 2);
  }

  // Just multiply all the pairs...
  UIntMaxResult result;
  for (auto x : numbers) {
    for (auto y : numbers) {
      // x,y have been chosen so that x*y can be a bit bigger than MAX,
      // but x,y < ~3 . sqrt(MAX) always.
      const bool expect_no_overflow = x * y <= max_value;
      check_product(T(x), T(y), expect_no_overflow);
      if (expect_no_overflow) {
        ++result.product_normal_count;
      } else {
        ++result.product_overflow_count;
      }
    }
  }
  // Now do some addition.
  const std::intmax_t max_value_signed = max_value;
  for (std::intmax_t ii = -10; ii < 10; ++ii) {
    for (std::intmax_t x : numbers) {
      std::intmax_t y = max_value_signed + ii - x;
      CHECK(y >= 0);
      if (y < 0 || y > max_value_signed) {
        continue;
      }
      const auto t_sum_opt = get_checked_sum(T(x), T(y));
      const std::intmax_t actual_sum = x + y;
      const bool normal = (actual_sum <= max_value_signed);
      if (t_sum_opt) {
        CHECK(normal);
        CHECK(t_sum_opt.value() == actual_sum);
        ++result.sum_normal_count;
      } else {
        CHECK(!normal);
        ++result.sum_overflow_count;
      }
    }
  }
  return result;
}

SCENARIO("Use uintmax to test sum, product for smaller int sizes") {
  const auto result_16_bits =
      test_checked_sum_and_product_with_uintmax<std::uint16_t>();
  CHECK(result_16_bits.sum_normal_count == 352);
  CHECK(result_16_bits.sum_overflow_count == 267);
  CHECK(result_16_bits.product_normal_count == 528);
  CHECK(result_16_bits.product_overflow_count == 496);

  const auto result_32_bits =
      test_checked_sum_and_product_with_uintmax<std::uint32_t>();
  CHECK(result_32_bits.sum_normal_count == 638);
  CHECK(result_32_bits.sum_overflow_count == 501);
  CHECK(result_32_bits.product_normal_count == 1711);
  CHECK(result_32_bits.product_overflow_count == 1653);
}

SCENARIO("test get_sum (or product) _or_throw") {
  const std::vector<std::uint16_t> numbers{0,    1,     2,     12,
                                           4124, 12313, 51235, 65535};
  std::stringstream errors;

  for (std::uint16_t number1 : numbers) {
    for (std::uint16_t number2 : numbers) {
      const std::uint64_t num1_64bit = number1;
      const std::uint64_t num2_64bit = number2;
      const std::uint64_t sum_64bit = num1_64bit + num2_64bit;
      const std::uint64_t product_64bit = num1_64bit * num2_64bit;
      const bool sum_valid = sum_64bit <= 65535;
      const bool product_valid = product_64bit <= 65535;

      try {
        const std::uint16_t sum = get_sum_or_throw(number1, number2);
        REQUIRE(static_cast<std::uint64_t>(sum) == sum_64bit);
        REQUIRE(sum_valid);
        const std::uint16_t product = get_product_or_throw(number1, number2);
        REQUIRE(static_cast<std::uint64_t>(product) == product_64bit);
        REQUIRE(product_valid);
      } catch (const IntegerOverflow& e) {
        errors << e.what() << " ";
      }
    }
  }
  // clang-format off
  CHECK(errors.str() == "(1 + 65535) (2 * 51235) (2 + 65535) (12 * 12313)"
    " (12 * 51235) (12 + 65535) (4124 * 4124) (4124 * 12313) (4124 * 51235)"
    " (4124 + 65535) (12313 * 12) (12313 * 4124) (12313 * 12313)"
    " (12313 * 51235) (12313 + 65535) (51235 * 2) (51235 * 12)"
    " (51235 * 4124) (51235 * 12313) (51235 + 51235) (51235 + 65535)"
    " (65535 + 1) (65535 + 2) (65535 + 12) (65535 + 4124) (65535 + 12313)"
    " (65535 + 51235) (65535 + 65535) ");
  // clang-format on
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
