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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <random>
#include <tkwsm/Common/DyadicFraction.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

namespace {
// What do we multiply by?
enum class Action { NORMAL = 0, LARGE, NORMAL_FRACTION, SMALL_FRACTION };

const double small_value = 1e-30;
const double large_value = 1.0 / small_value;
const double min_value = small_value * 1e-20;
const double max_value = 1.0 / min_value;

static void get_next_action(double last_value, Action& action) {
  if (last_value < small_value) {
    action = Action::LARGE;
    return;
  }
  if (last_value > large_value) {
    action = Action::SMALL_FRACTION;
  }
}

const DyadicFraction::UInt min_normal_value = 1;
const DyadicFraction::UInt min_large_value =
    std::numeric_limits<std::uint32_t>::max();

// Returns TRUE if we should multiply by an (n over K) fraction
// rather than an int.
// CHange the random number to something suitable to perform
// the given action.
static bool convert_random_number(Action action, std::uint64_t& rnd_number) {
  switch (action) {
    case Action::NORMAL:
      rnd_number >>= 40;
      rnd_number = std::min(rnd_number, min_normal_value);
      return false;

    case Action::LARGE:
      rnd_number = std::min(rnd_number, min_large_value);
      return false;

    case Action::NORMAL_FRACTION:
      rnd_number &= 0xffff;
      rnd_number = std::min(rnd_number, min_normal_value);
      return true;

    case Action::SMALL_FRACTION:
      rnd_number &= 0xf;
      ++rnd_number;
      return true;

    default:
      REQUIRE(false);
  }
  return true;
}

struct LtOperatorCounts {
  unsigned successes = 0;
  unsigned failures = 0;
  unsigned inconclusive = 0;

  // Test all-against-all with lt.
  explicit LtOperatorCounts(
      const std::vector<DyadicFraction>& nonzero_fractions, double epsilon) {
    for (const auto& frac1 : nonzero_fractions) {
      for (const auto& frac2 : nonzero_fractions) {
        const auto approx_val1 = frac1.get_double();
        const auto approx_val2 = frac2.get_double();
        // Are the results too close to be reliable, due to roundoff?
        const double ratio = approx_val1 / approx_val2;
        const double diff = std::abs(ratio - 1.0);
        if (diff < epsilon) {
          ++inconclusive;
          continue;
        }
        const bool lt_result = frac1 < frac2;
        const bool double_lt_result = approx_val1 < approx_val2;
        if (lt_result == double_lt_result) {
          ++successes;
        } else {
          ++failures;
        }
      }
    }
  }
};

struct MultiplicationCounts {
  unsigned successes = 0;
  unsigned failures = 0;

  // Test all-against-all.
  explicit MultiplicationCounts(
      const std::vector<DyadicFraction>& nonzero_fractions, double epsilon) {
    for (const auto& frac1 : nonzero_fractions) {
      for (const auto& frac2 : nonzero_fractions) {
        const auto approx_val1 = frac1.get_double();
        const auto approx_val2 = frac2.get_double();
        const auto approx_product = approx_val1 * approx_val2;
        REQUIRE(approx_product > 0.1 * min_value * min_value);
        REQUIRE(approx_product < 10.0 * max_value * max_value);
        auto frac1_copy = frac1;
        const auto approx_product_again = frac1_copy.mult(frac2).get_double();
        const double diff = std::abs(approx_product - approx_product_again);
        if (diff < epsilon * approx_product) {
          ++successes;
        } else {
          ++failures;
        }
      }
    }
  }
};

}  // namespace

static void test_trivial_multipliers(
    const std::vector<DyadicFraction>& nonzero_fractions) {
  const DyadicFraction zero;
  REQUIRE(zero.get_double() == 0.0);
  {
    auto zero_again = zero;
    zero_again.mult(12345);
    REQUIRE(zero_again.get_double() == 0.0);
  }
  {
    auto zero_again = zero;
    zero_again.mult_n_over_k(98765);
    REQUIRE(zero_again.get_double() == 0.0);
  }
  const DyadicFraction one(1);

  for (const auto& frac : nonzero_fractions) {
    const double approx_value = frac.get_double();
    {
      auto frac1 = frac;
      // Multiplication by 1 really should be the identity operation,
      // so doubles match EXACTLY.
      REQUIRE(frac1.mult(one).get_double() == approx_value);
      REQUIRE(frac1.mult_n_over_k(1024).get_double() == approx_value);
      REQUIRE(frac1.mult(zero).get_double() == 0.0);
      REQUIRE(frac1.mult(frac).get_double() == 0.0);
      REQUIRE(frac1.mult_n_over_k(100).get_double() == 0.0);
    }
    // Multiplying by powers of 2 is exactly reversible.
    auto frac1 = frac;
    unsigned power_of_two = 1;
    for (int nn = 0; nn < 5; ++nn) {
      power_of_two *= 2;
      REQUIRE(
          frac1.mult(power_of_two)
              .mult_n_over_k(1024 / power_of_two)
              .get_double() == approx_value);
    }
  }
}

SCENARIO("Create random DyadicFractions and check approximation with doubles") {
  const unsigned number_of_fractions = 100;
  const double diff_epsilon = 1e-9;
  std::vector<DyadicFraction> nonzero_fractions;
  nonzero_fractions.reserve(number_of_fractions);
  std::mt19937_64 r_engine;

  nonzero_fractions.emplace_back(DyadicFraction(1));
  // In this special case, we know it must be exact.
  REQUIRE(nonzero_fractions.back().get_double() == 1.0);
  double last_value = 1.0;

  const std::array<Action, 4> actions{
      Action::NORMAL, Action::LARGE, Action::NORMAL_FRACTION,
      Action::SMALL_FRACTION};

  while (nonzero_fractions.size() < number_of_fractions) {
    auto last_fraction = nonzero_fractions.back();
    auto rnd_number = r_engine();
    auto action = actions[rnd_number % actions.size()];

    // Ensure that our fraction doesn't wander off too much into stupidly
    // large or small territory.
    get_next_action(last_value, action);
    const bool use_fraction = convert_random_number(action, rnd_number);
    last_value *= static_cast<double>(rnd_number);

    if (use_fraction) {
      last_fraction.mult_n_over_k(rnd_number);
      last_value /= 1024.0;
    } else {
      last_fraction.mult(rnd_number);
    }
    REQUIRE(last_value > min_value);
    REQUIRE(last_value < max_value);
    const double approx_value = last_fraction.get_double();
    const double diff = std::abs(approx_value - last_value);

    CHECK(diff < last_value * diff_epsilon);

    nonzero_fractions.emplace_back(last_fraction);
    // Stop it slowly drifting away...
    last_value = approx_value;
  }
  {
    const LtOperatorCounts counts(nonzero_fractions, 1e-10);
    CHECK(counts.successes == 9848);
    CHECK(counts.failures == 0);
    CHECK(counts.inconclusive == 152);
  }
  {
    const MultiplicationCounts counts(nonzero_fractions, 1e-8);
    CHECK(counts.successes == 10000);
    CHECK(counts.failures == 0);
  }
  test_trivial_multipliers(nonzero_fractions);
}

SCENARIO("Large random products") {
  // Multiply by many random ints, and fractions;
  // then check the logs.
  unsigned remaining_mults = 200;
  unsigned remaining_pk_fracs = 50;

  DyadicFraction fraction(1);
  double logs_sum = 0.0;

  std::mt19937_64 r_engine;
  auto bits = r_engine();

  // Get an int between 0 and 2^k-1, by extracting the given number of bits.
  const auto get_bits = [&bits, &r_engine](unsigned num_bits) -> std::uint64_t {
    // OK, this introduces some bias, but who cares.
    if (bits == 0) {
      bits = r_engine();
    }
    std::uint64_t mask = 1;
    mask <<= num_bits;
    --mask;
    auto x = bits & mask;
    bits >>= num_bits;
    return x;
  };

  while (remaining_mults > 0 || remaining_pk_fracs > 0) {
    if (remaining_mults > 0) {
      --remaining_mults;
      auto x = get_bits(10);
      x += 2;
      const auto log_val = std::log(x);
      logs_sum += log_val;
      fraction.mult(x);
    }
    if (remaining_pk_fracs > 0) {
      --remaining_pk_fracs;
      auto x = get_bits(4);
      ++x;
      const auto log_val = std::log(x / 1024.0);
      logs_sum += log_val;
      fraction.mult_n_over_k(x);
    }
  }
  const double recalc_log = fraction.get_log();
  const double diff = std::abs(logs_sum - recalc_log);

  // Notice that taking the exponential would give
  // a number x too large for doubles to represent.
  CHECK(std::abs(logs_sum - 829.184) < 0.01);
  CHECK(diff < 1e-10);
}

SCENARIO("Log of large factorial") {
  // Find N! for large N.
  DyadicFraction fraction(1);
  double calc_log = 0.0;
  for (unsigned ii = 2; ii <= 1000; ++ii) {
    fraction.mult(ii);
    calc_log += std::log(ii);
  }
  const double factorial_log_approx = fraction.get_log();
  const double diff = std::abs(factorial_log_approx - calc_log);
  CHECK(std::abs(calc_log - 5912.13) < 0.01);
  CHECK(diff < 1e-10);
  CHECK(fraction.str() == "val=6076743920982394875 exp=8467");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
