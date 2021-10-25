#include "TestStatsStructs.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <sstream>

;

namespace tket {
namespace tsa_internal {
namespace tests {

void MinMaxAv::add(size_t result) {
  min = std::min(min, result);
  max = std::max(max, result);
  total += result;
}

void PartialTsaStatistics::add_problem_result(
    size_t initial_L, size_t final_L, size_t tokens, size_t swaps) {
  REQUIRE(final_L <= initial_L);
  REQUIRE(final_L + 2 * swaps >= initial_L);
  total_number_of_tokens += tokens;
  if (initial_L == 0) {
    CHECK(swaps == 0);
    l_decrease_percentages.add(100);
    powers.add(100);
    return;
  }
  ++number_of_problems;
  total_of_L += initial_L;
  const size_t l_decrease = initial_L - final_L;
  total_of_L_decreases += l_decrease;

  l_decrease_percentages.add((100 * (initial_L - final_L)) / initial_L);
  total_number_of_swaps += swaps;
  if (swaps == 0) {
    powers.add(0);
  } else {
    powers.add((50 * l_decrease) / swaps);
  }
}

std::string PartialTsaStatistics::str(size_t number_of_problems) const {
  REQUIRE(number_of_problems != 0);
  std::stringstream ss;
  ss << total_number_of_tokens << " tokens; " << total_of_L << " total L; "
     << total_number_of_swaps << " swaps.\nL-decr %: min "
     << l_decrease_percentages.min << ", max " << l_decrease_percentages.max
     << ", av " << l_decrease_percentages.total / number_of_problems
     << ".\nPower %: min " << powers.min << ", max " << powers.max << ", av "
     << powers.total / number_of_problems;
  return ss.str();
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
