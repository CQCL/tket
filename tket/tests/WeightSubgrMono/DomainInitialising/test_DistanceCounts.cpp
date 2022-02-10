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

#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <sstream>

#include "WeightSubgrMono/DomainInitialising/DistanceCounts.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// Vector entries are 1,2,3,4,5,6,7,8 so use 3 bits to represent.
static std::vector<std::size_t> get_distance_counts(
    std::uint64_t& encoding, unsigned size) {
  std::vector<std::size_t> result;
  while (result.size() < size) {
    auto value = encoding & 7;
    result.push_back(value + 1);
    encoding >>= 3;
  }
  return result;
}

// Returns true if reduced to empty.
static bool remove_top_zeros(std::vector<std::size_t>& counts) {
  for (;;) {
    if (counts.empty()) {
      return true;
    }
    if (counts.back() != 0) {
      break;
    }
    counts.pop_back();
  }
  return false;
}

// Reimplement DistanceCounts::test_against_target
static bool slow_pair_up_entries(
    std::vector<std::size_t>& p_counts, std::vector<std::size_t>& t_counts) {
  // Recursively destroy...
  if (remove_top_zeros(p_counts)) {
    return true;
  }
  if (remove_top_zeros(t_counts)) {
    return false;
  }
  // Both p,t are nonempty and have no top zeros.
  while (p_counts.size() < t_counts.size()) {
    // Higher level t, which cannot be paired off, are irrelevant.
    t_counts.pop_back();
  }
  REQUIRE(p_counts.size() >= t_counts.size());
  REQUIRE(!t_counts.empty());

  // Match the top elements.
  if (p_counts.back() >= t_counts.back()) {
    p_counts.back() -= t_counts.back();
  } else {
    t_counts.back() -= p_counts.back();
  }
  return slow_pair_up_entries(p_counts, t_counts);
}

SCENARIO("Exhaustive distance counts reductions") {
  std::vector<std::vector<std::size_t>> counts_list;
  for (unsigned size = 0; size <= 6; ++size) {
    for (std::uint64_t encoding = 0; encoding < 1000000; ++encoding) {
      auto encoding_copy = encoding;
      counts_list.emplace_back(get_distance_counts(encoding_copy, size));
      counts_list.emplace_back(counts_list.back());
      counts_list.back().push_back(0);
      if (encoding_copy > 0) {
        // Still some bits left, so we've calculated enough.
        break;
      }
      // Skip a few, for larger sizes.
      if (size >= 3) {
        encoding += 11;
      }
      if (size >= 4) {
        encoding += 13 * 7;
      }
      if (size >= 5) {
        encoding += 17 * 23;
      }
      if (size >= 6) {
        encoding += 19 * 71 * 61;
      }
    }
  }
  CHECK(counts_list.size() == 468);
  unsigned total_entries = 0;
  for (const auto& counts : counts_list) {
    total_entries += counts.size();
  }
  CHECK(total_entries == 1844);
  unsigned total_returned_true = 0;
  unsigned total_returned_false = 0;
  // Now calculate in 2 different ways.
  std::vector<std::size_t> p_counts_copy;
  std::vector<std::size_t> t_counts_copy;
  for (const auto& p_counts : counts_list) {
    for (const auto& t_counts : counts_list) {
      p_counts_copy = p_counts;
      t_counts_copy = t_counts;
      const bool success = slow_pair_up_entries(p_counts_copy, t_counts_copy);
      if (success) {
        ++total_returned_true;
      } else {
        ++total_returned_false;
      }
      CHECK(
          success ==
          DistanceCounts::test_against_target(p_counts_copy, t_counts_copy));
    }
  }
  CHECK(total_returned_true == 217168);
  CHECK(total_returned_false == 1856);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
