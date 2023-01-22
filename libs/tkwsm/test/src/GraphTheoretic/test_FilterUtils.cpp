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
#include <cstddef>
#include <map>
#include <sstream>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <tkwsm/GraphTheoretic/FilterUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

static std::string to_string(const FilterUtils::DegreeCounts& deg_counts) {
  std::stringstream ss;
  ss << "[ ";
  for (const auto& entry : deg_counts) {
    ss << entry.first << ":" << entry.second << " ";
  }
  ss << "]";
  return ss.str();
}

// We have two different implementations of degree sequence compatibility
// (because it's convenient, in different applications, for the input data
// to come in different formats and we don't want to waste time converting
// between them).
// We test that they agree with each other; both are used extensively
// for solving problems, so if either one had an error,
// we might expect it to show up in at least one specific problem.
SCENARIO("Test random degree sequences for compatibility") {
  const unsigned list_size = 100;
  const unsigned divisor = 30;
  const unsigned max_degree_minus_1 = 6;
  const unsigned max_count_minus_1 = 4;

  RNG rng;
  std::vector<FilterUtils::DegreeCounts> degree_counts_list(list_size);
  std::vector<std::vector<std::size_t>> raw_deg_seqs(list_size);

  std::map<unsigned, unsigned> degree_counts_map;

  for (unsigned ii = 0; ii < degree_counts_list.size(); ++ii) {
    auto& degree_counts = degree_counts_list[ii];
    degree_counts_map.clear();
    for (unsigned number_of_degrees = 1 + (ii / divisor); number_of_degrees > 0;
         --number_of_degrees) {
      size_t a0 = rng.get_size_t(max_degree_minus_1);
      size_t a1 = rng.get_size_t(max_count_minus_1);
      degree_counts_map[1 + a0] += 1 + a1;
    }
    unsigned size = 0;
    for (const auto& entry : degree_counts_map) {
      degree_counts.emplace_back(entry.first, entry.second);
      size += entry.second;
    }
    auto& raw_deg_seq = raw_deg_seqs[ii];
    raw_deg_seq.reserve(size);
    for (const auto& entry : degree_counts) {
      for (unsigned ii = 0; ii < entry.second; ++ii) {
        raw_deg_seq.push_back(entry.first);
      }
    }
    REQUIRE(std::is_sorted(raw_deg_seq.cbegin(), raw_deg_seq.cend()));
    REQUIRE(raw_deg_seq.size() == size);
  }
  for (unsigned ii = 0; ii < list_size; ++ii) {
    for (unsigned jj = 0; jj < list_size; ++jj) {
      const bool compat_with_vect =
          FilterUtils::compatible_sorted_degree_sequences(
              raw_deg_seqs[ii], raw_deg_seqs[jj]);

      if (ii == jj) {
        REQUIRE(compat_with_vect);
      }
      REQUIRE(
          compat_with_vect ==
          FilterUtils::compatible_sorted_degree_counts(
              degree_counts_list[ii], degree_counts_list[jj]));
    }
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
