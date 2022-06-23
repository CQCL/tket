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
#include <numeric>
#include <random>
#include <tkrng/RNG.hpp>

#include "WeightSubgrMono/Common/SetIntersection.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("Test set intersection with unsigned ints as bitsets") {
  const unsigned number_of_bits = 6;

  // Element[i] gives the representation of i as a set.
  std::vector<std::set<unsigned>> raw_sets(1u << number_of_bits);

  std::vector<std::vector<unsigned>> raw_vectors(raw_sets.size());
  std::vector<std::vector<std::pair<unsigned, unsigned>>> raw_vectors_with_junk(
      raw_sets.size());
  for (unsigned ii = 0; ii < raw_sets.size(); ++ii) {
    unsigned ii_copy = ii;
    unsigned element = 0;
    while (ii_copy != 0) {
      if ((ii_copy & 1u) != 0) {
        raw_sets[ii].insert(element);
      }
      ++element;
      ii_copy >>= 1;
    }
    raw_vectors[ii] = {raw_sets[ii].cbegin(), raw_sets[ii].cend()};

    raw_vectors_with_junk[ii].resize(raw_vectors[ii].size());
    for (unsigned jj = 0; jj < raw_vectors[ii].size(); ++jj) {
      raw_vectors_with_junk[ii][jj].first = raw_vectors[ii][jj];
      raw_vectors_with_junk[ii][jj].second = jj % 17;
    }
  }
  std::set<unsigned> calc_set;
  std::vector<unsigned> calc_vect;
  // Fill calc set with junk each time.
  const std::set<unsigned> junk_set{0, 2, 5, 7, 8, 2342, 56235};
  const std::vector<unsigned> junk_vect{123, 0, 34, 1, 34, 2, 2, 42, 5435};
  // Test all against all.
  // Do them in random order...extra paranoia!
  std::vector<unsigned> bitsets(raw_sets.size());
  std::iota(bitsets.begin(), bitsets.end(), 0);
  {
    RNG rng;
    rng.do_shuffle(bitsets);
  }

  unsigned disjoint_count = 0;

  for (unsigned bitset1 : bitsets) {
    for (unsigned bitset2 : bitsets) {
      const auto& set1 = raw_sets[bitset1];
      const auto& set2 = raw_sets[bitset2];
      const auto& vect2 = raw_vectors[bitset2];
      const unsigned final_bitset = bitset1 & bitset2;
      const auto& final_set = raw_sets[final_bitset];
      const auto& final_vect = raw_vectors[final_bitset];
      calc_vect = junk_vect;
      fill_intersection(set1, vect2, calc_vect);
      REQUIRE(final_vect == calc_vect);

      const bool are_disjoint = disjoint(set1, set2);
      REQUIRE(are_disjoint == (final_bitset == 0));
      if (are_disjoint) {
        ++disjoint_count;
      }
      calc_set = junk_set;
      fill_intersection_ignoring_second_elements(
          set1, raw_vectors_with_junk[bitset2], calc_set);
      REQUIRE(final_set == calc_set);
    }
  }

  // A bitset with k ones is disjoint from 2^{n-k} other bitsets.
  // (2 choices for each of the n-k other bits).
  // Therefore, each of these (n choose k) choices contributes 2^{n-k}.
  // So, the total comes from a Binomial expansion!
  //  sum_{k=0,1,...,n} (n C k).2^{n-k} = 2^n . (1+1/2)^n = 3^n.
  unsigned expected_disjoint_count = 1;
  for (unsigned jj = 1; jj <= number_of_bits; ++jj) {
    expected_disjoint_count *= 3;
  }
  REQUIRE(disjoint_count == expected_disjoint_count);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
