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
#include <limits>
#include <map>

#include "TokenSwapping/FilteredSwapSequences.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

SCENARIO("Trivial table lookup tests") {
  // Permutation hash 0 is the identity.
  for (unsigned edges_bitset = 0; edges_bitset < 50; ++edges_bitset) {
    const FilteredSwapSequences::SingleSequenceData identity_result(
        0, edges_bitset, 10);
    CHECK(identity_result.edges_bitset == 0);
    CHECK(identity_result.swaps_code == 0);
    CHECK(identity_result.number_of_swaps == 0);
  }

  // (0,1) is the first swap (index 0). So, just need to include that bit.
  for (unsigned edges_bitset = 1; edges_bitset < 50; edges_bitset += 2) {
    const FilteredSwapSequences::SingleSequenceData single_swap_result(
        2, edges_bitset, 10);
    CHECK(single_swap_result.edges_bitset == 0x1);
    CHECK(single_swap_result.swaps_code == 0x1);
    CHECK(single_swap_result.number_of_swaps == 1);
  }

  // Enact a non-identity permutation without edges; impossible!
  const vector<unsigned> nontrivial_permutation_hashes{2,  3,  4,  5,  6,
                                                       22, 33, 32, 42, 222};
  for (unsigned perm_hash : nontrivial_permutation_hashes) {
    const FilteredSwapSequences::SingleSequenceData impossible_result(
        perm_hash, 0x0, 10);
    CHECK(impossible_result.edges_bitset == 0);
    CHECK(impossible_result.swaps_code == 0);
    CHECK(
        impossible_result.number_of_swaps ==
        std::numeric_limits<unsigned>::max());
  }
}

SCENARIO("Random entries test") {
  // Note: the entries are definitely NOT real swap sequence codes,
  // they are just random nunmbers.

  const unsigned num_bits = 15;

  std::map<SwapConversion::SwapHash, FilteredSwapSequences::SingleSequenceData>
      original_entries;
  // Make a vector, with duplicates.
  vector<SwapConversion::SwapHash> codes_vect;

  RNG rng;

  for (unsigned nn = 0; nn < 1000; ++nn) {
    const auto num_swaps = rng.get_size_t(1, 6);
    SwapConversion::SwapHash code = 0;
    SwapConversion::EdgesBitset edges_bitset = 0;

    for (unsigned mm = 0; mm < num_swaps; ++mm) {
      const auto new_swap = rng.get_size_t(1, num_bits);
      code <<= 4;
      code |= new_swap;
      edges_bitset |= (1u << (new_swap - 1));
    }
    auto& entry = original_entries[code];
    entry.edges_bitset = edges_bitset;
    entry.swaps_code = code;
    entry.number_of_swaps = num_swaps;
    for (int kk = 0; kk < 3; ++kk) {
      codes_vect.push_back(code);
    }
  }
  rng.do_shuffle(codes_vect);

  FilteredSwapSequences filtered_sequences;
  REQUIRE(filtered_sequences.get_total_number_of_entries() == 0);
  filtered_sequences.initialise(codes_vect);
  REQUIRE(
      filtered_sequences.get_total_number_of_entries() ==
      original_entries.size());

  // Now, look up every single edge bitset in turn and check that it finds the
  // (joint) fewest number of swaps.
  const SwapConversion::EdgesBitset max_bitset = (1u << num_bits) - 1;
  for (SwapConversion::EdgesBitset bitset = 0; bitset <= max_bitset; ++bitset) {
    // By brute force, find the (joint) fewest number of swaps in a sequence
    // using only this bitset.
    SwapConversion::SwapHash fewest_swaps_code =
        std::numeric_limits<SwapConversion::SwapHash>::max();
    unsigned number_of_swaps = 10000;
    for (const auto& entry : original_entries) {
      if (entry.first > fewest_swaps_code) {
        break;
      }
      REQUIRE(entry.second.number_of_swaps <= number_of_swaps);
      // Is it a subset?
      if ((entry.second.edges_bitset & bitset) != entry.second.edges_bitset) {
        continue;
      }
      // We've found a better entry than what we've got.
      number_of_swaps = entry.second.number_of_swaps;
      fewest_swaps_code = entry.first;
    }

    for (unsigned max_num_swaps = 1; max_num_swaps < num_bits + 3;
         ++max_num_swaps) {
      const auto result =
          filtered_sequences.get_lookup_result(bitset, max_num_swaps);
      if (result.number_of_swaps <= max_num_swaps) {
        // It found an entry. It must be an existing entry.
        const auto& existing_entry = original_entries.at(result.swaps_code);
        REQUIRE(result.number_of_swaps == existing_entry.number_of_swaps);
        REQUIRE(result.edges_bitset == existing_entry.edges_bitset);
        REQUIRE(result.swaps_code == existing_entry.swaps_code);

        // ...and it must be valid...
        REQUIRE((result.edges_bitset & bitset) == result.edges_bitset);
        REQUIRE(result.number_of_swaps == number_of_swaps);
      } else {
        // No entry was found. It MUST be because none actually exist, subject
        // to the constraints.
        REQUIRE(number_of_swaps > max_num_swaps);
        // Must be a null result.
        REQUIRE(result.edges_bitset == 0);
        REQUIRE(result.swaps_code == 0);
        REQUIRE(result.number_of_swaps == std::numeric_limits<unsigned>::max());
      }
    }
  }
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
