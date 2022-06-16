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
#include <catch2/catch_test_macros.hpp>
#include <numeric>

#include "PermutationTestUtils.hpp"
#include "TokenSwapping/SwapConversion.hpp"
#include "TokenSwapping/SwapListOptimiser.hpp"
#include "TokenSwapping/SwapSequenceTable.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

// Extra redundant data in the table slows it down,
// but does not affect the returned results.
// But the stored swap sequences are used directly without further checks
// or optimisations, so they should be as close to optimal as possible.
static void test_irreducibility_of_codes(
    unsigned permutation_hash, const vector<SwapSequenceTable::Code>& codes,
    SwapListOptimiser& optimiser, SwapList& swap_list) {
  for (auto& code : codes) {
    swap_list.fast_clear();
    auto swap_sequence_hash_copy = code;
    while (swap_sequence_hash_copy != 0) {
      const Swap& swap =
          SwapConversion::get_swap_from_hash(swap_sequence_hash_copy & 0xF);
      swap_list.push_back(swap);
      swap_sequence_hash_copy >>= 4;
    }
    const auto initial_number_of_swaps = swap_list.size();

    // We don't yet have good theoretical results about order of passes,
    // so just try all of them.
    optimiser.optimise_pass_with_zero_travel(swap_list);
    REQUIRE(initial_number_of_swaps == swap_list.size());
    optimiser.optimise_pass_with_token_tracking(swap_list);
    REQUIRE(initial_number_of_swaps == swap_list.size());

    // This may reorder the swaps, without reducing.
    optimiser.optimise_pass_with_frontward_travel(swap_list);
    REQUIRE(initial_number_of_swaps == swap_list.size());

    // We'd LIKE to have a theorem assuring us that this pass isn't necessary
    // after the previous passes, but currently we don't.
    optimiser.optimise_pass_with_token_tracking(swap_list);
    REQUIRE(initial_number_of_swaps == swap_list.size());
    optimiser.optimise_pass_with_zero_travel(swap_list);
    REQUIRE(initial_number_of_swaps == swap_list.size());
  }
}

// All the swap sequences encoded in the vector should enact
// the given permutation.
static void test_correctness_of_codes(
    unsigned permutation_hash, const vector<SwapSequenceTable::Code>& codes) {
  REQUIRE(codes.size() >= 2);

  // Reconstruct the desired permutation from the hash.
  const auto expected_tokens =
      PermutationTestUtils::get_end_tokens_for_permutation(permutation_hash);

  // Element i is the token at vertex i.
  // We start with tokens 0,1,2,...,5 on vertices 0,1,2,...,5,
  // then perform the swaps.
  std::array<unsigned, 6> tokens;
  for (const auto& code : codes) {
    std::iota(tokens.begin(), tokens.end(), 0);
    auto swap_sequence_hash_copy = code;
    unsigned number_of_swaps = 0;
    while (swap_sequence_hash_copy != 0) {
      const Swap& swap =
          SwapConversion::get_swap_from_hash(swap_sequence_hash_copy & 0xF);
      swap_sequence_hash_copy >>= 4;
      std::swap(tokens[swap.first], tokens[swap.second]);
      ++number_of_swaps;
    }
    REQUIRE(number_of_swaps >= 1);

    // Actually, 16 is the maximum.
    CHECK(number_of_swaps <= 12);
    REQUIRE(tokens == expected_tokens);
  }
}

// The swap sequences encoded in the vector should not have
// any redundancies: if sequences S1, S2 have edge bitsets E1, E2
// (i.e., E(j) is the set of swaps used in S(j)), AND give the same permutation,
// then E1 != E2. (No point in having both).
// Also, if E1 is a subset of E2, then length(S2) < length(S1).
// (Otherwise, S2 would be a pointless entry: whenever S2 is possible,
// S1 is also possible, with an equal or smaller number of swaps).
static void test_redundancies(
    unsigned permutation_hash, const vector<SwapSequenceTable::Code>& codes) {
  vector<unsigned> edge_bitsets;
  edge_bitsets.reserve(codes.size());
  for (const auto& code : codes) {
    edge_bitsets.push_back(SwapConversion::get_edges_bitset(code));
  }
  // Crude quadratic algorithm to check which codes are redundant.
  // Don't rely on sorted codes.
  for (unsigned ii = 0; ii < codes.size(); ++ii) {
    for (unsigned jj = 0; jj < codes.size(); ++jj) {
      if (ii == jj) {
        continue;
      }
      const auto intersection = edge_bitsets[ii] & edge_bitsets[jj];
      const bool e1_subset_of_e2 = (intersection == edge_bitsets[ii]);
      const auto num_swaps1 = SwapConversion::get_number_of_swaps(codes[ii]);
      const auto num_swaps2 = SwapConversion::get_number_of_swaps(codes[jj]);

      if (e1_subset_of_e2 && num_swaps1 <= num_swaps2) {
        INFO(
            "For perm.hash "
            << permutation_hash << ", Code 1: 0x" << std::hex << codes[ii]
            << " only uses swaps from code 2: 0x" << codes[jj]
            << ", and uses the same or fewer swaps (" << std::dec << num_swaps1
            << " vs " << num_swaps2
            << "). Thus code 2 is pointless and could be removed.");
        CHECK(false);
      }
    }
  }
}

// Checks that all entries returned by the table do actually
// give the required permutation of vertices.
SCENARIO("Fixed table entries test") {
  const auto table = SwapSequenceTable::get_table();
  // const auto table = get_new_table();
  SwapListOptimiser optimiser;
  SwapList swap_list;
  unsigned total_entries = 0;
  for (const auto& entry : table) {
    REQUIRE(entry.first >= 2);
    test_correctness_of_codes(entry.first, entry.second);
    test_irreducibility_of_codes(
        entry.first, entry.second, optimiser, swap_list);
    test_redundancies(entry.first, entry.second);

    // No duplication. Not necessary, but a good test.
    CHECK(std::is_sorted(entry.second.cbegin(), entry.second.cend()));
    CHECK(
        std::adjacent_find(entry.second.cbegin(), entry.second.cend()) ==
        entry.second.cend());

    // NOTE: we should really also test that inverse mappings are not stored in
    // the table. This was previously true, but a negligibly small number of
    // entries have crept in. They're a bit fiddly to track down and remove, so
    // forget about them for now. (Confusion: within each permutation hash, e.g.
    // 32 corresponding to (012)(34)(5), the INVERSE mapping is (021)(34)(5).
    // This will have the same permutation hash, but of course vertices must be
    // RELABELLED. To find the inverse entry in the table, we cannot JUST
    // reverse the swaps, we also need to relabel them.
    /// TODO: test for, track down and remove redundant inverse entries.
    total_entries += entry.second.size();
  }
  CHECK(total_entries == 7939);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
