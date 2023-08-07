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

#include "PermutationTestUtils.hpp"

#include <catch2/catch_test_macros.hpp>
#include <vector>
#include <algorithm>

namespace tket {
namespace tsa_internal {
namespace tests {

std::array<unsigned, 6> PermutationTestUtils::get_end_tokens_for_permutation(
    unsigned permutation_hash) {
  REQUIRE(permutation_hash >= 2);
  std::vector<unsigned> digits;
  {
    unsigned perm_hash_copy = permutation_hash;
    while (perm_hash_copy != 0) {
      digits.push_back(perm_hash_copy % 10);
      perm_hash_copy /= 10;
    }
    REQUIRE(!digits.empty());
    REQUIRE(std::is_sorted(digits.cbegin(), digits.cend()));
    REQUIRE(digits[0] >= 2);
    std::reverse(digits.begin(), digits.end());
  }
  unsigned cycle_start_v = 0;
  std::array<unsigned, 6> tokens;
  // No significance to 9999, just a number>5 which stands out
  tokens.fill(9999);
  for (unsigned cycle_length : digits) {
    // We want to enact the cycle (a,b,c,d). Thus a->b, etc. is the vertex
    // mapping. Now "tokens" represents what happens IF the vertex mapping is
    // applied to [0,1,2,...]. Thus, whatever was INITIALLY at vertex "a" (the
    // number "a" itself) should end up at "b", i.e. tokens[b] == a.
    for (unsigned ii = 0; ii < cycle_length; ++ii) {
      const unsigned source_v = cycle_start_v + ii;
      const unsigned target_v = cycle_start_v + ((ii + 1) % cycle_length);
      REQUIRE(source_v != target_v);
      REQUIRE(source_v <= 5);
      REQUIRE(target_v <= 5);
      tokens[target_v] = source_v;
    }
    cycle_start_v += cycle_length;
  }
  REQUIRE(cycle_start_v <= 6);
  for (unsigned ii = cycle_start_v; ii < 6; ++ii) {
    tokens[ii] = ii;
  }
  for (unsigned tok : tokens) {
    REQUIRE(tok < 6);
  }
  return tokens;
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
