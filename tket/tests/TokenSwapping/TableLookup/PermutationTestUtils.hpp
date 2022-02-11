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

#pragma once

#include <array>

namespace tket {
namespace tsa_internal {
namespace tests {

// See CanonicalRelabelling.hpp for an explanation of the "permutation hash".

struct PermutationTestUtils {
  /** Given a permutation hash, return the final tokens after performing that
   *  mapping on the vertices 0,1,2,...,5 in the canonical way.
   *  @param permutation_hash A decimal number representing a permutation on
   * {0,1,...,5}.
   *  @return The numbers {0,1,2,...,5} giving the final tokens, if we perform
   * the permutation, with each start token label equalling the vertex label.
   */
  static std::array<unsigned, 6> get_end_tokens_for_permutation(
      unsigned permutation_hash);
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
