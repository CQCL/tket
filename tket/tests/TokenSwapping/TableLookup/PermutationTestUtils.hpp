#ifndef _TKET_TESTS_TokenSwapping_TableLookup_PermutationTestUtils_H_
#define _TKET_TESTS_TokenSwapping_TableLookup_PermutationTestUtils_H_

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
#endif
