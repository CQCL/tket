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

#pragma once
#include <cstdint>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Functions which are useful for checking for uint overflow, etc.
 * involving only bit shifts and int operations (no doubles).
 */
struct BitFunctions {
  /** How many zeros are there at the end of the binary expansion of x?
   * (So, we can right shift x by this amount and not lose any bits).
   * Returns 64 for x=0.
   * @param x The raw bits.
   * @return The number of consective zero bits, starting from the right (the
   * least significant bit). F(0)=64.
   */
  static unsigned get_number_of_rightmost_zero_bits(std::uint64_t x);

  /** The smallest n such that 2^n > x.
   * Thus, if you replace all bits to the right of the leftmost bit
   * with "1", then this is the total number of bits you will have.
   * @param x The raw bits x.
   * @return The smallest n such that 2^n > x.
   */
  static unsigned get_bit_length(std::uint64_t x);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
