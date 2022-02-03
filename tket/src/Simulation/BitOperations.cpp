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

#include "BitOperations.hpp"

#include <stdexcept>

#include "Utils/Assert.hpp"

namespace tket {
namespace tket_sim {
namespace internal {

// We want to push back to ExpansionData, but check if it can be combined
// into an existing element.
static void push_back(
    ExpansionData& result, SimUInt single_bit, unsigned left_shift_argument) {
  if (!result.empty() && left_shift_argument == result.back().second) {
    // Since the left shift arguments agree,
    // they can be combined into a single operation,
    // which is more efficient.
    result.back().first |= single_bit;
    return;
  }
  result.emplace_back(single_bit, left_shift_argument);
}

ExpansionData get_expansion_data(
    SimUInt forbidden_bits, unsigned number_of_free_bits) {
  ExpansionData result;
  SimUInt next_bit = 1;

  for (unsigned bits_count = 0; bits_count < number_of_free_bits;
       ++bits_count) {
    auto test_bit = next_bit;
    for (unsigned left_shift_arg = 0;; ++left_shift_arg) {
      if ((test_bit & forbidden_bits) == 0) {
        TKET_ASSERT(test_bit != 0);
        // A free space has been found.
        push_back(result, next_bit, left_shift_arg);
        forbidden_bits |= test_bit;
        break;
      }
      test_bit <<= 1;
    }
    next_bit <<= 1;
  }
  return result;
}

SimUInt get_expanded_bits(const ExpansionData& expansion_data, SimUInt bits) {
  SimUInt result = 0;
  for (const auto& entry : expansion_data) {
    SimUInt block = bits & entry.first;
    block <<= entry.second;
    result |= block;
  }
  return result;
}

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
