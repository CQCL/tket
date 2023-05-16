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

#include "tkwsm/InitPlacement/FastRandomBits.hpp"

#include <tkassert/Assert.hpp>
#include <tkrng/RNG.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

// Element[i] is the mask to extract i+1 bits.
static std::vector<std::uint64_t> get_masks() {
  std::vector<std::uint64_t> masks(64);
  masks[0] = 1;
  for (unsigned ii = 1; ii <= 63; ++ii) {
    masks[ii] = (masks[ii - 1] << 1) | masks[ii - 1];
  }
  return masks;
}

static const std::vector<std::uint64_t>& get_masks_global_ref() {
  static const auto masks_global = get_masks();
  return masks_global;
}

FastRandomBits::FastRandomBits() : m_bits(0), m_number_of_random_bits(0) {}

std::uint64_t FastRandomBits::get_random_bits(
    RNG& rng, unsigned number_of_bits) {
  TKET_ASSERT(number_of_bits >= 1);
  TKET_ASSERT(number_of_bits <= 64);
  if (number_of_bits <= m_number_of_random_bits) {
    const std::uint64_t bits_to_return =
        m_bits & get_masks_global_ref()[number_of_bits - 1];
    m_bits >>= number_of_bits;
    m_number_of_random_bits -= number_of_bits;
    return bits_to_return;
  }
  // There are not enough random bits, so we'll need to generate more.
  // But first, use the existing bits.
  std::uint64_t bits_to_return = m_bits;
  const unsigned number_of_extra_bits =
      number_of_bits - m_number_of_random_bits;
  bits_to_return <<= number_of_extra_bits;
  m_bits = rng();
  bits_to_return |= m_bits & get_masks_global_ref()[number_of_extra_bits - 1];
  m_bits >>= number_of_extra_bits;
  m_number_of_random_bits = 64 - number_of_bits;
  return bits_to_return;
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
