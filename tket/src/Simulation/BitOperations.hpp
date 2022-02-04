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

#include <cstdint>
#include <utility>
#include <vector>

namespace tket {
namespace tket_sim {
namespace internal {

// Used to represent qubit indices.
typedef std::uint_fast32_t SimUInt;

// Application: we have a mask like 000101100111001.
// Those positions with "1" are "forbidden";
// the other positions are "free": we are allowed to set them.
// In the example 000101100111001, there are 7 forbidden positions
// (ones) and 8 free positions (zeros).
// If we have a length 8 binary string "abcdefgh",
// we want to insert zeros and stretch it to become abc0d00ef000gh0
// (so that the forbidden positions are all zero and thus may be set
// by some other function by a simple OR).
//
// This is done by ANDing abcdefgh with 00000011, 00001100, 00010000,
// 11100000 in turn, shifting each result left by a certain number
// of bits, then ORing them together.
// In each element, "first" is the mask that we must AND with;
// "second" is the argument to the shift left operator.
typedef std::vector<std::pair<SimUInt, unsigned>> ExpansionData;

ExpansionData get_expansion_data(
    SimUInt forbidden_bits, unsigned number_of_free_bits);

// Given the bit string "abcdefgh", stretch it to become abc0d00ef000gh0
// according to the expansion data.
SimUInt get_expanded_bits(const ExpansionData& expansion_data, SimUInt bits);

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
