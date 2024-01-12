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
#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
class RNG;
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

/** We want to extract only a few random bits each time,
 * assuming that generating 64 random bits is relatively slow.
 * TODO move into Utils, or maybe RNG directly.
 */
class FastRandomBits {
 public:
  FastRandomBits();

  /** Will call the RNG object as little as possible. */
  std::uint64_t get_random_bits(RNG& rng, unsigned number_of_bits);

 private:
  std::uint64_t m_bits;
  unsigned m_number_of_random_bits;
};

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
