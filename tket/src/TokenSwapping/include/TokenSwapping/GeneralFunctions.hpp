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

// This is for "leftover" functions not specifically linked to token swapping
// which are candidates for being used and moved elsewhere,
// e.g. the main src/Utils folder.

#include <map>
#include <optional>
#include <set>

#include "Utils/Assert.hpp"

namespace tket {
namespace tsa_internal {

/** Returns the value in a map corresponding to a key, IF it exists,
 *  or an empty optional object if it does not.
 *  @param map The std::map object.
 *  @param key The key.
 *  @return The value if it exists, or an empty optional value if it doesn't.
 */
template <class K, class V>
std::optional<V> get_optional_value(const std::map<K, V>& map, const K& key) {
  const auto citer = map.find(key);
  if (citer == map.cend()) {
    return {};
  }
  return citer->second;
}

/** The key->value mapping is required to be bijective (reversible).
 *  @param map The std::map object.
 *  @return Another std::map, with the key->value mappings reversed.
 */
template <class K, class V>
std::map<V, K> get_reversed_map(const std::map<K, V>& map) {
  std::map<V, K> reversed_map;
  for (const auto& entry : map) {
    reversed_map[entry.second] = entry.first;
  }
  TKET_ASSERT(map.size() == reversed_map.size());
  return reversed_map;
}

/** Finds the rightmost "one" (least significant bit)
 * occurring in the binary expansion of x, an unsigned integer type.
 * Returns the bit, whilst also removing it from x.
 * @param x The original unsigned integer type, which will have one bit removed
 * (or remain at zero if already at zero).
 * @return The bit which was removed from x (or 0 if none was removed).
 */
template <class UINT>
static UINT get_rightmost_bit(UINT& x) {
  // Standard bit hack: decrementing 10000... gives 01111...
  // E.g., consider:
  //      x = 001101011010000
  //    x-1 = 001101011001111
  // ~(x-1) = 110010100110000
  // Therefore, AND x with ~(x-1).

  // No "if" statements; unsigned int wraparound is allowed.
  UINT y = x;
  --y;
  y = ~y;
  const UINT bit = (x & y);
  x ^= bit;
  return bit;
}

}  // namespace tsa_internal
}  // namespace tket
