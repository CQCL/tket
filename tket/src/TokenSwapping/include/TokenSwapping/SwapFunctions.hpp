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

#include <map>
#include <optional>
#include <utility>

#include "TokenSwapping/VectorListHybrid.hpp"

namespace tket {

typedef std::pair<size_t, size_t> Swap;
typedef VectorListHybrid<Swap> SwapList;
typedef SwapList::ID SwapID;

/** No distinction between (v1, v2) and (v2, v1).
 *  Will ensure that v1<v2. So, swaps will automatically be comparable using
 *  the standard ordering, so can be used as keys in sets and maps.
 *  @param v1 The first vertex.
 *  @param v2 The second vertex (throws if equal to the first).
 *  @return A swap object.
 */
Swap get_swap(size_t v1, size_t v2);

/** Do the swaps act on 4 different vertices?
 *  @param swap1 The first swap.
 *  @param swap2 The second swap.
 *  @return whether the swaps have no vertices in common
 *      (and thus, commute with each other).
 */
bool disjoint(const Swap& swap1, const Swap& swap2);

}  // namespace tket
