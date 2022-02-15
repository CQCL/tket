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

#include <set>

#include "Utils/RNG.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

/** Return a random subset of given size from the population {0,1,2,...,N}.
 *  @param rng A random number generator.
 *  @param sample_size The desired size of the returned set.
 *  @param population_size The number of elements in the population (an interval
 *    of nonnegative integers, starting at 0).
 *  @return A set of numbers. Throws upon invalid parameters.
 */
std::set<size_t> get_random_set(
    RNG& rng, size_t sample_size, size_t population_size);

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
