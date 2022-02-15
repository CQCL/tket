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

#include "GetRandomSet.hpp"

#include <numeric>

#include "Utils/Assert.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

std::set<size_t> get_random_set(
    RNG& rng, size_t sample_size, size_t population_size) {
  TKET_ASSERT(sample_size <= population_size);

  std::set<size_t> result;
  if (sample_size == 0 || population_size == 0) {
    return result;
  }
  if (sample_size < population_size / 2) {
    while (result.size() < sample_size) {
      result.insert(rng.get_size_t(population_size - 1));
    }
    return result;
  }
  std::vector<size_t> elems(population_size);
  std::iota(elems.begin(), elems.end(), 0);
  rng.do_shuffle(elems);
  for (const auto& elem : elems) {
    result.insert(elem);
    if (result.size() == sample_size) {
      return result;
    }
  }
  TKET_ASSERT(!"get_random_set: dropped out of loop");
  return result;
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
