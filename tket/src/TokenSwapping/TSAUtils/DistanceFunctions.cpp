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

#include "TokenSwapping/DistanceFunctions.hpp"

#include <numeric>
#include <stdexcept>

namespace tket {
namespace tsa_internal {

size_t get_total_home_distances(
    const VertexMapping& vertex_mapping,
    DistancesInterface& distances_calculator) {
  size_t sum_of_distances = 0;
  for (const auto& entry : vertex_mapping) {
    sum_of_distances += distances_calculator(entry.first, entry.second);
  }
  return sum_of_distances;
}

int get_move_decrease(
    const VertexMapping& vertex_mapping, size_t v1, size_t v2,
    DistancesInterface& distances) {
  const auto citer = vertex_mapping.find(v1);
  if (citer == vertex_mapping.cend()) {
    return 0;
  }
  const auto target = citer->second;
  const std::intmax_t v1_to_target = distances(v1, target);
  const std::intmax_t v2_to_target = distances(v2, target);
  return static_cast<int>(v1_to_target - v2_to_target);
}

int get_swap_decrease(
    const VertexMapping& vertex_mapping, size_t v1, size_t v2,
    DistancesInterface& distances) {
  return get_move_decrease(vertex_mapping, v1, v2, distances) +
         get_move_decrease(vertex_mapping, v2, v1, distances);
}

}  // namespace tsa_internal
}  // namespace tket
