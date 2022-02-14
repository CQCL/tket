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

#include "CyclicShiftCostEstimate.hpp"

#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

CyclicShiftCostEstimate::CyclicShiftCostEstimate(
    const std::vector<size_t>& vertices, DistancesInterface& distances) {
  TKET_ASSERT(vertices.size() >= 2);
  // We first work out the total distance v(0)->v(1)-> .. -> v(n) -> v(0).
  // If we snip out v(i)->v(i+1), the remaining path tells us how many swaps
  // we need. So, we must snip out the LARGEST distance(v(i), v(i+1)).
  size_t largest_distance = distances(vertices.back(), vertices[0]);
  size_t total_distance = largest_distance;

  if (vertices.size() == 2) {
    start_v_index = 0;
  } else {
    // The value i such that distance(v(i), v(i+1)) is largest.
    size_t v_index_with_largest_distance = vertices.size() - 1;
    for (size_t ii = 0; ii + 1 < vertices.size(); ++ii) {
      const auto distance_i = distances(vertices[ii], vertices[ii + 1]);
      TKET_ASSERT(distance_i > 0);
      total_distance += distance_i;
      if (distance_i < largest_distance) {
        largest_distance = distance_i;
        v_index_with_largest_distance = ii;
      }
    }
    // Now, remove the largest distance again...
    total_distance -= largest_distance;
    // We've snipped out (v[i], v[i+1]), so logically we start from v[i+1].
    start_v_index = (v_index_with_largest_distance + 1) % vertices.size();
  }
  // To enact an abstract cyclic shift [a,b,c,d],
  // choose abstract swaps (cd), (bc), (ab).
  // The number of CONCRETE swaps to enact an abstract swap (xy) is
  // 2.dist(x,y) - 1.
  // e.g., to swap x,y along the path [x,u,v,y], dist(x,y)=3,
  // we use 5 concrete vertex swaps (xu), (uv), (vy), (uv), (xu).
  // What we've currently stored is the sum of dist(x,y),
  // and clearly (sum)(-1) = -(Number of terms in the sum).
  estimated_concrete_swaps = 2 * total_distance;
  TKET_ASSERT(estimated_concrete_swaps > vertices.size() - 1);
  estimated_concrete_swaps -= vertices.size() - 1;
}

}  // namespace tsa_internal
}  // namespace tket
