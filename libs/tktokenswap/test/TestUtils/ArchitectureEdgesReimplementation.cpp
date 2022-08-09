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

#include "ArchitectureEdgesReimplementation.hpp"

#include <catch2/catch_test_macros.hpp>

namespace tket {
namespace tsa_internal {
namespace tests {

// This is just copied from Architecture.cpp,
// but we WANT it to remain fixed for testing purposes;
// do NOT keep in sync!
std::vector<std::pair<unsigned, unsigned>> get_square_grid_edges(
    unsigned dim_r, const unsigned dim_c, const unsigned layers) {
  // A trivial injective hash function on the cuboid.
  const auto vertex = [dim_r, dim_c, layers](
                          unsigned ver, unsigned hor, unsigned l) -> unsigned {
    REQUIRE(ver < dim_r);
    REQUIRE(hor < dim_c);
    REQUIRE(l < layers);
    return ver + dim_r * (hor + dim_c * l);
  };

  std::vector<std::pair<unsigned, unsigned>> edges;
  for (unsigned l = 0; l < layers; l++) {
    for (unsigned ver = 0; ver < dim_r; ver++) {
      for (unsigned hor = 0; hor < dim_c; hor++) {
        const auto n = vertex(ver, hor, l);
        if (hor != dim_c - 1) {
          const auto h_neighbour = vertex(ver, hor + 1, l);
          edges.push_back({n, h_neighbour});
        }
        if (ver != dim_r - 1) {
          const auto v_neighbour = vertex(ver + 1, hor, l);
          edges.push_back({n, v_neighbour});
        }
        if (l != layers - 1) {
          const auto l_neighbour = vertex(ver, hor, l + 1);
          edges.push_back({n, l_neighbour});
        }
      }
    }
  }
  return edges;
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
