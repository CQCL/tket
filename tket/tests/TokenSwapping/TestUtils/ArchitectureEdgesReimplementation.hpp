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

#include <utility>
#include <vector>

namespace tket {
namespace tsa_internal {
namespace tests {

// We would like to use the SquareGrid Architecture class,
// but the order of edges is not guaranteed (an implementation detail).
// Therefore, we copy the code to have a single, fixed ordering
// for testing purposes with token swapping.
// NOTE: the only important thing is the order of edges,
// NOT the specific vertex labels. The vertices will be relabelled
// in order of appearance by ArchitectureMapping.
std::vector<std::pair<unsigned, unsigned>> get_square_grid_edges(
    unsigned dim_r, const unsigned dim_c, const unsigned layers);

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
