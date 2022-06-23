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

#include "DebugFunctions.hpp"

#include <sstream>

namespace tket {
namespace tsa_internal {

std::string str(const VertexMapping& vertex_mapping) {
  std::stringstream ss;
  ss << "VM:";
  for (const auto& entry : vertex_mapping) {
    ss << " " << entry.first << "->" << entry.second << " ";
  }
  return ss.str();
}

std::string str(const SwapList& swaps) { return str(swaps.to_vector()); }

std::string str(const std::vector<Swap>& swaps) {
  std::stringstream ss;
  for (auto swap : swaps) {
    ss << " (" << swap.first << "," << swap.second << ") ";
  }
  return ss.str();
}

size_t get_swaps_lower_bound(
    const VertexMapping& vertex_mapping,
    DistancesInterface& distances_calculator) {
  // Each swap decreases the sum by at most 2 (and more likely 1 in many cases,
  // if the mapping is sparse), so we need  >= sum/2. But it's an integer of
  // course.
  return (get_total_home_distances(vertex_mapping, distances_calculator) + 1) /
         2;
}

}  // namespace tsa_internal
}  // namespace tket
