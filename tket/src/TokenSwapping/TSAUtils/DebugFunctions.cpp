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

#include "TokenSwapping/DebugFunctions.hpp"

#include <sstream>

namespace tket {
namespace tsa_internal {

// GCOVR_EXCL_START
std::string str(const VertexMapping& vertex_mapping) {
  std::stringstream ss;
  ss << "VM:";
  for (const auto& entry : vertex_mapping) {
    ss << " " << entry.first << "->" << entry.second << " ";
  }
  return ss.str();
}
// GCOVR_EXCL_STOP

std::string str(const SwapList& swaps) { return str(swaps.to_vector()); }

std::string str(const std::vector<Swap>& swaps) {
  std::stringstream ss;
  for (auto swap : swaps) {
    ss << " (" << swap.first << "," << swap.second << ") ";
  }
  return ss.str();
}

}  // namespace tsa_internal
}  // namespace tket
