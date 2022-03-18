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

#include "TokenSwapping/SwapFunctions.hpp"

#include <sstream>
#include <stdexcept>

namespace tket {

Swap get_swap(size_t v1, size_t v2) {
  if (v1 == v2) {
    std::stringstream ss;
    ss << "get_swap : for equal vertices v1 = v2 = v_" << v1;
    throw std::runtime_error(ss.str());
  }
  if (v1 < v2) {
    return std::make_pair(v1, v2);
  }
  return std::make_pair(v2, v1);
}

bool disjoint(const Swap& s1, const Swap& s2) {
  return s1.first != s2.first && s1.first != s2.second &&
         s1.second != s2.first && s1.second != s2.second;
}

}  // namespace tket
