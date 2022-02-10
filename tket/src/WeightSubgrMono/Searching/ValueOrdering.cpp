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

#include "WeightSubgrMono/Searching/ValueOrdering.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SharedData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

VertexWSM ValueOrdering::get_target_value(
    const std::set<VertexWSM>& possible_values, const SharedData& shared_data) {
  TKET_ASSERT(possible_values.size() >= 2);
  unsigned max_degree = 0;
  VertexWSM best_tv = 0;
  for (auto tv : possible_values) {
    const unsigned degree =
        shared_data.fixed_data.target_neighbours_data.get_degree(tv);
    if (degree > max_degree) {
      max_degree = degree;
      best_tv = tv;
    }
  }
  return best_tv;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
