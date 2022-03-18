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

#include "WeightSubgrMono/GraphTheoretic/TriangleCounts.hpp"

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/GraphTheoretic/DerivedGraphsUpdater.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

TriangleCounts::TriangleCounts(DerivedGraphsUpdater& updater)
    : m_updater(updater) {}

DerivedGraphStructs::Count TriangleCounts::get_count(VertexWSM v) {
  {
    const auto iter = m_data.find(v);
    if (iter != m_data.end()) {
      return iter->second;
    }
  }
  // It doesn't exist, so must be calculated.
  m_updater.fill_data_in_container(v);

  // The call above has now filled the data!
  return m_data.at(v);
}

void TriangleCounts::fill_count(VertexWSM v, DerivedGraphStructs::Count count) {
  m_data[v] = count;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
