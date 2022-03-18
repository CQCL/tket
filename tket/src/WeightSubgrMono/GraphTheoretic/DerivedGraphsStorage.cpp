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

#include "WeightSubgrMono/GraphTheoretic/DerivedGraphStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DerivedGraphStructs::Iter
DerivedGraphsStorage::get_new_neighbours_and_counts_iter() {
  m_stored_neighbours_and_counts.emplace_front();
  return m_stored_neighbours_and_counts.begin();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
