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

#include "EdgeSequence.hpp"

#include "Graphs/AdjacencyData.hpp"

using std::size_t;

namespace tket {
namespace graphs {
namespace tests {

EdgeSequence::EdgeSequence(AdjacencyData& _adj_data, RNG& _rng)
    : adjacency_data(_adj_data), rng(_rng) {}

bool EdgeSequence::add_edge(size_t i, size_t j) {
  if (adjacency_data.add_edge(i, j)) {
    edges.emplace_back(i, j);
    return true;
  }
  return false;
}

void EdgeSequence::clear() {
  adjacency_data.clear(0);
  edges.clear();
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
