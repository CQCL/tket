// Copyright 2019-2021 Cambridge Quantum Computing
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

#ifndef _TKET_AbstractGraph_H
#define _TKET_AbstractGraph_H

#include <utility>
#include <vector>

namespace tket::graphs {

/**
 * Abstract class for representing graphs.
 *
 * @tparam T type of nodes in the graph
 */
template <typename T>
class AbstractGraph {
 public:
  /** Construct an empty graph */
  AbstractGraph() : nodes_() {}

  /** Construct from list of nodes */
  explicit AbstractGraph(const std::vector<T>& nodes) : nodes_(nodes) {}

 protected:
  using Edge = std::pair<T, T>;

 private:
  std::vector<T> nodes_;
};

}  // namespace tket::graphs

#endif  // _TKET_AbstractGraph_H
