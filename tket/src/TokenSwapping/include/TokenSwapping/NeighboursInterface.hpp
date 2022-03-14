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

#include <cstdint>
#include <vector>

namespace tket {

/** What are the adjacent vertices to a given vertex on a graph?
 *  For larger, sparse graphs, it might
 *  calculate and store neighbours only when required.
 */
class NeighboursInterface {
 public:
  /** Returns the neighbours of the given vertex.
   *  The vector of neighbours is required to be stored internally.
   *  However, no guarantee that the reference will remain valid
   *  once another function call occurs.
   *  By default, throws (not implemented).
   *  (It's useful to be able to create a "null" object like this,
   *  because some algorithms don't actually need a neighbours object,
   *  but others do).
   *  @param vertex A vertex.
   *  @return A sorted list of all adjacent vertices, stored internally.
   */
  virtual const std::vector<size_t>& operator()(size_t vertex);

  virtual ~NeighboursInterface();
};

}  // namespace tket
