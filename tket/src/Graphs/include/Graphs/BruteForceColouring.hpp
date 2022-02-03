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

#include <cstddef>
#include <map>
#include <memory>

namespace tket {
namespace graphs {

class ColouringPriority;

/**
 * For colouring vertices in a single connected component of a graph.
 * Although it's brute force, because it goes down the vertices
 * in a fixed order, as specified by "ColouringPriority",
 * (like traversing a "game tree"), there is the potential for pruning,
 * and hence it can be much quicker than simply trying every possible colouring.
 * It only looks for ONE colouring.
 */
class BruteForceColouring {
 public:
  /**
   * All the calculation is done upon construction.
   *
   * @param priority This ColouringPriority object already contains all the data
   *  about vertices in this component, and a sensible traversal order.
   * @param suggested_number_of_colours A hint that you think this many colours
   * are needed. If you set it too high you may end up with a suboptimal
   * colouring, but it might be quicker.
   */
  BruteForceColouring(
      const ColouringPriority& priority,
      std::size_t suggested_number_of_colours = 0);

  /**
   * The colours found for this component (already calculated during
   * construction). It is a  vertex -> colour  mapping.
   */
  const std::map<std::size_t, std::size_t>& get_colours() const;

  ~BruteForceColouring();

 private:
  struct Impl;
  std::unique_ptr<Impl> m_pimpl;
};

}  // namespace graphs
}  // namespace tket
