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
#include "DerivedGraph.hpp"
#include "TriangleCounts.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DerivedGraphsUpdater;

/** The actual graphs are stored here, initialised
 * so that they know how to evaluate things lazily and propagate calculations
 * automatically as necessary; the caller need not worry.
 */
struct DerivedGraphs {
  DerivedGraph d2_graph;
  DerivedGraph d3_graph;

  /** This is a derived vertex property rather than a graph
   * (although one could also regard it as a totally disconnected
   * graph with vertex labels).
   */
  TriangleCounts triangle_counts;

  DerivedGraphs(DerivedGraphsUpdater& updater);
};


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
