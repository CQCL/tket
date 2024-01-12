// Copyright 2019-2024 Cambridge Quantum Computing
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
#include "GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** If we have vertices {0,1,2,...,N} instead of arbitrary integers,
 * then we can use std::vector instead of std::map in many places.
 * This is actually a significant speedup (>10%).
 * Thus, when we pass in a problem, we relabel internally to ensure
 * that this holds; and convert back to the original vertices when the user
 * requests a solution.
 * This just contains data to do the conversions.
 */
struct VertexRelabelling {
  // KEY: is the old vertex label.
  // VALUE: the new vertex label (contiguous: {0,1,2,...,N}).
  // If empty, it means that the old vertices were ALREADY
  // contiguous, so no need to relabel them.
  std::map<VertexWSM, unsigned> old_to_new_vertex_labels;

  // Element[i] is the old vertex label for new vertex i.
  // If empty, the old vertices were ALREADY contiguous.
  std::vector<VertexWSM> new_to_old_vertex_labels;

  unsigned number_of_vertices;

  // Uses the new vertex labels.
  GraphEdgeWeights new_edges_and_weights;

  explicit VertexRelabelling(GraphEdgeWeights edges_and_weights);

  unsigned get_new_label(VertexWSM v) const;

  VertexWSM get_old_label(unsigned v) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
