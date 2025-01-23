// Copyright Quantinuum
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

#include "tkwsm/GraphTheoretic/VertexRelabelling.hpp"

#include <sstream>
#include <stdexcept>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

VertexRelabelling::VertexRelabelling(GraphEdgeWeights edges_and_weights)
    : new_edges_and_weights(std::move(edges_and_weights)) {
  for (const std::pair<const EdgeWSM, WeightWSM>& entry :
       new_edges_and_weights) {
    const VertexWSM& v1 = entry.first.first;
    const VertexWSM& v2 = entry.first.second;
    if (v1 == v2) {
      throw std::runtime_error("Loop found in graph.");
    }
    old_to_new_vertex_labels[v1];
    old_to_new_vertex_labels[v2];
    const std::optional<WeightWSM> weight_opt =
        get_optional_value(new_edges_and_weights, std::make_pair(v2, v1));
    if (weight_opt) {
      const WeightWSM other_weight = weight_opt.value();
      if (entry.second != other_weight) {
        throw std::runtime_error("reversed edge has different weight");
      }
    }
  }
  if (old_to_new_vertex_labels.empty()) {
    throw std::runtime_error("Input graph has no edges");
  }
  number_of_vertices = old_to_new_vertex_labels.size();
  TKET_ASSERT(number_of_vertices >= 2);
  if (old_to_new_vertex_labels.cbegin()->first == 0 &&
      old_to_new_vertex_labels.crbegin()->first == number_of_vertices - 1) {
    // Already contiguous!
    old_to_new_vertex_labels.clear();
    return;
  }
  new_to_old_vertex_labels.reserve(number_of_vertices);
  for (std::pair<const VertexWSM, unsigned>& entry : old_to_new_vertex_labels) {
    const unsigned new_vertex = new_to_old_vertex_labels.size();
    new_to_old_vertex_labels.push_back(entry.first);
    entry.second = new_vertex;
  }
  TKET_ASSERT(new_to_old_vertex_labels.size() == number_of_vertices);
  TKET_ASSERT(old_to_new_vertex_labels.at(new_to_old_vertex_labels[0]) == 0);

  // Another copy is unavoidable (the vertices aren't contiguous,
  // we MUST relabel; but this means changing the KEYS of the map
  // new_edges_and_weights)
  const auto original_edges_and_weights = new_edges_and_weights;
  new_edges_and_weights.clear();
  for (const auto& entry : original_edges_and_weights) {
    const VertexWSM& old_v1 = entry.first.first;
    const VertexWSM& old_v2 = entry.first.second;

    new_edges_and_weights[get_edge(
        old_to_new_vertex_labels.at(old_v1),
        old_to_new_vertex_labels.at(old_v2))] = entry.second;
  }
}

unsigned VertexRelabelling::get_new_label(VertexWSM v) const {
  if (old_to_new_vertex_labels.empty()) {
    return v;
  }
  return old_to_new_vertex_labels.at(v);
}

VertexWSM VertexRelabelling::get_old_label(unsigned v) const {
  if (new_to_old_vertex_labels.empty()) {
    return v;
  }
  return new_to_old_vertex_labels.at(v);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
