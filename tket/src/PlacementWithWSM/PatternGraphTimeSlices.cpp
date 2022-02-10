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

#include "PlacementWithWSM/PatternGraphTimeSlices.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

static std::size_t get_time_for_this_interaction(
    const PatternGraphTimeSlices& time_slices,
    const std::set<VertexWSM>& vertices) {
  // This will be one more than a valid index.
  std::size_t current_index = time_slices.time_sliced_data.size();
  for (; current_index > 0; --current_index) {
    const auto& possibly_blocking_slice =
        time_slices.time_sliced_data[current_index - 1];
    // Does it actually block?
    for (const auto& pair : possibly_blocking_slice) {
      if (vertices.count(pair.first) != 0 || vertices.count(pair.second) != 0) {
        return current_index;
      }
    }
    // If we reached here, nothing blocks.
  }
  return current_index;
}

PatternGraphTimeSlices::PatternGraphTimeSlices(
    const std::vector<std::set<VertexWSM>>& gates_in_order) {
  for (const auto& entry : gates_in_order) {
    if (entry.size() < 2) {
      continue;
    }
    const auto time = get_time_for_this_interaction(*this, entry);
    if (time >= time_sliced_data.size()) {
      time_sliced_data.resize(time + 1);
    }
    auto& slice = time_sliced_data[time];
    auto citer = entry.cbegin();
    auto previous_vertex = *citer;
    for (;;) {
      ++citer;
      if (citer == entry.cend()) {
        break;
      }
      auto current_vertex = *citer;
      slice.emplace_back(previous_vertex, current_vertex);
      previous_vertex = current_vertex;
    }
  }
}

GraphEdgeWeights PatternGraphTimeSlices::get_weights(
    const std::vector<WeightWSM>& single_edge_weights_at_all_times) const {
  TKET_ASSERT(
      single_edge_weights_at_all_times.size() >= time_sliced_data.size());

  GraphEdgeWeights result;
  for (unsigned ii = 0; ii < time_sliced_data.size(); ++ii) {
    const auto& single_edge_weight = single_edge_weights_at_all_times[ii];
    for (const auto& pair : time_sliced_data[ii]) {
      auto& weight = result[pair];
      weight = get_sum_or_throw(weight, single_edge_weight);
    }
  }
  return result;
}

GraphEdgeWeights PatternGraphTimeSlices::get_weights(
    const WeightParameters& parameters) const {
  return get_weights(
      parameters.get_single_edge_weights(time_sliced_data.size()));
}

std::vector<WeightWSM>
PatternGraphTimeSlices::WeightParameters::get_single_edge_weights(
    unsigned size) const {
  TKET_ASSERT(time_zero_edge_weight > 0);
  TKET_ASSERT(final_time_edge_weight > 0);
  TKET_ASSERT(size > 1);
  std::vector<WeightWSM> result(size);
  result[0] = time_zero_edge_weight;
  result.back() = final_time_edge_weight;
  if (size == 2) {
    return result;
  }
  // Just do the simplest: linear interpolation.
  // Check for overflow.
  {
    const auto big_weight =
        get_sum_or_throw(time_zero_edge_weight, final_time_edge_weight);
    get_product_or_throw(big_weight, WeightWSM(2 * size));
  }

  // We increase from the smallest to the highest weight.
  // Let "initial_index" have the smaller weight,
  // "final_index" the larger weight.
  int index_step = 1;
  int initial_index = 0;
  unsigned final_index = size - 1;

  if (result[initial_index] > result[final_index]) {
    // Actually, this is the usual order:
    // early times have higher weights!
    index_step = -1;
    initial_index = final_index;
    final_index = 0;
  }
  const WeightWSM lower_weight = result[initial_index];
  const WeightWSM weight_diff = result[final_index] - lower_weight;

  for (unsigned counter = 0; counter + 2 < size; ++counter) {
    initial_index += index_step;
    result[initial_index] =
        lower_weight + (weight_diff * (counter + 1)) / (size - 1);
  }
  return result;
}

}  // namespace tket
