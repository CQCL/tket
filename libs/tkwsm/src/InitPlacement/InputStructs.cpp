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

#include "tkwsm/InitPlacement/InputStructs.hpp"

#include <sstream>
#include <stdexcept>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

PatternGraphData::PatternGraphData(DebugNoInputData) {}

// Interpolate the weights, to get a time->actual weight mapping.
// Element[i] gives the weight to use for time i.
static std::vector<WeightWSM> get_interpolated_weights(
    WeightWSM initial_gate_weight, WeightWSM final_gate_weight, unsigned size) {
  if (size < 1 || initial_gate_weight < final_gate_weight ||
      final_gate_weight == 0) {
    std::stringstream ss;
    ss << "WSM initial qubit placement: invalid input pattern data: size="
       << size << ", initial_gate_weight=" << initial_gate_weight
       << ", final_gate_weight=" << final_gate_weight;
    throw std::runtime_error(ss.str());
  }
  std::vector<WeightWSM> weights(size);
  weights[0] = initial_gate_weight;
  if (size == 1) {
    return weights;
  }
  weights.back() = final_gate_weight;
  if (size == 2) {
    return weights;
  }
  const WeightWSM final_index = size - 1;
  for (WeightWSM ii = 1; ii < final_index; ++ii) {
    // Just use arithmetic progession, rather than geometric progression.
    // A fun exercise to implement nearest-integer geometric progressions
    // using only int operations (no std::exp or std::log),
    // but not worth the trouble.

    // CHECK: when ii==0, should be K.init.weight; when ii==size-1,
    // should be K.(final weight).
    // Here, K = final_index.
    //  weights[ii] = final_gate_weight * ii
    //      + initial_gate_weight * (final_index-ii);
    // But, beware of overflow!
    const WeightWSM term1 = get_product_or_throw(final_gate_weight, ii);
    const WeightWSM term2 =
        get_product_or_throw<WeightWSM>(initial_gate_weight, final_index - ii);
    weights[ii] = get_sum_or_throw(term1, term2);

    // We want the nearest integer to w / final_index.
    // We cannot use std::round(double x) because we only want int operations.
    // If we set  k = w / final_index using C++ int division, and x
    // the exact value using real division, then
    // k <= x < k+1.
    // Which is smaller, x-k or k+1 - x?
    // Multiplying up, w - k*(final_index) or (k+1)*(final_index) - w ?
    // Should be a fancier way to do this.
    // Mathematically,
    //      w - k*(final_index) <= (k+1)*(final_index) - w
    //  <==>
    //    2w <= (2k+1)*final_index.
    const WeightWSM kk = weights[ii] / final_index;
    const WeightWSM lhs = get_product_or_throw<WeightWSM>(2, weights[ii]);
    const WeightWSM two_k = get_product_or_throw<WeightWSM>(2, kk);
    const WeightWSM two_k_plus_one = get_sum_or_throw<WeightWSM>(two_k, 1);
    const WeightWSM rhs = get_product_or_throw(two_k_plus_one, final_index);
    if (lhs <= rhs) {
      weights[ii] = kk;
    } else {
      weights[ii] = kk + 1;
    }
  }
  return weights;
}

PatternGraphData::PatternGraphData(
    const std::vector<std::pair<VertexWSM, VertexWSM>>& gate_sequence,
    const PatternGraphDataInput& input) {
  if (input.method == PatternGraphDataInput::ReorderingMethod::ORIGINAL_ORDER) {
    // Just use the index in gate_sequence for the time.
    const std::vector<WeightWSM> weights = get_interpolated_weights(
        input.initial_gate_weight, input.final_gate_weight,
        gate_sequence.size());
    for (unsigned ii = 0; ii < gate_sequence.size(); ++ii) {
      // C++ standard: initialised to zero on first use.
      auto& weight_sum = pattern_graph_weights[get_edge(
          gate_sequence[ii].first, gate_sequence[ii].second)];
      weight_sum = get_sum_or_throw(weight_sum, weights[ii]);
    }
    final_time = gate_sequence.size() - 1;
    return;
  }

  // We must reorder with parallel time slices.
  reordered_gates.reserve(gate_sequence.size());

  // KEY: a vertex
  // VALUE: the time when an edge with that vertex was last seen.
  // We need this for an O(N.log N) algorithm; simply moving each gate
  // backwards in time naively checking gates one-by-one
  // would be O(N^2).
  std::map<VertexWSM, unsigned> most_recent_time_map;

  GateTiming single_gate_data;

  for (const auto& pair : gate_sequence) {
    single_gate_data.time = 0;
    {
      const std::optional<unsigned> time1_opt =
          get_optional_value(most_recent_time_map, pair.first);
      if (time1_opt) {
        single_gate_data.time = 1 + time1_opt.value();
      }
    }
    {
      const std::optional<unsigned> time2_opt =
          get_optional_value(most_recent_time_map, pair.second);
      if (time2_opt) {
        single_gate_data.time =
            std::max<unsigned>(single_gate_data.time, 1 + time2_opt.value());
      }
    }
    single_gate_data.gate = get_edge(pair.first, pair.second);
    reordered_gates.emplace_back(single_gate_data);
    most_recent_time_map[pair.first] = single_gate_data.time;
    most_recent_time_map[pair.second] = single_gate_data.time;
  }
  final_time = 0;
  for (const auto& entry : most_recent_time_map) {
    final_time = std::max(final_time, entry.second);
  }
  const auto weights = get_interpolated_weights(
      input.initial_gate_weight, input.final_gate_weight, final_time + 1);
  for (const GateTiming& entry : reordered_gates) {
    auto& weight_sum = pattern_graph_weights[entry.gate];
    weight_sum = get_sum_or_throw(weight_sum, weights[entry.time]);
  }
}

TargetGraphData::TargetGraphData(DebugNoInputData) {}

void TargetGraphDataInput::check_validity() const {
  if (new_weight_multiplier <= 2 || max_number_of_new_edge_generations < 1 ||
      max_weight_multiplier < 2 || max_edge_density_percentage == 0 ||
      max_largest_to_smallest_final_weight_ratio < 2) {
    std::stringstream ss;
    ss << "TargetGraphData::Input: invalid values: new_weight_multiplier="
       << new_weight_multiplier << ", max_number_of_new_edge_generations="
       << max_number_of_new_edge_generations
       << ", max_edge_density_percentage=" << max_edge_density_percentage
       << ", max_largest_to_smallest_final_weight_ratio="
       << max_largest_to_smallest_final_weight_ratio;
    throw std::runtime_error(ss.str());
  }
}

// Go through all "V shapes"  x--y--z   and calculate the weight of
// a new edge x--z, storing it within "next_weights_to_add"
static void get_next_edge_weights_to_add(
    const GraphEdgeWeights& current_edges_and_weights,
    WeightWSM largest_allowed_weight, const std::vector<VertexWSM>& vertices,
    GraphEdgeWeights& next_weights_to_add, WeightWSM weight_multiplier) {
  next_weights_to_add.clear();
  for (VertexWSM root_v : vertices) {
    // Process all edges (v,x). Recall that every edge occurs twice,
    // so we definitely can start with (v,x) for smallest x.
    auto citer = current_edges_and_weights.lower_bound(
        std::make_pair(root_v, VertexWSM(0)));
    TKET_ASSERT(citer != current_edges_and_weights.cend());
    const EdgeWSM& edge = citer->first;
    TKET_ASSERT(edge.first == root_v);
    // Now, just consider all pairs (x,z) where (v,x) and (v,z) both exist.
    // So, we go through edges (v,x1), (v,x2), ..., (v,z).
    for (auto citer_x = citer; citer_x != current_edges_and_weights.cend();
         ++citer_x) {
      const EdgeWSM& edge_x = citer_x->first;
      if (edge_x.first != root_v) {
        break;
      }
      const VertexWSM& xx = edge_x.second;
      const WeightWSM& x_weight = citer_x->second;
      if (x_weight >= largest_allowed_weight) {
        // Only include initial excess weights (i.e., the first time
        // an added edge goes over the limit and is capped).
        // We don't want to combine an already-capped weight again,
        // as that's even worse.
        continue;
      }
      auto citer_z = citer_x;
      for (++citer_z; citer_z != current_edges_and_weights.cend(); ++citer_z) {
        const EdgeWSM& edge_z = citer_z->first;
        if (edge_z.first != root_v) {
          break;
        }
        const VertexWSM& zz = edge_z.second;
        const WeightWSM& z_weight = citer_z->second;
        TKET_ASSERT(xx < zz);
        if (z_weight >= largest_allowed_weight) {
          continue;
        }

        // Now we combine the weights of the 2 edges x--y, y--z to get
        // a new weight for the (possibly new) edge x--z.
        WeightWSM new_weight = get_product_or_throw(
            weight_multiplier, get_sum_or_throw(x_weight, z_weight));
        new_weight = std::min(new_weight, largest_allowed_weight);

        // It might not actually be new, of course!
        const EdgeWSM new_edge = get_edge(xx, zz);
        {
          const auto existing_weight_opt =
              get_optional_value(current_edges_and_weights, new_edge);
          if (existing_weight_opt) {
            new_weight = std::min(new_weight, existing_weight_opt.value());
          }
        }
        {
          const auto existing_weight_opt =
              get_optional_value(next_weights_to_add, new_edge);
          if (existing_weight_opt) {
            new_weight = std::min(new_weight, existing_weight_opt.value());
          }
        }
        next_weights_to_add[new_edge] = new_weight;
      }
    }
  }
}

TargetGraphData::TargetGraphData(
    GraphEdgeWeights original_target_weights,
    const TargetGraphDataInput& input) {
  try {
    input.check_validity();
    if (original_target_weights.empty()) {
      throw std::runtime_error("no input target edges!");
    }
    GetVerticesOptions options;
    options.allow_duplicate_edges = true;
    options.allow_edge_vertices_not_in_order = true;
    options.allow_zero_weights = false;

    // This checks the edges and weights for validity.
    sorted_vertices = get_vertices(original_target_weights, options);
    const unsigned number_of_vertices = sorted_vertices.size();

    WeightWSM largest_allowed_weight;
    {
      WeightWSM smallest_weight;
      set_maximum(smallest_weight);
      for (const auto& entry : original_target_weights) {
        smallest_weight = std::min(smallest_weight, entry.second);
      }

      largest_allowed_weight = get_product_or_throw<WeightWSM>(
          smallest_weight, input.max_largest_to_smallest_final_weight_ratio);
      TKET_ASSERT(smallest_weight > 0);
      TKET_ASSERT(largest_allowed_weight > smallest_weight);
    }

    GraphEdgeWeights new_weights;
    GraphEdgeWeights next_weights_to_add = std::move(original_target_weights);

    const auto add_reversed_edges = [&next_weights_to_add, &new_weights]() {
      for (const auto& entry : next_weights_to_add) {
        const EdgeWSM& edge = entry.first;
        const WeightWSM& weight = entry.second;
        new_weights[edge] = weight;
        new_weights[std::make_pair(edge.second, edge.first)] = weight;
      }
    };
    add_reversed_edges();

    // Double the number of edges in a complete graph,
    // as each edge is counted twice.
    const unsigned max_possible_map_size =
        number_of_vertices * (number_of_vertices + 1);
    unsigned max_weight_map_size = max_possible_map_size;
    unsigned number_of_generations;
    if (number_of_vertices <
        input.min_num_vertices_to_break_off_new_generations) {
      set_maximum(number_of_generations);
      set_maximum(max_weight_map_size);
    } else {
      number_of_generations = input.max_number_of_new_edge_generations;
      max_weight_map_size =
          (max_weight_map_size * input.max_edge_density_percentage) / 100;
    }

    // Now, repeatedly add new edges and weights.
    for (unsigned generation = 0; generation < number_of_generations;
         ++generation) {
      const auto current_size = new_weights.size();
      if (new_weights.size() >= max_weight_map_size) {
        break;
      }
      get_next_edge_weights_to_add(
          new_weights, largest_allowed_weight, sorted_vertices,
          next_weights_to_add, input.new_weight_multiplier);

      add_reversed_edges();
      TKET_ASSERT(current_size <= new_weights.size());
      if (current_size == new_weights.size()) {
        break;
      }
    }
    TKET_ASSERT(new_weights.size() % 2 == 0);
    TKET_ASSERT(new_weights.size() <= max_possible_map_size);

    // Calculate the implicit weight of all remaining unmentioned edges.
    // Also, remove all
    implicit_weight = 0;
    for (const auto& entry : new_weights) {
      const EdgeWSM& edge = entry.first;
      const WeightWSM& weight = entry.second;
      implicit_weight = std::max(implicit_weight, weight);
      if (edge.first < edge.second) {
        explicit_target_graph_weights.emplace(entry);
      }
    }
    implicit_weight = get_product_or_throw<WeightWSM>(
        implicit_weight, input.max_weight_multiplier);
  } catch (const std::exception& e) {
    std::stringstream ss;
    ss << "WSM initial qubit placement: constructing target: " << e.what();
    throw std::runtime_error(ss.str());
  }
}

WeightWSM TargetGraphData::get_edge_weight(VertexWSM tv1, VertexWSM tv2) const {
  const auto edge = get_edge(tv1, tv2);
  const auto weight_opt =
      get_optional_value(explicit_target_graph_weights, edge);
  if (weight_opt) {
    return weight_opt.value();
  }
  return implicit_weight;
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
