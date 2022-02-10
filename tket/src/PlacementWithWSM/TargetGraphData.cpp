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

#include "PlacementWithWSM/TargetGraphData.hpp"

#include <algorithm>
#include <stdexcept>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Common/SpecialExceptions.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {

using namespace WeightedSubgraphMonomorphism;

// If we have a path [v(1), v(2), v(3), ..., v(n)], with initial edge weights
// w(i), what should the cost of the new edge v(1) -- v(n) be?
//
// What would actually happen if qubits at v(1), v(n) needed a 2-qubit gate
// to be applied between them? Our model is that SWAP gates are inserted along
// the path to make the qubits adjacent, the gate is performed,
// and then the same SWAP gates are performed to move the qubits back again
// to the endpoints.
//
// In this model, if K is the number of primitive gates needed to make a SWAP
// (and all primitive 2-qubit gates along an edge are assumed to have
// the same fidelity, and no preferred edge direction),
// the cost is
//
//     2K . sum { w(i) : SWAP gates are performed }  +  w(j),
//
// where w(j) is the weight of the edge chosen to perform the non-SWAP gate.
// This is 2K.{sum w(i)} - (2K-1)w(j).
// Obviously we'd choose j to minimise this cost, so the answer is
//    2K.{sum w(i)} - (2K-1).(max w).
// (It's ironic that we perform the gate along the WORST edge
// in the path; we only use it once, unlike the others
// which are used multiple times for SWAPs).
//
// Of course, in practice we probably wouldn't do this; we might move the qubits
// to be adjacent, then leave them in place.
// (This is why we want a time decay factor in the pattern graph weights;
// as time goes on, it's less likely that qubits will be where they were
// initially, and the qubits will gradually drift away from their initial
// positions).
//
// So, the cost would instead be K.{sum w(i)} - (K-1).(max w).
// For K>1 this is a slightly odd cost function: if we choose paths [a,b], [b,c]
// giving the lowest costs and consider the concatenated path,
// we get a reverse triangle inequality:
//    (Cost along this path)(a,c) >= Cost(a,b) + Cost(b,c).
// So all the usual stuff with Dijkstra etc. doesn't apply;
// subpaths of optimal paths need not be optimal.
//
// We could also have an existing edge (v1,v2) being worse
// than taking a roundabout path [v1,a,b,c,...,v2].
//
// We'll do a little hack: to find the new costs from a source to many targets,
// we'll allow 2 hits on a target vertex (not just one) to get the cost.
// This is of course inaccurate, but hopefully accurate enough to be useful
// (this whole framework with weights etc. is only approximate anyway).
namespace {
// We've begun a search from a given source vertex
// to hit all target vertices; record our current status along the path.
struct PathData {
  unsigned path_length;
  WeightWSM sum_of_weights_so_far;
  WeightWSM max_weight_so_far;
  VertexWSM end_vertex;
};
}  // namespace

static void check_for_overflow(
    const TargetGraphData& tgraph_data,
    const TargetGraphData::Parameters& parameters) {
  WeightWSM total_weight = 0;
  for (const auto& entry : tgraph_data.final_data) {
    total_weight = get_sum_or_throw(total_weight, entry.second);
  }
  total_weight =
      get_product_or_throw(total_weight, WeightWSM(parameters.swap_gate_count));
  total_weight = get_product_or_throw(total_weight, WeightWSM(10));
}

static void erase_high_weights(GraphEdgeWeights& data, WeightWSM max_weight) {
  std::vector<EdgeWSM> edges_to_erase;
  for (const auto& entry : data) {
    if (entry.second > max_weight) {
      edges_to_erase.push_back(entry.first);
    }
  }
  for (const auto& edge : edges_to_erase) {
    data.erase(edge);
  }
}

static void cap_high_weights(GraphEdgeWeights& data, WeightWSM max_weight) {
  for (auto& entry : data) {
    entry.second = std::min(entry.second, max_weight);
  }
}

static void continue_setup(
    TargetGraphData& tgraph_data, const TargetGraphData::Parameters& parameters,
    const NeighboursData& initial_ndata) {
  TKET_ASSERT(parameters.max_edge_weight);
  check_for_overflow(tgraph_data, parameters);

  std::vector<PathData> all_paths;
  for (auto source_v : tgraph_data.sorted_vertices) {
    all_paths.resize(1);
    all_paths[0].path_length = 0;
    all_paths[0].sum_of_weights_so_far = 0;
    all_paths[0].max_weight_so_far = 0;
    all_paths[0].end_vertex = source_v;

    while (!all_paths.empty()) {
      const auto path = all_paths.back();
      all_paths.pop_back();

      bool should_extend =
          path.path_length < parameters.max_path_length_for_new_edges;

      if (path.path_length > 0) {
        bool should_add_data = true;
        const auto edge = get_edge(source_v, path.end_vertex);
        const bool already_hit_vertex = tgraph_data.final_data.count(edge) != 0;

        const bool edge_existed_initially =
            initial_ndata.get_edge_weight_opt(source_v, path.end_vertex)
                .has_value();

        if (edge_existed_initially) {
          // We're at a neighbour of the source vertex.
          if (path.path_length > 1) {
            // It's the second time we've hit the neighbour.
            should_extend = false;
            should_add_data =
                parameters
                    .replace_low_fidelity_primitive_gates_with_longer_paths;
          }
        } else {
          if (path.end_vertex == source_v) {
            // We've looped round...
            should_extend = false;
            should_add_data = false;
          } else {
            // We've hit a new vertex, not an original neighbour.
            if (already_hit_vertex) {
              should_extend = false;
            }
          }
        }
        // NOW, what do we do?
        if (should_add_data) {
          const WeightWSM new_weight =
              parameters.swap_gate_count * path.sum_of_weights_so_far -
              (parameters.swap_gate_count - 1) * path.max_weight_so_far;
          auto& weight_entry = tgraph_data.final_data[edge];
          if (already_hit_vertex) {
            weight_entry = std::min(weight_entry, new_weight);
          } else {
            weight_entry = new_weight;
          }
        }
      }
      if (should_extend) {
        const auto& new_neighbours_data =
            initial_ndata.get_neighbours_and_weights(path.end_vertex);
        for (const auto& entry : new_neighbours_data) {
          const auto& new_v = entry.first;
          if (new_v == source_v) {
            continue;
          }
          const auto& new_w = entry.second;
          all_paths.emplace_back();
          all_paths.back().end_vertex = new_v;
          all_paths.back().max_weight_so_far =
              std::max(path.max_weight_so_far, new_w);
          all_paths.back().path_length = path.path_length + 1;
          all_paths.back().sum_of_weights_so_far =
              path.sum_of_weights_so_far + new_w;
        }
      }
    }
  }
  // Now, deal with weights over the limit.
  if (parameters.remove_high_edge_weights) {
    erase_high_weights(
        tgraph_data.final_data, parameters.max_edge_weight.value());
  } else {
    cap_high_weights(
        tgraph_data.final_data, parameters.max_edge_weight.value());
  }
}

TargetGraphData::TargetGraphData(GraphEdgeWeights data, Parameters parameters)
    : final_data(std::move(data)) {
  if (!(!final_data.empty() && parameters.swap_gate_count > 0 &&
        parameters.max_path_length_for_new_edges > 0 &&
        parameters.max_edge_weight_largest_weight_ratio > 5 &&
        parameters.max_edge_weight_smallest_weight_ratio > 5)) {
    throw std::runtime_error("Invalid inputs for TargetGraphData");
  }
  WeightWSM max_weight = 0;
  WeightWSM min_weight;
  set_maximum(min_weight);
  for (const auto& entry : final_data) {
    max_weight = std::max(max_weight, entry.second);
    if (entry.second > 0) {
      min_weight = std::min(min_weight, entry.second);
    }
  }
  TKET_ASSERT(!is_maximum(min_weight));
  if (max_weight == 0) {
    throw std::runtime_error("Invalid max weight 0 for TargetGraphData");
  }
  WeightWSM max_weight_to_use;
  set_maximum(max_weight_to_use);
  if (parameters.max_edge_weight) {
    if (parameters.max_edge_weight.value() <= max_weight) {
      throw std::runtime_error(
          "Invalid parameters.max_edge_weight.value() for TargetGraphData");
    }
    max_weight_to_use =
        std::min(max_weight_to_use, parameters.max_edge_weight.value());
  }
  {
    const auto weight_prod_opt = get_checked_product(
        max_weight, parameters.max_edge_weight_largest_weight_ratio);
    if (weight_prod_opt) {
      max_weight_to_use = std::min(max_weight_to_use, weight_prod_opt.value());
    }
  }
  {
    const auto weight_prod_opt = get_checked_product(
        min_weight, parameters.max_edge_weight_smallest_weight_ratio);
    if (weight_prod_opt) {
      max_weight_to_use = std::min(max_weight_to_use, weight_prod_opt.value());
    }
  }
  const NeighboursData initial_ndata(final_data);
  sorted_vertices = initial_ndata.get_nonisolated_vertices_expensive();
  parameters.max_edge_weight = max_weight_to_use;
  continue_setup(*this, parameters, initial_ndata);
}

}  // namespace tket
