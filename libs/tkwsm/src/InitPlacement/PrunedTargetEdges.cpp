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

#include "tkwsm/InitPlacement/PrunedTargetEdges.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

// Given all assignments PV->TV, keep exactly those target edges
// which need to exist for a valid monomorphism (using the implicit
// edge weight if th eedge is not explicitly given),
// discarding all other target edges.
static GraphEdgeWeights get_only_used_target_edges(
    const NeighboursData& pattern_ndata,
    const NeighboursData& explicit_target_ndata,
    WeightWSM implicit_target_weight,
    const std::vector<unsigned>& assigned_target_vertices) {
  GraphEdgeWeights new_target_graph_data;
  TKET_ASSERT(!assigned_target_vertices.empty());
  TKET_ASSERT(
      assigned_target_vertices.size() ==
      pattern_ndata.get_number_of_nonisolated_vertices());
  for (unsigned pv1 = 0; pv1 < assigned_target_vertices.size(); ++pv1) {
    const unsigned tv1 = assigned_target_vertices[pv1];
    const auto& neighbours_and_weights =
        pattern_ndata.get_neighbours_and_weights(pv1);
    // Only consider p-edges (pv1, pv2) with pv1 < pv2.
    for (auto citer = std::lower_bound(
             neighbours_and_weights.cbegin(), neighbours_and_weights.cend(),
             std::make_pair(VertexWSM(pv1 + 1), WeightWSM(0)));
         citer != neighbours_and_weights.cend(); ++citer) {
      // We have an edge (pv1, pv2).
      const unsigned tv2 = assigned_target_vertices.at(citer->first);
      const auto t_edge = get_edge(tv1, tv2);
      TKET_ASSERT(new_target_graph_data.count(t_edge) == 0);

      // We NEED a corresponding target edge to exist in the new graph.
      const auto existing_weight_opt =
          explicit_target_ndata.get_edge_weight_opt(tv1, tv2);
      if (existing_weight_opt) {
        new_target_graph_data[t_edge] = existing_weight_opt.value();
      } else {
        new_target_graph_data[t_edge] = implicit_target_weight;
      }
    }
  }
  return new_target_graph_data;
}

static unsigned get_max_number_of_new_target_edges(
    const NeighboursData& explicit_target_ndata,
    unsigned number_of_current_used_target_edges,
    const TargetEdgePruningParameters& parameters) {
  const unsigned total_number_of_possible_tv =
      explicit_target_ndata.get_number_of_nonisolated_vertices();
  const unsigned number_of_complete_target_edges =
      (total_number_of_possible_tv * (total_number_of_possible_tv - 1)) / 2;
  TKET_ASSERT(
      number_of_complete_target_edges >= number_of_current_used_target_edges);
  if (number_of_complete_target_edges == number_of_current_used_target_edges) {
    return 0;
  }

  unsigned max_target_edges =
      number_of_current_used_target_edges +
      get_product_or_throw(
          parameters.max_additional_number_of_target_edges_factor_per_kilo,
          number_of_current_used_target_edges) /
          1024;

  const unsigned min_unused_target_edges =
      (get_product_or_throw(
           parameters
               .min_implicit_unused_number_of_target_edges_factor_per_kilo,
           number_of_complete_target_edges) -
       number_of_current_used_target_edges) /
      1024;

  const unsigned max_edges_from_this_constraint =
      number_of_complete_target_edges - min_unused_target_edges;
  max_target_edges = std::min(max_target_edges, max_edges_from_this_constraint);
  if (max_target_edges <= number_of_current_used_target_edges) {
    return 0;
  }
  return max_target_edges - number_of_current_used_target_edges;
}

// A list of (weight, t-edge) pairs; will be sorted by increasing weight,
// since we prefer to disallow larger target weights if we can.
typedef std::vector<std::pair<WeightWSM, std::pair<VertexWSM, VertexWSM>>>
    WeightAndTEdgeList;

// We'll prioritise new t-edges between TWO used TVs over those
// touching only ONE used TV, and also lower weights over highrer weights.

// Get all extra t-edges, both explcit and implicit, between TWO used TV.
static void get_extra_edges_adjoining_two_used_tv_unsorted(
    const NeighboursData& explicit_target_ndata,
    const std::vector<unsigned>& original_used_tv_sorted,
    const GraphEdgeWeights& new_target_graph_data,
    WeightAndTEdgeList& explicit_t_edges_data_joining_two_tv,
    std::vector<std::pair<VertexWSM, VertexWSM>>&
        implicit_t_edges_joining_two_tv) {
  for (auto citer_tv1 = original_used_tv_sorted.cbegin();
       citer_tv1 != original_used_tv_sorted.cend(); ++citer_tv1) {
    // Add edges (tv1, tv2) where tv2 > tv1 and tv1, tv2 are BOTH
    // used vertices, whether or not the edge is explicit.
    for (auto citer_tv2 = citer_tv1 + 1;
         citer_tv2 != original_used_tv_sorted.cend(); ++citer_tv2) {
      const EdgeWSM t_edge = get_edge(*citer_tv1, *citer_tv2);
      if (new_target_graph_data.count(t_edge) != 0) {
        continue;
      }
      const auto explicit_weight_opt =
          explicit_target_ndata.get_edge_weight_opt(*citer_tv1, *citer_tv2);
      if (explicit_weight_opt) {
        explicit_t_edges_data_joining_two_tv.emplace_back(
            explicit_weight_opt.value(), t_edge);
      } else {
        implicit_t_edges_joining_two_tv.emplace_back(t_edge);
      }
    }
  }
}

// Get all extra explcit t-edges adjoining exactly ONE used TV.
static WeightAndTEdgeList
get_explicit_extra_edges_adjoining_one_used_tv_unsorted(
    const NeighboursData& explicit_target_ndata,
    const std::vector<unsigned>& original_used_tv_sorted) {
  WeightAndTEdgeList t_edges_joining_one_tv;
  for (auto citer_tv1 = original_used_tv_sorted.cbegin();
       citer_tv1 != original_used_tv_sorted.cend(); ++citer_tv1) {
    // Now, for TV2 we restrict to explicit neighbours of TV1.
    // Also, can just take TV2 > TV1.
    const auto& explicit_neighbours_and_weights =
        explicit_target_ndata.get_neighbours_and_weights(*citer_tv1);

    for (auto citer_other = std::lower_bound(
             explicit_neighbours_and_weights.cbegin(),
             explicit_neighbours_and_weights.cend(),
             std::make_pair(VertexWSM(*citer_tv1 + 1), WeightWSM(0)));
         citer_other != explicit_neighbours_and_weights.cend(); ++citer_other) {
      const auto& tv2 = citer_other->first;
      const auto& explicit_weight = citer_other->second;

      // Now, (tv1, tv2) is an EXPLICIT edge. Has it already been considered?
      // We can save a little time by not searching the whole range.
      if (!std::binary_search(
              citer_tv1 + 1, original_used_tv_sorted.cend(), tv2)) {
        // tv2 is unused, as we want.
        // Thus (tv1, tv2) automatically is unused so far,
        // because used edges have BOTH endpoints used (by definition!)
        t_edges_joining_one_tv.emplace_back(
            explicit_weight, get_edge(*citer_tv1, tv2));
      }
    }
  }
  return t_edges_joining_one_tv;
}

GraphEdgeWeights get_new_target_graph_data(
    const NeighboursData& pattern_ndata,
    const NeighboursData& explicit_target_ndata,
    WeightWSM implicit_target_weight,
    const std::vector<unsigned>& assigned_target_vertices,
    const TargetEdgePruningParameters& parameters) {
  GraphEdgeWeights new_target_graph_data{get_only_used_target_edges(
      pattern_ndata, explicit_target_ndata, implicit_target_weight,
      assigned_target_vertices)};

  // How many target edges do we want at the end?
  unsigned max_number_of_new_target_edges = get_max_number_of_new_target_edges(
      explicit_target_ndata, new_target_graph_data.size(), parameters);

  if (max_number_of_new_target_edges == 0) {
    return new_target_graph_data;
  }

  // Which TV are used in the ORIGINAL solution?
  std::vector<unsigned> original_used_tv = assigned_target_vertices;
  std::sort(original_used_tv.begin(), original_used_tv.end());
  TKET_ASSERT(
      std::adjacent_find(original_used_tv.cbegin(), original_used_tv.cend()) ==
      original_used_tv.cend());
  TKET_ASSERT(
      original_used_tv.back() <
      explicit_target_ndata.get_number_of_nonisolated_vertices());

  // Now, get data for additional edges to add.
  WeightAndTEdgeList explicit_t_edges_data_joining_two_tv;
  std::vector<std::pair<VertexWSM, VertexWSM>> implicit_t_edges_joining_two_tv;

  get_extra_edges_adjoining_two_used_tv_unsorted(
      explicit_target_ndata, original_used_tv, new_target_graph_data,
      explicit_t_edges_data_joining_two_tv, implicit_t_edges_joining_two_tv);

  // Add the explicit two vertex edges first, in order of increasing weight.
  std::sort(
      explicit_t_edges_data_joining_two_tv.begin(),
      explicit_t_edges_data_joining_two_tv.end());

  for (const auto& entry : explicit_t_edges_data_joining_two_tv) {
    new_target_graph_data[entry.second] = entry.first;
    --max_number_of_new_target_edges;
    if (max_number_of_new_target_edges == 0) {
      return new_target_graph_data;
    }
  }

  // Now, the explicit t-edges with only one used vertex.
  WeightAndTEdgeList t_edges_joining_one_tv =
      get_explicit_extra_edges_adjoining_one_used_tv_unsorted(
          explicit_target_ndata, original_used_tv);
  std::sort(t_edges_joining_one_tv.begin(), t_edges_joining_one_tv.end());
  for (const auto& entry : t_edges_joining_one_tv) {
    new_target_graph_data[entry.second] = entry.first;
    --max_number_of_new_target_edges;
    if (max_number_of_new_target_edges == 0) {
      return new_target_graph_data;
    }
  }

  // Finally, the implicit edges between two used TV.
  for (const auto& t_edge : implicit_t_edges_joining_two_tv) {
    new_target_graph_data[t_edge] = implicit_target_weight;
    --max_number_of_new_target_edges;
    if (max_number_of_new_target_edges == 0) {
      break;
    }
  }
  return new_target_graph_data;
}

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
