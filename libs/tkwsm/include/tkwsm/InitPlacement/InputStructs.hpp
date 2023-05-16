// Copyright 2019-2023 Cambridge Quantum Computing
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
#include <optional>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {

// Handy for testing, to allow constructing without input data.
struct DebugNoInputData {};

struct PatternGraphDataInput {
  // The initial gate weight should be higher than the final gate weight.
  // Weights for other gates will be found by interpolation,
  // with integer-only operations which are thus fully portable.
  // The values should be large enough that we get a good spread
  // of discrate values, but not so large that integer overflow is likely.
  WeightWSM initial_gate_weight = 100;
  WeightWSM final_gate_weight = 20;

  enum class ReorderingMethod { TIME_SLICES_OF_PARALLEL_GATES, ORIGINAL_ORDER };
  ReorderingMethod method = ReorderingMethod::TIME_SLICES_OF_PARALLEL_GATES;
};

/** Convert a sequence of gates into a weighted pattern graph,
 * to pass into WSM.
 * Note that this uses input vertex numbers; no relabelling.
 */
struct PatternGraphData {
  struct GateTiming {
    // Uses original input vertices, no relabelling.
    EdgeWSM gate;
    unsigned time;
  };

  // The gates may be reordered according to time slices.
  // If so, record them here. If not (i.e., we just leave the original
  // input gate sequence), leave this empty.
  std::vector<GateTiming> reordered_gates;

  GraphEdgeWeights pattern_graph_weights;

  // The maximum time occurring in "reordered_gates"
  unsigned final_time = 0;

  explicit PatternGraphData(
      const std::vector<std::pair<VertexWSM, VertexWSM>>& gate_sequence,
      const PatternGraphDataInput& input = {});

  /** For testing, just allow construction without any data.
   * @param dummy_object A dummy class object, used merely to indicate that
   * we're constructing the object without input data (i.e., not filling
   * anything).
   */
  explicit PatternGraphData(DebugNoInputData dummy_object);
};

struct TargetGraphDataInput {
  // How many primitive 2-qubit gates make up a SWAP gate?
  // unsigned number_of_gates_for_swap = 3;

  // When we have a path x-y-z, with no edge currently between x,z,
  // the new edge x-z will be created with weight
  // M * (weight(x,y) + weight(y,z)) for some M.
  // This is the value of M.
  unsigned new_weight_multiplier = 3;

  // If we don't have very many target vertices,
  // just fill in all the edges without worry, without using
  // the below parameters which restrict the number of new edges
  // explicitly added.
  unsigned min_num_vertices_to_break_off_new_generations = 10;

  // Only takes effect if we have enough vertices.
  // (I.e., if we have fewer than
  // min_num_vertices_to_break_off_new_generations vertices,
  // then it behaves as if this parameter is set to +infinity).
  // Each new generation takes time O(V.d^2) to calculate, where V is the
  // number of vertices, and d is the maximum degree. So for n generations,
  // time is O(V.(d^2 + (d+1)^2 + .. + (d+n-1)^2)
  //    = O(V.max(d,n)^3).
  // Obviously this is scary if V=1000 and d,n~V.
  // In future, if it's too slow for large graphs,
  // one could try a Monte Carlo approach, where we just add
  // a limited number of random edges.
  unsigned max_number_of_new_edge_generations = 5;

  // Only takes effect if we have enough vertices.
  unsigned max_edge_density_percentage = 25;

  // Only takes effect if we have enough vertices.
  // Once we've computed the final explicit edge weight,
  // any other edges are implicitly assumed to have k*(max weight).
  // This is the value of k.
  unsigned max_weight_multiplier = 3;

  // Once a weight gets too large, cap it at the maximum value
  unsigned max_largest_to_smallest_final_weight_ratio = 50;

  // Throws if not valid.
  void check_validity() const;
};

/** Convert a given architecture (physical qubit layout) with weighted edges
 * (weight being proportional to error probability of a primitive gate
 * applied along that edge) to a COMPLETE target weighted graph
 * (although some edges are left unmentioned, so regarded as "implicitly"
 * existing and having an edge weight).
 * This then gives a WSM problem (although actually NOT solved with
 * the main WSM routines, but rather with MCCT - Monte Carlo Complete Target -
 * since it's more similar to Travelling Salesperson problems, requiring random
 * jumping etc. similar to Simulated Annealing, rather than usual
 * WSM problems, where graph theory and local vertex properties etc.
 * is more important).
 */
struct TargetGraphData {
  // Experimentation is advisable, but the default values
  // should be sensible.

  // Any target edge not explicitly listed will have this weight
  // (we don't want a dense representation of a complete target graph
  // with 1000 vertices!)
  WeightWSM implicit_weight = 0;

  // The list of all target vertices, sorted.
  std::vector<VertexWSM> sorted_vertices;

  // Again, the target graph will be COMPLETE; any edges not listed in here
  // will be treated as though they have "implicit_weight".
  // This uses the original vertices; no relabelling.
  // The keys will be edges (v1,v2) with v1<v2.
  GraphEdgeWeights explicit_target_graph_weights;

  explicit TargetGraphData(
      GraphEdgeWeights original_target_weights,
      const TargetGraphDataInput& input = {});

  explicit TargetGraphData(DebugNoInputData dummy_object);

  WeightWSM get_edge_weight(VertexWSM tv1, VertexWSM tv2) const;
};

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
